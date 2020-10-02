import argparse
import strformat
import os
import d4
import hts
import tables
import ./utils
import ./typeenum
import algorithm
import times


proc summarize_block(d4s: var seq[Cover], chrom:string, start:uint32, stop:uint32, percentiles: seq[float], outs:var seq[D4], vals: var seq[seq[int32]], out_vals: var seq[seq[int32]]) =
  # extract data for each sample from start to stop, calculate per-base
  # percentiles and write to outs.
  # vals and out_vals just let us avoid allocations
  for i in 0..vals.high:
    if vals[i].len != int(stop) - int(start):
      vals[i].setLen(int(stop) - int(start))
      echo "changing length:", vals[i].len, " start:", start, " stop:", stop

    d4s[i].values(chrom, start, vals[i])

  if out_vals[0].len != vals[0].len:
    for i in 0..<out_vals.len:
      out_vals[i].setLen(vals[0].len)

  let n_samples = d4s.len.float

  # we just keep here which sample we'll stop at when getting to each
  # percentile
  var int_sample = newSeqUninitialized[uint32](percentiles.len)
  for i, p in percentiles:
    int_sample[i] = uint32(max(0, int(n_samples * p - 0.5)))

  var H : array[1024, uint32]
  for i in 0..vals[0].high:
    zeroMem(H[0].addr.pointer, sizeof(H[0]) * H.len)
    # fill the histogram
    for j in 0..d4s.high:
      H[min(H.high, vals[j][i])].inc

    var cutoff_index = 0
    var isamp = int_sample[cutoff_index]

    # extract the percentiles from the histogram.
    # this is faster than doing binary search for each value as we iterate over
    # array only once and most values are near the left end.
    var high_set = false
    block cumsum_block:
      var cumsum:uint32
      for k, n in H:
        cumsum += n
        # use while in case we had a chunk of samples that span some percentiles.
        while cumsum > isamp:
          out_vals[cutoff_index][i] = k
          cutoff_index += 1
          if cutoff_index > int_sample.high:
            high_set = true
            break cumsum_block
          isamp = int_sample[cutoff_index]
    # if we ran out of samples, just set it to highest value.
    if not high_set:
      out_vals[^1][i] = H.len.int32

  for i in 0..out_vals.high:
    outs[i].write(chrom, start, out_vals[i])


proc generate_backgrounds(dps: var seq[Cover], output_dir: string, percentiles: seq[int], fai:Fai) =
  let block_size = 1_000_000'u32
  let info = dps[0].chromosomes(fai)
  var chroms = newSeqOfCap[tuple[name:string, length:int]](info.len)
  for k, v in info: chroms.add((k, v.int))

  discard outputDir.existsOrCreateDir
  var percentiles = sorted(percentiles)

  var pctiles = newSeq[float](percentiles.len)
  for i, p in percentiles:
    pctiles[i] = p.float / 100.0

  var outs = newSeq[D4](percentiles.len)
  for i in 0..<percentiles.len:
    doAssert outs[i].open(&"{outputDir}/seqcover_p{percentiles[i]}.d4", mode="w")
    outs[i].set_chromosomes(chroms)
  var vals = newSeq[seq[int32]](dps.len)
  for v in vals.mitems: v.setLen(block_size)

  var out_vals = newSeq[seq[int32]](outs.len)
  for i in 0..outs.high:
    out_vals[i] = newSeqUninitialized[int32](vals[0].len)


  for (chrom, length) in info.pairs:
    dps.check_same_lengths(chrom, length, fai)
    let t0 = cpuTime()
    var j = 0
    for i in countup(0'u32, length, block_size):
      if j mod 10 == 0 and i > 0:
        echo &"{chrom}:{min(i + blocksize, length)} time:{cpuTime() - t0:.2f} bases/second: {min(i + blocksize, length).float / (cpuTime() - t0):.0f}"
      dps.summarize_block(chrom, i, min(i + block_size, length), pctiles, outs, vals, out_vals)
      j.inc

  for o in outs.mitems:
    o.close()

proc generate_background_main*() =
  let p = newParser("seqcover generate-background"):
    option("-o", "--output-dir", help="directory for output", default="seqcover-backgrounds")
    option("-f", "--fasta", help="indexed fasta required for bed.gz files")
    option("-p", "--percentile", default="5", help="percentile to extract from backgrounds")
    arg("samples", nargs= -1, help="per-base bed.gz files or d4 files or a glob of either to generate background")

  var argv = commandLineParams()
  if len(argv) > 0 and argv[0] == "generate-background":
    argv = argv[1..argv.high]
  if len(argv) == 0:
    argv.add("--help")

  var fai:Fai

  var opts = p.parse(argv)
  if opts.help:
    quit 0

  let percentiles = @[parseInt(opts.percentile)]

  if opts.fasta != "":
    doAssert fai.open(opts.fasta), "[seqcover] error opening fasta file:" & opts.fasta

  let paths = get_glob_samples(opts.samples)
  echo paths
  if paths.len < 4:
    quit "[seqcover] need at least 4 samples to generate a background"
  if paths.len < 10:
    stderr.write_line "[seqcover] warning: creating background with fewer than 10 samples might give unexpected results"

  var covers = read_dps(paths)
  covers.generate_backgrounds(opts.output_dir, percentiles, fai)

  echo opts


when isMainModule:
  generate_background_main()
