import argparse
import strformat
import os
import d4
import tables
import ./utils
import algorithm
import times

proc summarize_block(d4s: var seq[D4], chrom:string, start:uint32, stop:uint32, percentiles: seq[float], outs:var seq[D4], vals: var seq[seq[int32]]) =
  # extract data for each sample from start to stop, calculate per-base
  # percentiles and write to outs.
  for i in 0..vals.high:
    if vals[i].len != int(stop) - int(start):
      vals[i].setLen(int(stop) - int(stop))
      echo "changing length:", vals[i].len, " start:", start, " stop:", stop

    d4s[i].values(chrom, start, vals[i])

  let n_samples = d4s.len.float

  # we just keep here which sample we'll stop at when getting to each
  # percentile
  var int_sample = newSeq[int](percentiles.len)
  for i, p in percentiles:
    int_sample[i] = max(0, int(n_samples * p - 0.5))

  var out_vals = newSeq[seq[int32]](outs.len)
  for i in 0..outs.high:
    out_vals[i] = newSeqUninitialized[int32](vals[0].len)

  var column = newSeq[int32](d4s.len)

  var H : array[256, uint32]
  for i in 0..vals[0].high:
    zeroMem(H[0].addr.pointer, sizeof(H[0]) * H.len)
    # fill the histogram
    for j in 0..d4s.high:
      H[min(H.high, vals[j][i])].inc

    for k, isamp in int_sample:
      var tot = 0'u32
      var m = 0
      for ns in H:
        tot += ns
        if tot > isamp.uint32:
          #echo "cutoff:", n, ".. ", m, " vs ", column[n], " tot:", tot
          break
        m.inc
      out_vals[k][i] = m.int32 # column[n]

  for i in 0..out_vals.high:
    outs[i].write(chrom, start, out_vals[i])


proc generate_backgrounds(d4s: var seq[D4], output_dir: string, percentiles: seq[int]) =
  let block_size = 1_000_000'u32
  let info = d4s[0].chromosomes
  var chroms = newSeq[tuple[name:string, length:int]](info.len)
  for k, v in info: chroms.add((k, v.int))

  discard outputDir.existsOrCreateDir
  var percentiles = sorted(percentiles)

  var pctiles = newSeq[float](percentiles.len)
  for i, p in percentiles:
    pctiles[i] = p.float / 100.0

  var outs = newSeq[D4](percentiles.len)
  for i in 0..<percentiles.len:
    doAssert outs[i].open(&"{outputDir}/seqcover_p{percentiles[i]}", mode="w")
    outs[i].set_chromosomes(chroms)
  var vals = newSeq[seq[int32]](d4s.len)
  for v in vals.mitems: v.setLen(block_size)

  for (chrom, length) in info.pairs:
    d4s.check_same_lengths(chrom, length)
    let t0 = cpuTime()
    var j = 0
    for i in countup(0'u32, length, block_size):
      d4s.summarize_block(chrom, i, min(i + block_size, length), pctiles, outs, vals)
      if j mod 5 == 0:
        echo &"{chrom}:{min(i + blocksize, length)} time:{cpuTime() - t0:.2f} bases/second: {min(i + blocksize, length).float / (cpuTime() - t0):.0f}"
      j.inc
      if i > 40_000_000: break
    break

  for o in outs.mitems:
    o.close()

proc main() =
  let p = newParser("seqcover background"):
    option("-o", "--output-dir", help="directory for output", default="d4-backgrounds")
    arg("samples", nargs= -1, help="d4 files or a glob of d4 files to generate background")


  let percentiles = @[10, 50, 90]
  var argv = commandLineParams()
  if len(argv) > 0 and argv[0] == "report":
    argv = argv[1..argv.high]
  if len(argv) == 0:
    argv.add("--help")

  var opts = p.parse(argv)
  if opts.help:
    quit 0

  let paths = get_glob_samples(opts.samples)
  echo paths
  if paths.len < 4:
    quit "[seqcover] need at least 4 samples to generate a background"
  if paths.len < 10:
    stderr.write_line "[seqcover] warning: creating background with fewer than 10 samples might give unexpected results"

  var d4s = read_d4s(paths)
  d4s.generate_backgrounds(opts.output_dir, percentiles)



  echo opts


when isMainModule:
  main()
