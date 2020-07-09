import argparse
import strformat
import os
import d4
import tables
import ./utils
import algorithm
import times

proc summarize_block(d4s: var seq[D4], chrom:string, start:uint32, stop:uint32, percentiles: seq[float], outs:var seq[D4]) =
  var vals = newSeq[seq[int32]](d4s.len)
  for i in 0..d4s.high:
    #d in d4s.mitems:
    vals.add(d4s[i].values(chrom, start, stop))

  # we just keep here which sample we'll stop at when getting to each
  # percentile
  let n_samples = d4s.len.float
  var ipd = newSeq[int](percentiles.len)
  for i, p in percentiles:
    ipd[i] = max(0, int(n_samples * p - 0.5))

  var column = newSeq[int32](d4s.len)
  var out_vals = newSeq[seq[int32]](outs.len)
  for i in 0..outs.high:
    out_vals[i] = newSeqUninitialized[int32](vals[0].len)

  for i in 0..vals[0].high:
    for j in 0..d4s.high:
      column[j] = vals[j][i]
    # TODO: use histogram...
    sort(column)
    for k, n in ipd:
      out_vals[k][i] = column[n]

  for i in 0..out_vals.high:
    outs[i].write(chrom, start, out_vals[i])


proc generate_backgrounds(d4s: var seq[D4], output_dir: string, percentiles: seq[int]) =
  let block_size = 1_000_000'u32
  let info = d4s[0].chromosomes
  var chroms = newSeq[tuple[name:string, length:int]](info.len)
  for k, v in info: chroms.add((k, v.int))

  discard outputDir.existsOrCreateDir

  var pctiles = newSeq[float](percentiles.len)
  for i, p in percentiles:
    pctiles[i] = p.float / 100.0

  var outs = newSeq[D4](percentiles.len)
  for i in 0..<percentiles.len:
    doAssert outs[i].open(&"{outputDir}/seqcover_p{percentiles[i]}", mode="w")
    outs[i].set_chromosomes(chroms)

  for (chrom, length) in info.pairs:
    d4s.check_same_lengths(chrom, length)
    let t0 = cpuTime()
    var j = 0
    for i in countup(0'u32, length, block_size):
      d4s.summarize_block(chrom, i, min(i + block_size, length), pctiles, outs)
      if j mod 5 == 0:
        echo &"{chrom}:{min(i + blocksize, length)} time:{cpuTime() - t0:.2f} bases/second: {min(i + blocksize, length).float / (cpuTime() - t0):.0f}"
      j.inc

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
