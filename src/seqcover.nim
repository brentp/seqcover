import argparse
import strformat
import tables
import sequtils
import algorithm
import json
import hts/fai
import tables
import os
import d4
import ./seqcoverpkg/utils
import ./seqcoverpkg/transcript
import ./seqcoverpkg/background
import ./seqcoverpkg/seqcover_html

proc get_pctile(path:string): int =
  let parts = path.split("_")
  if parts.len == 0:
    raise newException(OSError, &"[seqcover] path {path} doesn't have expected naming. use files generated by seqcover background")
  let last = parts[^1].split(".")[0]
  if last[0] != 'p':
    raise newException(OSError, &"[seqcover] path {path} doesn't have expected naming. use files generated by seqcover background")
  return parseInt(last[1..last.high])


proc read_backgrounds(dir:string): TableRef[int, D4] =
  if dir == "": return
  result = newTable[int, D4]()
  if not dirExists(dir):
    raise newException(OSError, "[seqcover] directory {dir} not found")
  for path in (&"{dir}/seqcover_*.d4").walkFiles:
    var d:D4
    if not d.open(path):
      quit &"error reading {path}"
    result[get_pctile(path)] = d

  if result.len == 0:
    var msg = &"[seqcover] no background d4 files found in {dir}"
    raise newException(OSError, msg)
  stderr.write_line &"[seqcover] read {result.len} background percentiles"

proc report_main() =
  let p = newParser("seqcover report"):
    option("--background", default="", help="optional path to d4 file created with seqcover background")
    option("--genes", default="", help="comma-delimited list of genes for initial report")
    option("--fasta", default="", help="required path to fai indexed fasta file")
    option("-r", "--report-path", default="seqcover_report.html", help="path to html report to be written")
    option("-t", "--transcripts-file", default="", help="path to transcript file for use if no internet connection (can be made with the save-transcripts option)")
    option("--extend-intron", default="10",help="The number of nucleotides to extend into the intron")
    flag("--hg19", help="coordinates are in hg19/GRCh37 (default is hg38).")
    arg("samples", nargs= -1, help="d4 files, bed files or a glob of d4 or bed files")

  var argv = commandLineParams()
  if len(argv) > 0 and argv[0] == "report":
    argv = argv[1..argv.high]
  if len(argv) == 0:
    argv.add("--help")

  var opts = p.parse(argv)
  if opts.help:
    quit 0
  if opts.fasta == "":
    echo p.help
    stderr.write_line "[seqcover] --fasta argument is required."
    quit 1

  var fa:Fai
  if not fa.open(opts.fasta):
    quit "[seqcover] couldn't open fasta file"

  var extend_intron = parseInt(opts.extend_intron) 
  if extend_intron < 0:
    echo p.help
    stderr.write_line "--extend-intron must be >= 0"
    quit 1

  var backgrounds = read_d4s_to_table(@[opts.background])
  var sample_d4s = read_d4s_to_table(opts.samples)
  stderr.write_line &"[seqcover] read {sample_d4s.len} sample coverage files"
  var gpt: seq[GenePlotData]
  
  var genes = get_genes(opts.genes.split(","), hg19=opts.hg19, transcript_file=opts.transcripts_file)
      
  for gene in genes:
      var u = gene.transcripts.union
      echo &"{u.chr}\t{u.txstart - 500}\t{u.txend + 500}"
      var pd = gene.plot_data(sample_d4s, backgrounds, extend=extend_intron.uint32, fai=fa, max_gap=50)
      gpt.add(pd)

  gpt.sort(proc(a, b: GenePlotData): int = cmp(a.symbol, b.symbol))

  write_html(opts.report_path, gpt)
  stderr.write_line &"[seqcover] wrote report to:{opts.report_path}"

proc save_transcripts_main() =
  let p = newParser("seqcover save-transcripts"):
    option("--genes", default="", help="comma-delimited list of genes for initial report")
    option("-o", "--output-path", default="transcripts.json", help="path to transcript file to be written")
    flag("--hg19", help="coordinates are in hg19/GRCh37 (default is hg38).")
  
  var argv = commandLineParams()
  if len(argv) > 0 and argv[0] == "save-transcripts":
    argv = argv[1..argv.high]
  if len(argv) == 0:
    argv.add("--help")

  var opts = p.parse(argv)
  if opts.help:
    quit 0
  
  let g = get_genes(opts.genes.split(","), hg19=opts.hg19)
  writeFile(opts.output_path, $(pretty(%*g)))
  

proc main() =
  type pair = object
    fn: proc()
    description: string

  var dispatcher = {
    "generate-background": pair(fn:generate_background_main, description: "generate background file(s) from a set of samples"),
    "report": pair(fn:report_main, description: "create an HTML report from a set of sample coverage files"),
    "save-transcripts": pair(fn:save_transcripts_main, description: "create a json-file with transcripts that can be used as input for report if you cannot access mygene.info")
  }.toOrderedTable

  var args = commandLineParams()
  if len(args) > 0 and args[0] in dispatcher:
    dispatcher[args[0]].fn()
    return

  if len(args) == 0 or args[0] in ["-h", "--help"]:
    stdout.write_line "Commands: "
    for k, v in dispatcher:
      echo &"  {k:<19}:   {v.description}"
  else:
    echo &"unknown program '{args[0]}'"
    quit ""


when isMainModule:
  main()
