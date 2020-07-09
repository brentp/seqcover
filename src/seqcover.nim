import argparse
import os

proc main() =
  let p = newParser("seqcover report"):
    option("--server-port", default="", help="optional port on which to start local server for interactive viewing")
    option("--backgrounds", default="", help="optional glob for a set of background samples (e.g. bg/*.d4)")
    option("--genes", default="", help="comma-delimited list of genes for initial report")
    arg("samples", nargs= -1, help="d4 files or a glob of d4 files")

  var argv = commandLineParams()
  if len(argv) > 0 and argv[0] == "report":
    argv = argv[1..argv.high]
  if len(argv) == 0:
    argv.add("--help")

  var opts = p.parse(argv)
  if opts.help:
    quit 0

  echo opts


when isMainModule:
  main()
