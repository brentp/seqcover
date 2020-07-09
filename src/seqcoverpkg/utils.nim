import os
import d4


proc get_glob_samples*(paths: seq[string]): seq[string] =
  result = newSeqOfCap[string](paths.len)
  for i, p in paths:
    var n = 0
    for w in p.walkFiles:
      n.inc
      result.add(w)
    if n == 0:
      raise newException(OSError, "[seqcover]: file:" & p & " not found")


proc read_d4s*(paths: seq[string]): seq[D4] =
  result = newSeq[D4](paths.len)
  for i, p in paths:
    if not result[i].open(p):
      raise newException(OSError, "[seqcover] couldn't open d4 file:" & p)

proc check_same_lengths*(d4s: seq[D4], chrom: string, length: uint32) =
  for d in d4s:
    doAssert d.chromosomes[chrom] == length, "[seqcover] differing chromosome lengths among samples"
