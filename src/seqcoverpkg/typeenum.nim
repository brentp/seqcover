import d4
import hts
import strutils

type CoverType* {.pure.} = enum
  D4
  BED

type Cover* = object
  case kind:CoverType
  of CoverType.D4:
    d4*: D4
  of CoverType.BED:
    bed*: BGZI

type iv = object
  start:int32
  stop:int32
  depth:int32

proc open_dp*(path:string): Cover =
  if path.endsWith(".d4"):
    result = Cover(kind:CoverType.D4)
    doAssert result.d4.open(path), "seqcover: error opening d4 file:" & path
  else:
    result = Cover(kind:CoverType.BED)
    doAssert result.bed.open(path), "seqcover: error opening bed.gz file:" & path

proc close*(c:var Cover) =
  case c.kind
  of CoverType.D4:
    c.d4.close()
  of CoverType.BED:
    discard c.bed.close

proc parse(s:string): iv {.noInit.} =
  var i = -1
  for tok in s.split('\t', maxsplit=4):
    i.inc
    if i == 0: continue
    if i == 1:
      result.start = parseInt(tok).int32
    if i == 2:
      result.stop = parseInt(tok).int32
    if i == 3:
      result.depth = parseInt(tok).int32
      break

proc values*(c:var Cover, chrom:string, start:uint32, values: var seq[int32]) =
  case c.kind
  of CoverType.D4:
    c.d4.values(chrom, start, values)
  of CoverType.BED:
    zeroMem(values[0].addr, values.len * sizeof(int32))
    for s in c.bed.query(chrom, start.int64, start.int64 + values.len.int64):
      let iv = s.parse
      for p in max(0, iv.start - start.int32)..<(min(values.len.int32, iv.stop)):
        values[p] = iv.depth

proc chromosomes*(c:var Cover, fai:Fai): OrderedTableRef[string, uint32] =
  case c.kind
  of CoverType.D4:
    result = c.d4.chromosomes
  of CoverType.BED:
    result = newOrderedTable[string, uint32]()
    for i in 0..<fai.len:
      var cname = fai[i]
      result[cname] = fai.chrom_len(cname).uint32

when isMainModule:
  import os

  var x = paramStr(1)

  let y = if x == "bed":
    Cover(kind:CoverType.BED)
  else:
    Cover(kind:CoverType.D4)

  echo y
