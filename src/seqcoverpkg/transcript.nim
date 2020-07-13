import binaryheap
import strformat


type Transcript* = object
  cdsstart*: int
  cdsend*: int
  `chr`*: string
  position*: seq[array[2, int]]
  strand*: int
  transcript*: string
  txstart*:int
  txend*:int

type Gene* = object
  symbol*: string
  description*: string
  transcripts*: seq[Transcript]

proc union*(trs:seq[Transcript]): Transcript =
  result = trs[0]
  var H = newHeap[array[2, int]] do (a, b: array[2, int]) -> int:
    if a[0] == b[0]: return a[1] - b[1]
    return a[0] - b[0]

  for t in trs:
    if t.`chr` != result.`chr`: continue
    if t.cdsstart < result.cdsstart:
      result.cdsstart = t.cdsstart
      result.cdsend = t.cdsend
      result.txstart = t.txstart
      result.txend = t.txend
    for ex in t.position:
      H.push(ex)

  result.position = newSeqOfCap[array[2, int]](4)
  var A = H.pop
  var B: array[2, int]
  while H.size > 0:
    B = H.pop

    if A == B: continue

    if B[0] > A[1]:
      result.position.add(A)
      A = B
    else:
      A[1] = B[1]

  if result.position.len == 0 or result.position[^1] != A:
    result.position.add(A)
proc `$`*(t:Transcript): string =
  result = &"Transcript{system.`$`(t)}"

proc `$`*(g:Gene): string =
  result = &"Gene{system.`$`(g)}"

