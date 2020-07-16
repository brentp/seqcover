import binaryheap
import json
import sequtils
import strformat
import strutils
import d4
export d4
import tables

type Transcript* = object
  cdsstart*: int
  cdsend*: int
  `chr`*: string
  position*: seq[array[2, int]]
  strand*: int
  transcript*: string
  txstart*:int
  txend*:int

proc UTR5*(t:Transcript): array[2, int] =
  if t.strand >= 0:
    return [t.txstart, t.cdsstart]
  return [t.cdsend, t.txend]

proc UTR3*(t:Transcript): array[2, int] =
  if t.strand >= 0:
    return [t.cdsend, t.txend]
  return [t.txstart, t.cdsstart]

proc UTR_left(t:Transcript): array[2, int] {.inline.} =
  result = [t.txstart, t.cdsstart]

proc UTR_right(t:Transcript): array[2, int] {.inline.} =
  result = [t.cdsend, t.txend]

type Gene* = object
  symbol*: string
  description*: string
  transcripts*: seq[Transcript]

type GenePlotData* = object
  x*: seq[uint32]
  depths*: TableRef[string, seq[uint32]]
  g*: seq[uint32]
  symbol*: string
  description*: string
  unioned_transcript*: Transcript
  transcripts*: seq[Transcript]

proc union*(trs:seq[Transcript]): Transcript =
  result = trs[0]
  result.transcript = "union"
  ## TODO: heap is not needed here. just sort and check i vs i+1
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

proc translate*(u:Transcript, o:Transcript, extend:uint32|uint=10): Transcript =
  ## given a unioned transcript, translate the positions in u to plot
  ## coordinates and genomic coordinates.

  # exons and UTRs
  var extend = extend.int
  result.transcript = o.transcript
  result.strand = o.strand
  result.`chr` = o.`chr`

  result.txstart = (o.txstart - u.txstart) + extend
  result.cdsstart = (o.cdsstart - u.txstart) + extend

  # todo: this in n^2 (but n is small. iterate over uexons first and calc
  # offsets once)?
  for i, o_exon in o.position:
    var u_off = extend + (result.cdsstart - result.txstart)
    # increase u_off until we find the u_exon that encompasses this one.
    var u_i = 1
    while u_i < o.position.len:
      let u_exon = u.position[u_i]
      if u_exon[0] >= o_exon[1]: break
      doAssert u_exon[0] <= o_exon[0] and u_exon[1] >= o_exon[1]

      u_off += (u.position[u_i - 1][1] - u.position[u_i - 1][0])
      u_off += min(2 * extend, u_exon[0] - u.position[u_i - 1][1])
      u_i += 1
    u_i -= 1

    # handle o exon starting after start of u-exon
    if o.position.len > 0:
      u_off += o.position[u_i][0] - u.position[u_i][0]

    result.position.add([u_off, u_off + (o_exon[1] - o_exon[0])])

  result.cdsend = result.position[^1][1] + (o.cdsend - o.position[^1][1])
  result.txend = (o.txend - o.cdsend) + result.cdsend

type plot_coords* = object
  x*: seq[uint32]
  depths*: TableRef[string, seq[int32]]
  g*: seq[uint32]

proc `%`*[T](table: TableRef[string, T]): JsonNode =
  result = json.`%`(table[])

proc exon_plot_coords*(tr:Transcript, dps:TableRef[string, D4], extend:uint32=10): plot_coords =
  ## extract exonic depths for the transcript, extending into intron and
  ## up/downstream. This handles conversion to plot coordinates by removing
  ## introns. g: is the actual coordinates.
  var chrom = tr.`chr`
  var dp: D4
  var found = false
  for k, v in dps:
    dp = v
    found = true
    break
  if found and (chrom notin dp.chromosomes):
    if chrom[0] != 'c' and "chr" & chrom in dp.chromosomes:
      chrom = "chr" & chrom
    elif chrom[0] == 'c' and chrom.len > 3 and chrom[1] == 'h' and chrom[2] == 'r' and chrom[3..chrom.high] in dp.chromosomes:
      chrom = chrom[3..chrom.high]
    else:
      raise newException(KeyError, "chromosome not found:" & chrom)

  result.depths = newTable[string, seq[int32]]()

  let left = max(0, tr.txstart - extend.int)
  let right = tr.txend + extend.int

  var lpos: array[2,int]
  var rpos: array[2,int]
  if tr.strand >= 0:
    lpos = tr.UTR5
    lpos[0] -= extend.int
    rpos = tr.UTR3
    rpos[1] += extend.int
  else:
    lpos = tr.UTR3
    lpos[1] += extend.int
    rpos = tr.UTR5
    rpos[0] -= extend.int

  result.g = toSeq(max(0, lpos[0]).uint32..< lpos[1].uint32)
  #result.x = toSeq(0'u32 ..< result.g.len.uint32)
  for sample, dp in dps.mpairs:
      result.depths[sample] = dp.values(chrom, result.g[0], result.g[^1])

  var exons = tr.position
  exons.add(rpos)

  for p in exons:

    #var lastx = result.x[^1] + 1
    var lastg = result.g[^1] + 1

    #result.x.add(result.x[^1])
    result.g.add(result.g[^1])

    let left = max(lastg, p[0].uint32 - extend)
    let right = max(left, p[1].uint32 + extend)
    let size = right - left
    if size == 0: continue
    result.g.add(toSeq(left..<right))

    #result.x.add(toSeq((lastx..<(lastx + size))))
    #echo chrom, " ", left, "-", right

    for sample, dp in dps.mpairs:
      result.depths[sample].add(int32.low)
      result.depths[sample].add(dp.values(chrom, left, right))


proc plot_data*(g:Gene): GenePlotData =
  result.description = g.description
  result.symbol = g.symbol
  result.transcripts = g.transcripts
  result.unioned_transcript = g.transcripts.union


proc `$`*(t:Transcript): string =
  result = &"Transcript{system.`$`(t)}"

proc `$`*(g:Gene): string =
  result = &"Gene{system.`$`(g)}"



