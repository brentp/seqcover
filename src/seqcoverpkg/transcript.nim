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

# see: https://github.com/nim-lang/Nim/issues/15025
proc `%`*(a:array[2, int]): JsonNode =
  result = newJArray()
  result.add(newJint(a[0]))
  result.add(newJint(a[1]))


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

type plot_coords* = object
  x*: seq[uint32]
  depths*: TableRef[string, seq[int32]]
  g*: seq[uint32]

type GenePlotData* = object
  plot_coords*: plot_coords
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
      #doAssert u_exon[0] <= o_exon[0] and u_exon[1] >= o_exon[1], $(u, o) & $(u_exon, o_exon)

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


proc `%`*[T](table: TableRef[string, T]): JsonNode =
  result = json.`%`(table[])

proc get_chrom(chrom:string, dp:D4): string =
  ## add or remove "chr" to match chromosome names.
  if chrom in dp.chromosomes: return chrom
  if chrom[0] != 'c' and ("chr" & chrom) in dp.chromosomes:
    result = "chr" & chrom
  elif chrom[0] == 'c' and chrom.len > 3 and chrom[1] == 'h' and chrom[2] == 'r' and chrom[3..chrom.high] in dp.chromosomes:
    result = chrom[3..chrom.high]
  else:
    raise newException(KeyError, "chromosome not found:" & chrom)

proc exon_plot_coords*(tr:Transcript, dps:TableRef[string, D4], extend:uint32=10, utrs:bool=true): plot_coords =
  ## extract exonic depths for the transcript, extending into intron and
  ## up/downstream. This handles conversion to plot coordinates by removing
  ## introns. g: is the actual coordinates.
  var chrom = tr.`chr`
  var dp: D4
  for k, v in dps:
    dp = v
    break
  if dps.len > 0: chrom = chrom.get_chrom(dp)
  let left = max(0, tr.txstart - extend.int)

  result.depths = newTable[string, seq[int32]]()

  var  lutr = tr.UTR_left
  lutr[0] = max(0, lutr[0] - extend.int)
  var rutr = tr.UTR_right
  rutr[1] += extend.int

  if utrs:
    result.g = toSeq(lutr[0].uint32 ..< lutr[1].uint32)
    result.x = toSeq(0'u32 ..< result.g.len.uint32)
    for sample, dp in dps.mpairs:
        result.depths[sample] = dp.values(chrom, result.g[0], result.g[^1])

  var lastx:uint32
  var lastg:uint32

  for i, p in tr.position:

    lastx = result.x[^1] + 1
    lastg = result.g[^1] + 1

    # insert value for missing data to interrupt plot
    if i > 0:
      result.x.add(result.x[^1])
      result.g.add(result.g[^1])

    let left = max(lastg, p[0].uint32 - (if i == 0: 0'u32 else: extend))
    let right = max(left, p[1].uint32 + (if i == tr.position.high: 0'u32 else: extend))
    let size = right - left
    if size <= 0: continue

    result.g.add(toSeq(left..<right))
    result.x.add(toSeq((lastx..<(lastx + size))))

    for sample, dp in dps.mpairs:
      if i > 0:
        result.depths[sample].add(int32.low)
      result.depths[sample].add(dp.values(chrom, left, right))

  if utrs:
    lastx = result.x[^1] + 1
    lastg = result.g[^1] + 1
    let left = max(lastg, rutr[0].uint32)
    let right = max(left, rutr[1].uint32 + extend)
    let size = right - left
    if size > 0:
      result.g.add(toSeq(left..<right))
      result.x.add(toSeq((lastx..<(lastx + size))))

    for sample, dp in dps.mpairs:
      result.depths[sample].add(dp.values(chrom, left, right))

  doAssert result.x.len == result.g.len


proc plot_data*(g:Gene, d4s:TableRef[string, D4], extend:uint32=10, utrs:bool=true): GenePlotData =
  result.description = g.description
  result.symbol = g.symbol
  result.transcripts = g.transcripts
  result.unioned_transcript = g.transcripts.union
  result.plot_coords = result.unioned_transcript.exon_plot_coords(d4s, extend, utrs)

proc `$`*(t:Transcript): string =
  result = &"Transcript{system.`$`(t)}"

proc `$`*(g:Gene): string =
  result = &"Gene{system.`$`(g)}"
