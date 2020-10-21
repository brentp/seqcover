import binaryheap
import json
import sequtils
import strformat
import strutils
import d4
export d4
import tables
import ./typeenum
import hts/fai

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


proc `$`*(t:Transcript): string =
  result = &"Transcript{system.`$`(t)}"

proc coding*(t:Transcript): bool =
  result = t.cdsstart != t.cdsend

type Gene* = object
  symbol*: string
  description*: string
  transcripts*: seq[Transcript]

proc `$`*(g:Gene): string =
  result = &"Gene{system.`$`(g)}"

type plot_coords* = object
  x*: seq[uint32]
  depths*: TableRef[string, seq[int32]]
  background_depths*: TableRef[string, seq[int32]]
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
    if t.coding and t.cdsstart < result.cdsstart:
      result.cdsstart = t.cdsstart

    if t.coding and t.cdsend > result.cdsend:
      result.cdsend = t.cdsend

    if t.txstart < result.txstart:
      result.txstart = t.txstart
    if t.txend > result.txend:
      result.txend = t.txend
    for ex in t.position:
      H.push(ex)

  result.position = newSeqOfCap[array[2, int]](4)
  if H.size == 0: return
  var A = H.pop
  var B: array[2, int]
  while H.size > 0:
    B = H.pop

    if A == B: continue

    if B[0] > A[1]:
      result.position.add(A)
      A = B
    else:
      A[1] = max(A[1], B[1])

  if result.position.len == 0 or result.position[^1] != A:
    if result.position.len > 0 and result.position[^1][0] <= A[0] and result.position[^1][1] >= A[1]: return
    # final interval needs to be merged/extended
    if result.position.len > 0 and result.position[^1][1] >= A[0]:
      result.position[^1][1] = max(A[1], result.position[^1][1])
    else:
      result.position.add(A)

proc find_offset*(u:Transcript, pos:int, extend:int, max_gap:int): int =
  doAssert pos >= u.txstart and pos <= u.txend, &"error can't translate position {pos} outside of unioned transcript ({u.txstart}, {u.txend})"
  if u.position.len == 0 or pos < u.position[0][0]:
    return pos - u.txstart

  result = u.position[0][0] - u.txstart

  for i, exon in u.position:
    if exon[0] > pos: break

    # add previous intron
    if i > 0:
      let intron = min(2 * extend + max_gap, exon[0] - u.position[i - 1][1])
      result += intron

    # add full size of this exon
    if exon[1] <= pos:
      let exon_bases = exon[1] - exon[0]
      result += exon_bases
    else:
      # exon contains current position
      result += (pos - exon[0])
      break

  if pos > u.position[^1][1]:
    result += pos - u.position[^1][1]


proc translate*(u:Transcript, o:Transcript, extend:uint32, max_gap:uint32): Transcript =
  ## given a unioned transcript, translate the positions in u to plot
  ## coordinates and genomic coordinates.

  var extend = extend.int
  result.transcript = o.transcript
  result.strand = o.strand
  result.`chr` = o.`chr`

  # needs some padding on left to go upstream, and use same param as for
  # introns, but here we limit to 1000 bases.
  let left = min(1000, extend.int)

  doAssert o.txstart >= u.txstart, $o
  result.txstart = left + u.find_offset(o.txstart, extend.int, max_gap.int)
  result.cdsstart = left + u.find_offset(o.cdsstart, extend.int, max_gap.int)

  # todo: this in n^2 (but n is small. iterate over uexons first and calc offsets once)?
  for i, o_exon in o.position:
    let l_off = left + u.find_offset(o_exon[0], extend.int, max_gap.int)
    let r_off = left + u.find_offset(o_exon[1], extend.int, max_gap.int)
    doAssert r_off - l_off == o_exon[1] - o_exon[0], $(r_off, l_off, o_exon, i)

    result.position.add([l_off, r_off])

  result.cdsend = left + u.find_offset(o.cdsend, extend.int, max_gap.int)
  result.txend = left + u.find_offset(o.txend, extend.int, max_gap.int)


proc `%`*[T](table: TableRef[string, T]): JsonNode =
  result = json.`%`(table[])

proc get_chrom(chrom:string, dp:var Cover, fai:Fai): string =
  ## add or remove "chr" to match chromosome names.
  var chroms = dp.chromosomes(fai)
  if chrom in chroms: return chrom
  if chrom[0] != 'c' and ("chr" & chrom) in chroms:
    result = "chr" & chrom
  elif chrom[0] == 'c' and chrom.len > 3 and chrom[1] == 'h' and chrom[2] == 'r' and chrom[3..chrom.high] in chroms:
    result = chrom[3..chrom.high]
  else:
    raise newException(KeyError, "chromosome not found:" & chrom)

proc exon_plot_coords*(tr:Transcript, dps:TableRef[string, Cover], backgrounds:TableRef[string, Cover], extend:uint32, max_gap:uint32, fai:Fai, utrs:bool=true): plot_coords =
  ## extract exonic depths for the transcript, extending into intron and
  ## up/downstream. This handles conversion to plot coordinates by removing
  ## introns. g: is the actual genomic coordinates.
  var chrom = tr.`chr`
  var dp: Cover
  for k, v in dps:
    dp = v
    break
  if dps.len > 0: chrom = chrom.get_chrom(dp, fai)
  let left = max(0, tr.txstart - min(1000, extend.int))

  result.depths = newTable[string, seq[int32]]()
  result.background_depths = newTable[string, seq[int32]]()

  var right = if tr.position.len > 0: tr.position[0][0].int else: tr.txstart
  var stop = tr.txend + min(1000, extend.int)

  if utrs:
    result.g = toSeq(left.uint32 ..< right.uint32)
    result.x = toSeq(0'u32 ..< result.g.len.uint32)
    when defined(debug):
      stderr.write_line &"utr: <adding> {result.g[0]} ..< {result.g[^1]}"
    for sample, dp in dps.mpairs:
      var x = newSeqUninitialized[int32](result.g[^1] - result.g[0])
      dp.values(chrom, result.g[0], x)
      result.depths[sample] = x
    for lvl, dp in backgrounds.mpairs:
      var x = newSeqUninitialized[int32](result.g[^1] - result.g[0])
      dp.values(chrom, result.g[0], x)
      result.background_depths[lvl] = x

  var lastx:uint32
  var lastg:uint32

  for i, p in tr.position:

    lastx = result.x[^1] + 1
    lastg = result.g[^1] + 1
    when defined(debug):
      stderr.write_line "\nbefore i:", $i, &" exon:{p} lastx: {lastx} lastg: {lastg}"


    # maxes and mins prevent going way past end of gene with huge extend value.
    let left = max(lastg, p[0].uint32 - (if i == 0: 0'u32 else: min(p[0].uint32, extend)))
    let right = min(stop.uint32, max(left, p[1].uint32 + (if i == tr.position.high: 0'u32 else: extend)))

    let isize = right.int - left.int
    if isize <= 0: continue
    let size = isize.uint32

    # insert value for missing data to interrupt plot
    if i > 0:
      let gap = min(max_gap, left - lastg)
      when defined(debug):
        stderr.write_line &"[seqcover] gap: {gap}"
      if gap > 0:
        result.x.add(lastx-1)
        result.g.add(lastg-1)
        lastx += gap
        result.x.add(lastx-1)
        result.g.add(lastg-1)
        for sample, dp in dps.mpairs:
          result.depths[sample].add(@[int32.low, int32.low])
        for lvl, dp in backgrounds.mpairs:
          result.background_depths[lvl].add(@[int32.low, int32.low])

      when defined(debug):
        stderr.write_line &"<adding> gap at {lastg} (x:{lastx})"

    when defined(debug):
      stderr.write_line "i:", $i, &"exon:{p}", &" <adding> {left} ..< {right}"
    result.g.add(toSeq(left..<right))
    result.x.add(toSeq((lastx..<(lastx + size))))

    var x = newSeqUninitialized[int32](right - left)
    for sample, dp in dps.mpairs:
      dp.values(chrom, left, x)
      result.depths[sample].add(x)
    for lvl, dp in backgrounds.mpairs:
      dp.values(chrom, left, x)
      result.background_depths[lvl].add(x)

  if utrs:
    lastx = result.x[^1] + 1
    lastg = result.g[^1] + 1
    let left = max(lastg, tr.cdsend.uint32)
    let right = max(left, stop.uint32) # already added extend to stop
    let size = right - left
    if size > 0:
      result.g.add(toSeq(left..<right))
      result.x.add(toSeq((lastx..<(lastx + size))))

      var x = newSeqUninitialized[int32](right - left)
      for sample, dp in dps.mpairs:
        dp.values(chrom, left, x)
        result.depths[sample].add(x)
      for lvl, dp in backgrounds.mpairs:
        dp.values(chrom, left, x)
        result.background_depths[lvl].add(x)

  doAssert result.x.len == result.g.len


proc plot_data*(g:Gene, d4s:TableRef[string, Cover], backgrounds:TableRef[string, Cover], extend:uint32, max_gap:uint32, fai:Fai, utrs:bool=true): GenePlotData =
  result.description = g.description
  result.symbol = g.symbol
  result.unioned_transcript = g.transcripts.union

  result.plot_coords = result.unioned_transcript.exon_plot_coords(d4s, backgrounds, extend, max_gap, fai, utrs)
  for t in g.transcripts:
    if t.`chr` != result.unioned_transcript.`chr`: continue
    result.transcripts.add(result.unioned_transcript.translate(t, extend=extend, max_gap=max_gap))

  result.unioned_transcript = result.unioned_transcript.translate(result.unioned_transcript, extend=extend, max_gap=max_gap)

