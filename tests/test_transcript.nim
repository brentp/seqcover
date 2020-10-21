import unittest
import seqcoverpkg/transcript
import seqcoverpkg/utils
import hts/fai
import seqcoverpkg/typeenum
import tables
import json

suite "transcript suite":
  test "simple plot exon coords":

    #exon_plot_coords*(tr:Transcript, dps:TableRef[string, D4], extend:uint32=10): tuple[x:seq[uint32], depths: TableRef[string, seq[int32]], g: seq[uint32]] =
    var tr = Transcript(txStart: 10, txEnd: 100, cdsStart: 15, cdsEnd: 95, position: @[[15, 25], [85, 91]])
    var dps = newTable[string, Cover]()
    var bgs = newTable[string, Cover]()
    var fa:Fai
    var res = tr.exon_plot_coords(dps, bgs, 10, 10, fa)
    check res.x.len == res.g.len
    for name, depth in res.depths.pairs:
      check res.x.len == depth.len

  test "translate":

    var u = Transcript(txstart:55, cdsstart:85, position: @[[85, 96], [122, 137]], cdsend: 137, txend: 152)
    var o = u #Transcript(txstart:55, cdsstart:85, position: @[[85, 96], [122, 137]], cdsend: 137, txend: 152)


    var t = u.translate(o, 10, 100)
    check t.position == @[[40, 51], [77, 92]]

    check t.txstart == 10
    check t.txend == 107

    o = Transcript(txstart:55, cdsstart:85, position: @[[89, 94], [132, 137]], cdsend: 137, txend: 152)
    t = u.translate(o, 10, 100)
    check t.position ==   @[[44, 49], [87, 92]]
    check t.txstart == 10
    check t.txend == 107

  test "spanning exon":
    var g = Gene(symbol: "HNRNPA1", description: "heterogeneous nuclear ribonucleoprotein A1", transcripts: @[
    Transcript(cdsstart: 74591, cdsend: 78097, chr: "12", position: @[[78041, 78101], [78332, 80871]], strand: 1, transcript: "NM_002136", txstart: 74509, txend: 80871),
    Transcript(cdsstart: 80871, cdsend: 80871, chr: "12", position: @[[77595, 77751], [78041, 78101], [80445, 80516],  [80736, 80871]], strand: 1, transcript: "NR_135167", txstart: 74509, txend: 80871)
    #Transcript(cdsstart: 80871, cdsend: 80871, chr: "12", position: @[[77595, 77751], [78041, 78101], [80445, 80516], [80736, 80871]], strand: 1, transcript: "NR_135167", txstart: 74509, txend: 80871)
    ])

    var u = g.transcripts.union

    var tr = u.translate(g.transcripts[0], 10, 100)

    check tr == Transcript(cdsstart: 92, cdsend: 3428, chr: "12", position: @[[3372, 3432], [3552, 6091]], strand: 1, transcript: "NM_002136", txstart: 10, txend: 6091)

  test "no exons":

    var g = Gene(symbol: "CCDC39", description: "coiled-coil domain containing 39", transcripts: @[Transcript(cdsstart: 180614920, cdsend: 180679380, chr: "3", position: @[[180614007, 180615077], [180616280, 180616363], [180616515, 180616695], [180616825, 180616966], [180619258, 180619365], [180619810, 180619970], [180631468, 180631592], [180641992, 180642201], [180644119, 180644257], [180647078, 180647243], [180648164, 180648359], [180651400, 180651533], [180652162, 180652266], [180654761, 180654953], [180659451, 180659580], [180659676, 180659769], [180660569, 180660728], [180661860, 180662007], [180663866, 180663986], [180679290, 180679489]], strand: -1, transcript: "NM_181426", txstart: 180614007, txend: 180679489)])
    var u = g.transcripts.union

    for t in g.transcripts:
      var v = u.translate(t, 10, 100)
      check v.`chr` == g.transcripts[0].`chr`

