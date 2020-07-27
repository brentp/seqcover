import unittest
import seqcoverpkg/transcript
import seqcoverpkg/utils

suite "translate suite":
  test "bug case":
    var o = Transcript(cdsstart: 15321505, cdsend: 15335513, chr: "X", position: @[[15335500, 15335554]], strand: -1, transcript: "NM_020473", txstart: 15319450, txend: 15335554)
    var u = Transcript(cdsstart: 15321505, cdsend: 15331930, chr: "X", position: @[[15319450, 15321772], [15324664, 15324871], [15325019, 15325152], [15325913, 15326046], [15331215, 15331992], [15335500, 15335554]], strand: -1, transcript: "union", txstart: 15319450, txend: 15335554)

    var extend = 100

    var r = Transcript()
    r.txstart = (o.txstart - u.txstart) + min(1000, extend.int)

    var o_exon = o.position[0]
    echo find_offset(o_exon, r, o, u, extend, extend)

