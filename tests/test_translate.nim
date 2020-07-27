import unittest
import seqcoverpkg/transcript
import seqcoverpkg/utils

var extend = 100
var max_gap = 100

suite "translate suite":
  test "bug case":
    var u = Transcript(cdsstart: 21505, cdsend: 31930, chr: "X", position: @[[19450, 21772], [24664, 24871], [25019, 25152]], strand: -1, transcript: "union", txstart: 19450, txend: 35554)
    var t = Transcript(cdsstart: 21505, cdsend: 31930, chr: "X", position: @[[19450, 21772], [24664, 24871], [25019, 25152]], strand: -1, transcript: "NM_002641", txstart: 19450, txend: 35554)
    echo u.translate(t, extend.uint32, max_gap.uint32)

  #[
  test "offset first exon":
    var u = Transcript(cdsstart: 21505, cdsend: 31930, chr: "X", position: @[[19450, 21772], [24664, 24871], [25019, 25152], [25913, 26046], [31215, 31992], [35500, 35554]], strand: -1, transcript: "union", txstart: 19450, txend: 35554)

    echo "first exon:", u.position[0][1] - u.position[0][0]
    echo u.find_offset(u.position[0][1], extend, max_gap)
    echo u.find_offset(u.position[1][0], extend, max_gap)

  ]#
