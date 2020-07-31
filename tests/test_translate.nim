import unittest
import seqcoverpkg/transcript
import seqcoverpkg/utils

var extend = 100
var max_gap = 100

suite "translate suite":
  test "bug case":
    var u = Transcript(cdsstart: 21505, cdsend: 31930, chr: "X", position: @[[19450, 21772], [24664, 24871], [25019, 25152]], strand: -1, transcript: "union", txstart: 19450, txend: 35554)
    var t = Transcript(cdsstart: 21505, cdsend: 31930, chr: "X", position: @[[19450, 21772], [24664, 24871], [25019, 25152]], strand: -1, transcript: "NM_002641", txstart: 19450, txend: 35554)
    var tr = u.translate(t, extend.uint32, max_gap.uint32)
    echo tr
    check tr.position == @[[100, 2422], [2722, 2929], [3077, 3210]]


  test "offset first exon":
    var u = Transcript(cdsstart: 21505, cdsend: 31930, chr: "X", position: @[[19450, 21772], [24664, 24871], [25019, 25152], [25913, 26046], [31215, 31992], [35500, 35554]], strand: -1, transcript: "union", txstart: 19450, txend: 35554)

    var l = u.find_offset(u.position[0][1], extend, max_gap)
    check l == 2322
    var r = u.find_offset(u.position[1][0], extend, max_gap)
    check r == l + 2 * extend + max_gap

  test "weird kcnq2":

    var o = Transcript(cdsstart: 33744, cdsend: 33900, chr: "20", position: @[[33677, 33903]], strand: -1, transcript: "NM_172109", txstart: 33677, txend: 72655)
    var u = Transcript(cdsstart: 06643, cdsend: 72463, chr: "20", position:
       @[[207, 7375],
         [8412, 8536],
         [13449, 13581],
         [14087, 14193],
         [14902, 15126],
         [19618, 19672],
         [24176, 24206],
         [28366, 28435],
         [31339, 31369],
         [33677, 33903],
         [38624, 38720],
         [39597, 39708],
         [42405, 42531],
         [44658, 44834], [45237, 45364], [46746, 46837], [72167, 72655]], strand: -1, transcript: "union", txstart: 207, txend: 72655)

    var ot = u.translate(o, extend.uint32, max_gap.uint32)
    var ut = u.translate(u, extend.uint32, max_gap.uint32)

    check ot.position[0] == ut.position[9]
    check ot.position[0][1] == ot.cdsend + 3
    check ot.position[0][0] - ot.cdsstart == o.position[0][0] - o.cdsstart

  test "translation doesn't change exon length and adjusts intron length correctly":

    var u = Transcript(cdsstart: 9551, cdsend: 26688, chr: "20", position: @[[8675, 9033], [9497, 9617], [13168, 13226], [13257, 13383], [13984, 14058], [15280, 15312], [18516, 18632]], strand: 1, transcript: "union", txstart: 8675, txend: 27449)

    var ut = u.translate(u, extend.uint32, max_gap.uint32)

    for i, e in ut.position:
      check e[1] - e[0] == u.position[i][1] - u.position[i][0]

      if i > 0:
        var l = ut.position[i-1]
        var r = e
        check r[0] - l[1] == min(2 * extend + max_gap, u.position[i][0] - u.position[i - 1][1])

    echo "ut:", ut


