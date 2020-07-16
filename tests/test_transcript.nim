import unittest
import seqcoverpkg/transcript
import seqcoverpkg/utils
import tables
import json

suite "transcript suite":
  test "simple plot exon coords":

    #exon_plot_coords*(tr:Transcript, dps:TableRef[string, D4], extend:uint32=10): tuple[x:seq[uint32], depths: TableRef[string, seq[int32]], g: seq[uint32]] =
    var tr = Transcript(txStart: 10, txEnd: 100, cdsStart: 15, cdsEnd: 95, position: @[[15, 25], [85, 91]])
    echo "UTR5:", tr.UTR5
    echo "UTR3:", tr.UTR3
    var dps = newTable[string, D4]()
    var res = tr.exon_plot_coords(dps, 10)
    check res.x.len == res.g.len


  test "generate data":
    var trs = @[Transcript(cdsstart: 15321505, cdsend: 15331930, chr: "X", position: @[[15319450, 15321772], [15324664, 15324871], [15325019, 15325152], [15325913, 15326046], [15331215, 15331992], [15335500, 15335554]], strand: -1, transcript: "union", txstart: 15319450, txend: 15335554),
                Transcript(cdsstart: 51662817, cdsend: 51807429, chr: "12", position: @[[51591232, 51591359], [51662763, 51663093], [51684173, 51684292], [51686367, 51686457], [51687090, 51687219], [51688757, 51688849], [51689004, 51689096], [51699569, 51699791], [51701143, 51701207], [51702772, 51702914], [51705416, 51705623], [51706421, 51706715], [51721545, 51721908], [51745902, 51746035], [51751354, 51751593], [51762502, 51762676], [51765670, 51766027], [51768864, 51769335], [51769867, 51769985], [51770528, 51770683], [51774188, 51774362], [51780648, 51780771], [51786541, 51786826], [51788694, 51788748], [51789280, 51789418], [51790397, 51790502], [51794370, 51794641], [51806281, 51812864]], strand: 1, transcript: "union", txstart: 51591235, txend: 51812864),
                Transcript(cdsstart: 63406643, cdsend: 63472463, chr: "20", position: @[[63400207, 63407375], [63408412, 63408536], [63413449, 63413581], [63414087, 63414193], [63414902, 63415126], [63419618, 63419672], [63424176, 63424206], [63428366, 63428435], [63431339, 63431369], [63433677, 63433903], [63438624, 63438720], [63439597, 63439708], [63442405, 63442531], [63444658, 63444834], [63445237, 63445364], [63446746, 63446837], [63472167, 63472655]], strand: -1, transcript: "union", txstart: 63400207, txend: 63472655)
                ]


    var d4s = read_d4s_to_table(@["./d4s/HG00105.final.d4",
                                  "./d4s/HG00103.final.d4",
                                  "./d4s/HG00112.final.d4",
                                  "./d4s/HG00109.final.d4",
                                  "./d4s/HG00100.final.d4"])

    var coords:seq[plot_coords]
    for t in trs:
      var n = 0
      for p in t.position:
        n += p[1] - p[0]
      stderr.write_line "size:", $n
      coords.add(t.exon_plot_coords(d4s, 10))
    #echo ((%coords))


  test "translate":

    var u = Transcript(txstart:55, cdsstart:85, position: @[[85, 96], [122, 137]], cdsend: 137, txend: 152)
    var o = Transcript(txstart:55, cdsstart:85, position: @[[85, 96], [122, 137]], cdsend: 137, txend: 152)

    var t = u.translate(o)
    check t.position == @[[40, 51], [71, 86]]
    check t.txstart == 10
    check t.txend == 101

    o = Transcript(txstart:55, cdsstart:85, position: @[[89, 94], [132, 137]], cdsend: 137, txend: 152)
    t = u.translate(o)
    check t.position ==  @[[44, 49], [81, 86]]
    check t.txstart == 10
    check t.txend == 101
    #echo cast[seq[byte]](coords[0].g)