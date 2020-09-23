var assert = require('assert');
let data = require("./data.js")
let transcript = require("../src/seqcoverpkg/transcript.js")

let tr = new transcript.Transcript(data.tr_data)

describe('Transcript', function () {
  describe('cdsstart', function () {
    it('should equal cdsstart', function () {
        assert.equal(tr.cdsstart, data.tr_data.cdsstart)
    });
  });

  describe('transcript tests', function () {
      it('inferred genomic coords should match actual', function () {
          let parts = tr.parts()
          parts.forEach((p, i) => {
              //  u:Transcript(cdsstart: 15321505, cdsend: 15335554, chr: "X", position: @[[15319450, 15321772], [15324664, 15324871], [15325019, 15325152], [15325913, 15326046], [15331215, 15331992], [15335500, 15335554]], strand: -1, transcript: "union", txstart: 15319450, txend: 15335554)
              let s = p.hoverinfo(data.xs, data.gs)
              if(i == 0){
                  assert.equal(s.start, 15319450)
                  assert.equal(s.stop, 15321505)
              } else if (i == parts.length) {
                  assert.equal(s.start, 15335554)
                  assert.equal(s.stop, 15335554)

              } else {

                  if(s.type == transcript.FeatureType.CDS) {
                      assert(s.start >= 15321505, s.start.toString() + " " + i.toString())
                      assert(s.stop <= 15335554, s.stop.toString() + " " + i.toString())
                  }
                  if(s.type == transcript.FeatureType.EXON) {
                      assert(s.start >= 15319450, s.start.toString() + " " + i.toString())
                      assert(s.stop <= 15335554, s.stop.toString() + " " + i.toString())
                  }
              }

          })

      });
      it('traces should work', function () {
          let traces = tr.traces(0, data.xs, data.gs)
          assert.equal(traces.length, 3)
          assert.equal(traces[0].x.length, 2, "only need start and end")
          assert.equal(traces[0].y[0], 0, "y offset should be 0")
          assert.equal(traces[0].genome_x.length, 2, "only need start and end")

          assert.equal(traces[0].name, tr.name)

      })

      it('stats returns expected values', function() {
          var depths = {
              "sampleA": Array(60).fill(10),
              "sampleB": Array(60).fill(0),
          }
          var background_depths = {
              "p10": Array(60).fill(5),
              "p90": Array(60).fill(90),
          }

          background_depths.p10[4] = 500;
          //{ sampleA: { lt_background: 0, low: 0, mean: 10, median: 10 },
          //  sampleB: { lt_background: 20, low: 20, mean: 0, median: 0 } }

          let stats = tr.stats([{start:20, stop:40}], depths, background_depths, 5);
          assert.equal(stats.sampleA.median, 10)
          assert.equal(stats.sampleB.median, 0)

          assert.equal(stats.sampleA.low, 0)
          assert.equal(stats.sampleB.low, 20)

          assert.equal(stats.sampleA["lt-background"], 0)
          assert.equal(stats.sampleB["lt-background"], 20)



      })
  });

});
