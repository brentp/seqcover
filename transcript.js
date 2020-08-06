"use strict";

function binary_search(A, v) {
    var result = 0;
    var j = A.length
    while (j != 0) {
        let step = j >> 1
        let pos = result + step;
        if (A[pos] < v) {
            result = pos + 1
            j -= step + 1
        } else {
            j = step
        }

    }
    return result
}

// enum
const FeatureType = Object.freeze({
    EXON:   Symbol("exon"),
    CDS:  Symbol("CDS"),
    UTR: Symbol("UTR")
})

class Feature {
  constructor(start, stop, type, transcript) {
      this.start = start
      this.stop = stop
      this.type = type
      this.transcript = transcript
  }
  get coding() {
      return this.type == FeatureType.CDS
  }


  shape(xs, gs) {
      // xs is plot coordinates, gs is genome coordinates
      // get genomic coordinate by finding index of plot-coord
      // and then looking it up in genomic array
      var start = gs[binary_search(xs, this.start)]
      var stop = gs[binary_search(xs, this.stop)]
      return {
          "transcript": this.transcript.data.transcript,
          "strand": this.transcript.strand,
          "start": start,
          "stop": stop,
          "type": this.type.toString()
      }
  }
}

class Transcript {
    constructor(data) {
        this.data = data
    }
    get cdsstart() {
        return this.data.cdsstart
    }
    get cdsend() {
        return this.data.cdsend
    }
    get chr() {
        return this.data.chr
    }

    get exons() {
        return this.data.position
    }

    get position() {
        return this.data.position
    }
    get strand() {
        return this.data.strand == -1 ? "-" : "+"
    }

    get name() {
        return this.data.transcript
    }

    get txstart() {
        return this.data.txstart
    }
    get txend() {
        return this.data.txend
    }

    parts() {
        // return CDS,exon,UTR in an array. exon and CDS are potentially
        // (often) duplicated.
        var that = this
        var result = []
        result.push(new Feature(this.data.txstart, this.data.cdsstart, FeatureType.UTR, that))
        this.data.position.forEach((exon, i) => {
            result.push(new Feature(exon[0], exon[1], FeatureType.EXON, that))
            if(exon[1] < this.data.cdsstart || exon[0] > this.data.cdsend) {
                // continue
            } else {
                result.push(new Feature(Math.max(this.data.cdsstart, exon[0]), Math.min(this.data.cdsend, exon[1]), FeatureType.CDS, that))
            }

        })
        result.push(new Feature(this.data.cdsend, this.data.txend, FeatureType.UTR, that))
        return result.filter(f => f.stop - f.start > 0)

    }
}

if(require.main === module) {
    var tr_data = {"cdsstart":2065,"cdsend":3986,"chr":"X","position":[[10,2332],[2402,2609],[2679,2812],[2882,3015],[3085,3862],[3932,3986]],"strand":-1,"transcript":"union","txstart":10,"txend":3986}
    var tr = new Transcript(tr_data)

    // xs and gs data for testing.
    let data = require("./data.js")

    var tr = new Transcript(tr_data)
    console.log(tr.cdsstart)
    tr.parts().forEach(p => console.log(p.shape(data.xs, data.gs)))

}
