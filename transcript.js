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
    EXON: "exon",
    CDS:  "CDS",
    UTR: "UTR",
    TRANSCRIPT: "transcript"
})

const aesthetics = {
    TRANSCRIPT_COLOR: "rgb(65, 65, 65)",
    TRANSCRIPT_WIDTH: 8,
    EXON_COLOR: "rgb(105,105,105)",
    EXON_WIDTH: 14,
    CDS_COLOR:'rgb(155, 155, 155)',
    CDS_WIDTH: 20,
}

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


  hoverinfo(xs, gs) {
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

    overlaps(position) {
        // return parts() that overlap with this position
        let that = this
        var result = []
        if(position < this.txstart || position > this.txend) { return result}
        if(position < this.cdsstart) { result.push(new Feature(this.data.txstart, this.data.cdsstart, FeatureType.UTR, that)) }
        this.data.position.forEach((exon, i) => {
            if (exon[0] > position || exon[1] < position) {
                return 
            }
            result.push(new Feature(exon[0], exon[1], FeatureType.EXON, that))
            if(exon[1] < that.data.cdsstart || exon[0] > that.data.cdsend) {
                // continue
            } else {
                var f = new Feature(Math.max(this.data.cdsstart, exon[0]), Math.min(this.data.cdsend, exon[1]), FeatureType.CDS, that)
                if(f.stop - f.start > 0){
                    result.push(f)
                }
            }
        })
        if(position >= this.cdsend) {
            result.push(new Feature(this.data.cdsend, this.data.txend, FeatureType.UTR, that))
        }
        return result;
    }

    traces(y_offset, xs, gs) {

        function get_genomic_coord(x) {
            if(isNaN(x)){ return NaN }
            return gs[binary_search(xs, x)]
        }

        var transcript_trace = {name: this.data.transcript, x: [this.data.txstart, this.data.txend], y:[y_offset, y_offset],
             type: "scatter", mode: "lines", showlegend: false,
             line: { color: aesthetics.TRANSCRIPT_COLOR, width: aesthetics.TRANSCRIPT_WIDTH}}

        let parts = this.parts()

        var exon_trace = {name: this.data.transcript + " exons", x: [], y:[],
             type: "scatter", mode: "lines", showlegend: false,
             line: { color: aesthetics.EXON_COLOR, width: aesthetics.EXON_WIDTH}}
        parts.filter(p => p.type == FeatureType.EXON).forEach(e => {
            if((exon_trace.x.length) > 0) {
                exon_trace.x.push(NaN)
                exon_trace.y.push(y_offset)
            }
            exon_trace.x.push(e.start, e.stop)
            exon_trace.y.push(y_offset, y_offset)
        })

        var cds_trace = {name: this.data.transcript + " CDS", x: [], y:[],
             type: "scatter", mode: "lines", showlegend: false,
             line: { color: aesthetics.CDS_COLOR, width: aesthetics.CDS_WIDTH}}
        parts.filter(p => p.type == FeatureType.CDS).forEach(c => {
            if((cds_trace.x.length) > 0) {
                cds_trace.x.push(NaN)
                cds_trace.y.push(y_offset)
            }
            cds_trace.x.push(c.start, c.stop)
            cds_trace.y.push(y_offset, y_offset)
        })

        var result = [transcript_trace, exon_trace, cds_trace]
        result.forEach(trace =>  {
            trace.genome_x = trace.x.map(x => get_genomic_coord(x))
        })
        return result

    }
}


try {
    // node.js stuff to allow testing
    exports.Transcript = Transcript
    exports.FeatureType = FeatureType


    if(require.main === module) {
        // xs and gs data for testing.
        let data = require("./test/data.js")
        let tr = new Transcript(data.tr_data)

        tr.parts().forEach(p => console.log(p.hoverinfo(data.xs, data.gs)))

        console.log(tr.traces(0))

    }
} catch (e) {
    // browser
}
