var nan = NaN; // hack to support json dumped with NaN values.
const gene_plot_height = 500
const color_list = [
    "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b",
    "#e377c2", "#7f7f7f", "#bcbd22", "#17becf", "#7CB5EC", "#434348",
    "#90ED7D", "#F7A35C", "#8085E9", "#F15C80", "#E4D354", "#2B908F",
    "#F45B5B", "#91E8E1", "#4E79A7", "#F28E2C", "#E15759", "#76B7B2",
    "#59A14F", "#EDC949", "#AF7AA1", "#FF9DA7", "#9C755F", "#BAB0AB",
]
const HOVER_TEMPLATE = '<b>position</b>:%{text}<br><b>depth</b>:%{y}<br>(debug) x: %{x}';


// https://stackoverflow.com/a/4198132
function getHashParams(h=false) {

    var hashParams = {};
    var e,
        a = /\+/g,  // Regex for replacing addition symbol with a space
        r = /([^&;=]+)=?([^&;]*)/g,
        d = function (s) { return decodeURIComponent(s.replace(a, " ")); },
        q = h ? h : window.location.hash.substring(1)

    while (e = r.exec(q))
       hashParams[d(e[1])] = d(e[2]);

    return hashParams;
}

function setHashParams(obj) {
    var old = getHashParams();
    for(var k in obj){
        if(obj[k] == null){
            delete old[k];
        } else {
          old[k] = obj[k]
        }
    }
    var str = "#";
    for(var k in old){
        str += `${k}=${old[k]}&`
    }

    window.location.hash = str.slice(0, -1); // drop final '&'
    return old;

}

window.onhashchange = function(event) {
    let old_hp = getHashParams(new URL(event.oldURL).hash.substring(1))
    // console.log(old_hp)
    let new_hp = getHashParams()
    update_selected_gene_header(new_hp.gene)
}

function update_selected_gene_header(gene) {
    document.getElementById("header-selected-gene").innerHTML = gene
}

function update_header_metadata(n_samples, n_genes) {
    document.getElementById("sample-count").innerHTML = n_samples
    document.getElementById("gene-count").innerHTML = n_genes
}

// https://stackoverflow.com/questions/8069315/create-array-of-all-integers-between-two-numbers-inclusive-in-javascript-jquer/8069367
function range(start, end, step = 1) {
    const len = Math.floor((end - start) / step) + 1
    return Array(len).fill().map((_, idx) => start + (idx * step))
}

function get_gene_plot_layout(gene) {
    let step = Math.round(gene.plot_coords.x.length / 10)
    let layout = {
        grid: { rows: 3, columns: 1, },
        autosize: true,
        height: gene_plot_height,
        margin: { t: 10, r: 0, b: 30 },
        hovermode: 'x',
        xaxis: {
            tickvals: range(gene.plot_coords.x[0], gene.plot_coords.x[gene.plot_coords.x.length - 10], step),
            ticktext: range(gene.plot_coords.g[0], gene.plot_coords.g[gene.plot_coords.x.length - 10], step),
            showgrid: false,
            ticklen: 5,
        },
        yaxis: { title: "Depth", domain: [0.28, 1], },
        yaxis2: {
            range: [0.5, -0.5],
            showlegend: false,
            zeroline: false,
            tickvals: [0],
            ticktext: ["Merged<br>Transcript"],
            tickangle: 40,
            domain: [0.25, 0.28],
        },
        yaxis3: {
            range: [0, 2],
            zeroline: false,
            showlegend: false,
            domain: [0.0, 0.25],
            tickangle: 40,
            ticktext: gene.transcripts.map(t => t.data.transcript),
        },
        hoverdistance: 100000,
        hoverinfo: "none",
        showlegend: false,
    }

    return (layout)
}

//Get a by position trace per sample of depth
function get_depth_trace(gene) {

    var traces = [];
    var low_lvl = 10000;
    var low_dp;
    for (var slvl in gene.plot_coords.background_depths) {
        let lvl = slvl.replace(/^\D+/g, "")

        let dp = gene.plot_coords.background_depths[slvl].map(function (v) { return v < -1000 ? NaN : v })
        if (lvl < low_lvl) { low_dp = dp; low_lvl = lvl }

        var trace = {
            x: gene.plot_coords.x,
            text: gene.plot_coords.g,
            y: dp,
            type: "scatter",
            mode: "lines",
            name: `Percentile: ${lvl}`,
            tracktype: "background",
            line: { width: 1, /* dash: "dot", */ color: 'rgb(200, 200, 200)' },
            yaxis: "y",
            hoverinfo: "skip",
        }
        traces.push(trace)
    }

    var i = -1;

    for (var sample in gene.plot_coords.depths) {
        i += 1
        var dp = gene.plot_coords.depths[sample]
        dp = dp.map(function (v) { return v < -1000 ? NaN : v })
        var color = color_list[i];
        var trace = {
            x: gene.plot_coords.x, text: gene.plot_coords.g, y: dp,
            type: 'scattergl', mode: 'lines', name: sample, line: { width: 0.36, color: color_list[i] },
            hovertemplate: HOVER_TEMPLATE,
            hoverdistance: 1000,
            hoverinfo: "none",
            yaxis: "y",
            tracktype: "sample",
        };

        traces.push(trace);

        // extract parts of this sample that are below the threshold and plot
        // in a similar trace with wider line.
        if (low_dp) {
            var low_trace = {
                x: [], y: [], text: [], type: 'scatter', mode: "lines", name: sample, line: { width: 3.0, color: color_list[i] },
                hoverinfo: "none",
                connectgaps: false,
                yaxis: "y",
            };

            dp.forEach((d, i) => {
                if (d < low_dp[i]) {
                    low_trace.x.push(gene.plot_coords.x[i])
                    low_trace.y.push(d)
                } else {
                    if (low_trace.x.length > 0 && !isNaN(low_trace.x[low_trace.x.length - 1])) {
                        low_trace.x.push(gene.plot_coords.x[i])
                        low_trace.y.push(null)
                    }
                }

            })
            traces.push(low_trace)

        }
    };

    return (traces)

};

//Get the transcript shapes for the current region
function get_transcript_traces(gene, transcripts, plotRef) {
    var traces = []
    // use negative values for y_offset so we can tell from y-value if we
    // are in depth plot (which must be positive) or in transcript plot

    transcripts.forEach((transcript, y_offset) => {
        var trs = transcript.traces(-(y_offset), gene.plot_coords.x, gene.plot_coords.g)
        trs.forEach(t => {
            t.yaxis = plotRef
            traces.push(t)
        })
    })
    return traces
};


function plot_per_base_depth(gene) {
    setHashParams({'gene': gene.symbol})
    let gene_layout = get_gene_plot_layout(gene)
    let depth_traces = get_depth_trace(gene)
    let unioned_transcript_traces = get_transcript_traces(gene, [gene.unioned_transcript], "y2")
    unioned_transcript_traces.forEach(t => depth_traces.push(t))

    let transcript_traces = get_transcript_traces(gene, gene.transcripts, "y3")
    transcript_traces.forEach(t => depth_traces.push(t))
    // each transcript is centered on the negative integers.
    gene_layout.yaxis3.range = [0.5, -(1 + transcript_traces.length / 3)]
    gene_layout.yaxis3.tickvals = [...Array(gene.transcripts.length).keys()].map(f => -f)
    gene_layout.yaxis3.ticktext = gene.transcripts.map(t => t.data.transcript)

    let d = document.getElementById("gene_plot")
    let p = Plotly.newPlot(d, depth_traces, gene_layout)
    d.on("plotly_hover", data => {
        // handle_hover(data, depth_traces, gene, gene_layout)
        Plotly.react("gene_plot", depth_traces, gene_layout)
        data.event.stopPropagation()
    }).on("plotly_unhover", data => {
        if (gene_layout.shapes != undefined && gene_layout.shapes.length > 0) {
            gene_layout.shapes.pop()
        }
        Plotly.react("gene_plot", depth_traces, gene_layout)
        data.event.stopPropagation()
    })
    return p
}

function mean(data) {
    var result = 0
    var n = 0
    for (var i = 0; i < data.length; i++) {
        let d = data[i]
        if (d >= 0) {
            n += 1
            result += d
        }
    }
    return result / n
}

function handle_hover(data, depth_traces, gene, gene_layout) {
    //  this function handles hover events in the transcript plots. it
    //  finds the correct transcript and feature (exon/transcript/UTR, etc)
    //  and calculates stats, and updates the gene_layout to draw the
    //  region shape.
    //
    let ax = data.yaxes[0]._attr;
    // don't handle hover in depth plot
    if (ax == "yaxis") { return false; }
    // this is the x in plot coordinates.
    var x = Math.round(data.xvals[0])
    var y = Math.round(Math.abs(data.yvals[0])) // can now use this as an index to get the transcript.
    var transcript = ax == "yaxis2" ? gene.unioned_transcript : gene.transcripts[y]


    if (transcript == undefined) { return false }

    let overlaps = transcript.overlaps(x)
    if (overlaps.length > 1) {
        if (overlaps.some(p => p.type == FeatureType.CDS)) {
            overlaps = overlaps.filter(p => p.type == FeatureType.CDS)
        } else if (overlaps.some(p => p.type == FeatureType.EXON)) {
            overlaps = overlaps.filter(p => p.type == FeatureType.EXON)
        }
    }
    if (overlaps.length == 0) {
        overlaps.push(new Feature(transcript.txstart, transcript.txend, FeatureType.TRANSCRIPT))
    }
    if (gene_layout.shapes == undefined) { gene_layout.shapes = [] }
    gene_layout.shapes.push({
        type: "rect",
        xref: "x",
        yref: "paper",
        y0: 0,
        y1: 1,
        fillcolor: "#cccccc",
        opacity: 0.2,
        x0: overlaps[0].start,
        x1: overlaps[0].stop,
    })

    var start_idx = binary_search(gene.plot_coords.x, overlaps[0].start)
    var stop_idx = binary_search(gene.plot_coords.x, overlaps[0].stop)

    let low_depth_cutoff = 7; // TODO: get this from user-form input.
    var selection_stats = transcript.stats([{ start: start_idx, stop: stop_idx }], gene.plot_coords.depths, gene.plot_coords.background_depths, low_depth_cutoff)
    // console.log(selection_stats)

    var tx_stats = transcript.stats([{ start: 0, stop: gene.plot_coords.x.length }], gene.plot_coords.depths, gene.plot_coords.background_depths, low_depth_cutoff)

    // example of how to get stats for all CDSs. NOTE: this should only be done
    // once per gene and result cached.
    var cds = transcript.parts().filter(p => p.type == FeatureType.CDS)
    var cds_stats = transcript.stats(cds, gene.plot_coords.depths, gene.plot_coords.background_depths, low_depth_cutoff)
    // console.log("CDS:", cds_stats)

    var gstart = gene.plot_coords.g[start_idx]
    var gstop = gene.plot_coords.g[stop_idx]

    var columns = [{ title: "sample" }, { title: "transcript mean" }, { title: "CDS mean" }, { title: "selection mean" },
    { title: "transcript bases < lower" }, { title: "CDS bases < lower" }, { title: "selection bases < lower" }];
    var table = []

    //var means = {}
    //hoverInfo.innerHTML = `${gene.unioned_transcript.chr}:${gstart}-${gstop}<br><ul>`
    for (var sample in gene.plot_coords.depths) {
        var row = [sample, tx_stats[sample].mean.toFixed(2), cds_stats[sample].mean.toFixed(2), selection_stats[sample].mean.toFixed(2), tx_stats[sample].low, cds_stats[sample].low, selection_stats[sample].low]
        table.push(row)
        //let depths = gene.plot_coords.depths[sample];
        //means[sample] = mean(depths.slice(start_idx, stop_idx)).toPrecision(4)
    }
    if (datatable) {
        datatable.clear()
        datatable.rows.add(table)
        datatable.draw()
    } else {
        datatable = jQuery('table#stats_table').DataTable({ data: table, columns: columns })
    }
}

function highlight_sample(sample) {
    setHashParams({'sample': sample})
    let d = document.getElementById("gene_plot")
    let vals = d.data.map((t, i) => {
        if (t.tracktype == "background") {
            return [1, 1.5]
        }
        if (t.tracktype != 'sample') {
            return [undefined, undefined]
        } else {
            if (t.name == sample) {
                return [1, 1.5]
            } else {
                return [0.25, 0.36]
            }
        }
    })

    Plotly.restyle(d, {'opacity': vals.map(i => i[0]), 'line.width': vals.map(i => i[1]), 'hovertemplate': vals.map(i => i[1] > 0.8 ? HOVER_TEMPLATE: null)})
}

function tie_heatmap_to_line_plot() {
    // clicks on heatmap x-axis labels highlights sample in line plot
    // https://issue.life/questions/44297012
    var xticks = jQuery("#heatmap_plot .xaxislayer-above > .xtick > text")
    xticks.css({ "cursor": "pointer" })
    var d = document.getElementById("gene_plot")

    xticks.each(function (i) {
        var item = jQuery(this);
        item.css({ "fill": color_list[i] })
        item.attr('pointer-events', 'all')
        item.on("click", function (e) {
            var n = jQuery(this)
            var undo = n.css("font-weight") == "800";
            xticks.attr('class', 'unselected_label')
            var sample = n.text()
            var vals = null
            if (undo) {
                vals = d.data.map((t, i) => [1, (t.tracktype == 'sample' || t.tracktype == 'background') ? 0.36 : undefined])
                Plotly.restyle(d, {'opacity': vals.map(i => i[0]), 'line.width': vals.map(i => i[1]), 'hovertemplate': d.data.map((t) => t.tracktype == "sample" ? HOVER_TEMPLATE: null)})
                setHashParams({'sample': null})
            } else {
                n.attr('class', 'selected_label')
                highlight_sample(sample)
            }
            // see: https://codepen.io/etpinard/pen/RQQqzq?editors=0010
        }).on("unhover", function (e) {
        })
    })
}

function draw_heatmap() {

    let stat_metric = jQuery("#metric_select").val()

    var low_depth_cutoff = 7;
    var z = []
    var x = [] // samples
    var y = [] // genes

    var is_cds = stat_metric.startsWith("cds")
    var metric = stat_metric.split("_")[1]
    for (var i = 0; i < plot_data.length; i++) {
        var g = plot_data[i]
        y.push(g.symbol)

        var transcript = g.unioned_transcript
        var stats = null
        if (is_cds) {
            // console.log("CDS")
            var cds = transcript.parts().filter(p => p.type == FeatureType.CDS)
            stats = transcript.stats(cds, g.plot_coords.depths, g.plot_coords.background_depths, low_depth_cutoff)
        } else {
            // console.log("transcript")
            stats = transcript.stats([{ start: 0, stop: g.plot_coords.x.length }], g.plot_coords.depths, g.plot_coords.background_depths, low_depth_cutoff)
        }

        var row = []

        for (var s in stats) {
            if (i == 0) {
                x.push(s)
            }
            row.push(stats[s][metric])
        }
        // console.log(row)
        z.push(row)
    }
    // heatmap draws from bottom up by default
    z = z.reverse()
    y = y.reverse()

    let hlayout = {
        title: "",
        margin: { t: 10 },
        height: 20 * y.length + 60,
        paper_bgcolor: '#f8f9fa',
        xaxis: { fixedrange: true, },
        yaxis: { fixedrange: true, },
        font: { size: 14 },
    }
    let heatmap_config = {
        displaylogo: false,
        toImageButtonOptions: {
            format: 'svg', // one of png, svg, jpeg, webp
            filename: 'seqcover_summary',
            height: 18 * y.length + 60,
            width: 1200,
            scale: 1 // Multiply title/legend/axis/canvas sizes by this factor
        },
        // responsive: true,
    }
    let trace = [{ x: x, y: y, z: z, type: 'heatmap', hoverongaps: false, colorscale: "Cividis" }]
    // let trace = [{ x: x, y: y, z: z, type: 'heatmap', xgap: 1, ygap: 1, hoverongaps: false, colorscale: "Cividis" }]

    //https://plotly.com/javascript/reference/heatmap/
    let p = document.getElementById('heatmap_plot')

    let yticks = jQuery("#heatmap_plot .yaxislayer-above > .ytick > text")
    let selected = false
    yticks.each(function (i) {
        let elem = jQuery(this)
        if (elem.css("font-weight") == "800") {
            selected = i
        }
    })

    Plotly.newPlot(p, trace, hlayout, heatmap_config)
    p.removeAllListeners("plotly_click")
    p.on('plotly_click', function (click_data) {
        let selected_gene_idx = click_data.points[0].pointIndex[0]
        let gene_idx = plot_data.length - selected_gene_idx - 1
        let selected_sample_idx = click_data.points[0].pointIndex[1]
        let sample = click_data.points[0].x
        let gene = click_data.points[0].y

        // y-axis ticks
        let yticks = jQuery("#heatmap_plot .yaxislayer-above > .ytick > text")
        yticks.attr('class', 'unselected_label')
        yticks.each(function (i) {
            if (i == selected_gene_idx) {
                jQuery(this).attr('class', 'selected_label')
            }
        })
        // x-axis ticks
        let xticks = jQuery("#heatmap_plot .xaxislayer-above > .xtick > text")
        xticks.attr("class", "unselected_label")
        xticks.each(function (i) {
            if (i == selected_sample_idx) {
                jQuery(this).attr("class", "selected_label")
            }
        })

        if(selected_gene_idx != selected) {
            // only redraw if it's a new gene
            plot_per_base_depth(plot_data[gene_idx]).then(highlight_sample(sample))
            selected = selected_gene_idx;
            setHashParams({'gene': gene, 'sample': sample});
        } else {
            highlight_sample(sample)
        }

    })
    p.on('plotly_afterplot', function () {
        // initialize first gene as selected; or re-use selected (when changing metric)
        let yticks = jQuery("#heatmap_plot .yaxislayer-above > .ytick > text")
        if (!selected) {
            jQuery(yticks[plot_data.length - 1]).attr("class", "selected_label")
        } else {
            jQuery(yticks[selected]).attr("class", "selected_label")
        }

        // x-axis ticks
        tie_heatmap_to_line_plot()
        // y-axis ticks
        Plotly.d3.selectAll(".yaxislayer-above").selectAll('text')
            .on("click", function (d) {
                yticks.attr("class", "unselected_label")
                jQuery(this).attr('class', 'selected_label')
                var hash = setHashParams({'gene': d.text})
                plot_per_base_depth(plot_data[plot_data.length - d.x - 1])
                if("sample" in hash){ highlight_sample(hash.sample) }
            })
    })
}

jQuery(document).ready(function () {
    update_header_metadata(Object.keys(plot_data[0].plot_coords.depths).length, plot_data.length)

    if (Object.keys(plot_data[0].plot_coords.background_depths).length == 0) {
        jQuery('#cds_lt_bg').remove()
    }

    for(var g of plot_data) {
        if (g.unioned_transcript.constructor.name == "Object") {
            g.unioned_transcript = new Transcript(g.unioned_transcript)
            g.transcripts = g.transcripts.map(t => new Transcript(t))
        }
    }

    jQuery('#metric_select').on('change', function() { draw_heatmap() }).trigger('change')
    var p = getHashParams();
    var i = 0;
    if ("gene" in p) {
        var gene = p.gene;
        update_selected_gene_header(gene)
        for(var j = 0; j < plot_data.length; j++){
            if(plot_data[j].symbol == gene){
                i = j
                break
            }
        }

    }
    plot_per_base_depth(plot_data[i]).then(function() {
        // NOTE we should put these events directly in the plot fn so we don't
        // have code duplication, but here for now...
        let yticks = jQuery("#heatmap_plot .yaxislayer-above > .ytick > text")
        yticks.attr('class', 'unselected_label')
        yticks.each(function (k) {
            if (i == plot_data.length - k - 1) {
                jQuery(this).attr('class', 'selected_label')
            }
        })
        if("sample" in p) {
            highlight_sample(p.sample)
            jQuery("#heatmap_plot .xaxislayer-above > .xtick > text").each(function(i, v) {
                if(v.innerHTML == p.sample){
                    jQuery(v).attr('class', 'selected_label')
                }
            })

        }
    })
})
