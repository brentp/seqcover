    var nan = NaN; // hack to support json dumped with NaN values.
    const gene_plot_height = 500
    const plot_colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf"]

    var initalIndex = 0
    var currentRegion = "";
    var nonCdsChecked = false;
    var utrsChecked = false;
    var showTranscripts = true;

    /*
    ------------------------------------------------------------------------------------------------------------------
                                                      Plot Functions
    ------------------------------------------------------------------------------------------------------------------
    */

    //by base depth plot layout
    function get_gene_plot_layout(gene) {

        var mid = Math.round(gene.plot_coords.x.length / 2)

        var layout = {
            grid: {
                rows: 3,
                columns: 1,
            },
            autosize: true,
            height: gene_plot_height,
            margin: { t: 10, r: 0, b: 30 },
            xaxis: {
                tickmode: "array",
                tickvals: [gene.plot_coords.x[0], gene.plot_coords.x[mid], gene.plot_coords.x[gene.plot_coords.x.length - 10]],
                ticktext: [gene.plot_coords.g[0], gene.plot_coords.g[mid], gene.plot_coords.g[gene.plot_coords.x.length - 10]],
                // title: "Chromosome " + String(gene.unioned_transcript.chr).replace("chr", "")
            },
            yaxis: {
                title: "Depth",
                domain: [0.55, 1]
            },
            yaxis2: {
                title: "Merged<br>Transcripts",
                range: [0.5, -0.5],
                showlegend: false,
                zeroline: false,
                showticklabels: false,
                ticktext: gene.unioned_transcript.transcript,
                domain: [0.0, 0.3],

            },

            hovermode: 'x unified',
            hoverdistance: 100000,
            showlegend: false,
            legend: {
                xanchor: "right",
                yanchor: "top",
                y: 1,
                x: 1,
                orientation: "h",
                borderwidth: 1,
                bordercolor: '#eeeeee'
            },
        };

        //Add the 3rd layout if show transcripts is on
        if (showTranscripts) {
            layout.yaxis.domain = [0.55, 0.90]
            layout.yaxis2.domain = [0.4, 0.55]
             layout["yaxis3"] = {
                range: [0, 2],
                zeroline: false,
                showlegend: false,
                domain: [0.0, 0.40]
            }

        };

        return (layout)

    };


    //Get a by position trace per sample of depth
    function get_depth_trace(gene) {

        var traces = [];

        for (sample in gene.plot_coords.depths) {
            var dp = gene.plot_coords.depths[sample]
            dp = dp.map(function (v) { return v < -1000 ? NaN : v })

            var trace = {
                x: gene.plot_coords.x, text: gene.plot_coords.g, y: dp,
                type: 'scatter', mode: 'lines', name: sample, line: { width: 1 },
                hovertemplate: '<b>position</b>:%{text}<br><b>depth</b>:%{y}<br>(debug) x: %{x}',
                hoverinfo: "text",
                yaxis: "y",
            };

            traces.push(trace);
        };

        return (traces)

    };


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


    //Get a map of transcript information to use as a customdata imput for plotly hover data
    function get_custom_data(transcript, exonStart, exonEnd, exonNumber, cdsExon, feature, xs, gposition_list) {

        custdata = {
            "tx_id": transcript.transcript,
            "strand": transcript.strand > 0 ? "+" : "-",
            "tx_start": gposition_list[binary_search(xs, transcript.txstart)],
            "tx_end": gposition_list[binary_search(xs, transcript.txend)],
            "cds_start": gposition_list[binary_search(xs, transcript.cdsstart)],
            "cds_end": gposition_list[binary_search(xs, transcript.cdsend)],
            "exon_start": gposition_list[binary_search(xs, exonStart)],
            "exon_end": gposition_list[binary_search(xs, exonEnd)],
            "exon_count": exonNumber,
            "CDS": cdsExon,
            "Feature": feature
        }

        return (custdata)

    };


    //Get UTR information/shape for a transcript
    function getUTRInfo(gene, symbol, transcript, utrOffset, yValue, customData, shapes, x_values, y_values, text_symbol, plotRef) {

        //Get Start UTRs
        let utr1 = {
            type: "rect", x0: transcript.txstart, x1: transcript.cdsstart,
            y0: yValue - utrOffset, y1: yValue + utrOffset,
            line_color: "black", fillcolor: "cyan", yref: plotRef
        }


        //Get End UTR
        let utr2 = {
            type: "rect", x0: transcript.cdsend, x1: transcript.txend,
            y0: yValue - utrOffset, y1: yValue + utrOffset,
            line_color: "cyan", fillcolor: "cyan", yref: plotRef
        }


        // If CDS region is at the end of the transcript, the whole region is a UTR
        if (transcript.txend <= transcript.cdsstart) {

            //Add UTR info to customData
            ///UTR1
            customData.push(get_custom_data(transcript, transcript.txstart,
                transcript.cdsstart, "UTR", "UTR", "UTR", gene.plot_coords.x, gene.plot_coords.g))
            text_symbol.push(symbol)

            //Add UTRs to shapes
            ///UR1
            shapes.push(utr1)

            // Add UTR plot coordinates for hover info
            ///UTR1
            x_values.push((transcript.txstart + ((transcript.cdsstart - transcript.txstart) / 2)))
            y_values.push(yValue)

            // IF CDS regions ends at the beginning of the transcript, the whole region is a UTR
        } else if (transcript.txstart >= transcript.cdsend) {

            //Add UTR info to customData
            ///UTR2
            customData.push(get_custom_data(transcript, transcript.cdsend,
                transcript.txend, "UTR", "UTR", "UTR", gene.plot_coords.x, gene.plot_coords.g))
            text_symbol.push(symbol)

            //Add UTRs to shapes
            ///UTR2
            shapes.push(utr2)

            // Add UTR plot coordinates for hover info
            ///UTR2
            x_values.push((transcript.cdsend + ((transcript.txend - transcript.cdsend) / 2)))
            y_values.push(yValue)

        } else {

            //Add UTR info to customData
            ///UTR1
            customData.push(get_custom_data(transcript, transcript.txstart,
                transcript.cdsstart, transcript.strand > 0 ? "5'-UTR" : "3'-UTR",
                "UTR", "UTR", gene.plot_coords.x, gene.plot_coords.g))
            text_symbol.push(symbol)
            ///UTR2
            customData.push(get_custom_data(transcript, transcript.cdsend,
                transcript.txend, transcript.strand > 0 ? "3'-UTR" : "5'-UTR",
                "UTR", "UTR", gene.plot_coords.x, gene.plot_coords.g))
            text_symbol.push(symbol)

            //Add UTRs to shapes
            ///UR1
            shapes.push(utr1)
            ///UTR2
            shapes.push(utr2)

            // Add UTR plot coordinates for hover info
            ///UTR1
            x_values.push((transcript.txstart + ((transcript.cdsstart - transcript.txstart) / 2)))
            y_values.push(yValue)
            ///UTR2
            x_values.push((transcript.cdsend + ((transcript.txend - transcript.cdsend) / 2)))
            y_values.push(yValue)

        };

        return ({ "custData": customData, "shapeData": shapes, "xValues": x_values, "yValues": y_values, "symbolArray": text_symbol })
    };


    //Get CDS Exons for a transcript
    function getCdsExons(gene, symbol, transcript, exon_2d_array, cdsExonCount, exonOffset, yValue, customData, shapes, x_values, y_values, text_symbol, plotRef) {

        // Get CDS exon numbers
        let exonNumbers = transcript.strand > 0 ? pv.range(1, cdsExonCount + 1) : pv.range(cdsExonCount, 0, -1)
        let exonNumberOffset = 0

        //Add exon shapes
        for (index in exon_2d_array) {

            let exonPos = exon_2d_array[index]

            //If exon in CDS region
            if (exonPos[1] > transcript.cdsstart & exonPos[0] < transcript.cdsend) {

                //Get Exon Positions
                let exonStart = exonPos[0] > transcript.cdsstart ? exonPos[0] : transcript.cdsstart
                let exonEnd = exonPos[1] < transcript.cdsend ? exonPos[1] : transcript.cdsend

                //Add Shapes
                shapes.push({
                    type: "rect", x0: exonStart, x1: exonEnd, y0: yValue - exonOffset,
                    y1: yValue + exonOffset, line_color: "black", fillcolor: "black", yref: plotRef
                })

                //Get exon region info
                let feature = "Coding Exon"
                let codingRegion = "True"
                if (exonStart != exonPos[0] | exonEnd != exonPos[1]) {
                    feature = "Partial Coding Exon"
                    codingRegion = "Partial"
                };

                //Add exon info to customData
                /// Add the full exon start and end position, even if the exon is partially in non-CDS region
                customData.push(get_custom_data(transcript, exonPos[0], exonPos[1], exonNumbers[(parseInt(index) - exonNumberOffset)],
                    codingRegion, feature, gene.plot_coords.x, gene.plot_coords.g))
                text_symbol.push(symbol)

                // Add exon plot coordinates for hover info
                //If non-CDS exons are being shown, the hover x coordinate will change
                if (nonCdsChecked) {
                    //Check which coordinate, if any, are partial Non-CDS
                    if (exonStart != exonPos[0]) {
                        x_values.push(exonStart)

                    } else if (exonEnd != exonPos[1]) {
                        x_values.push(exonEnd)

                    } else {
                        x_values.push((exonStart + ((exonEnd - exonStart) / 2)))

                    };

                } else {
                    x_values.push((exonStart + ((exonEnd - exonStart) / 2)))

                };

                y_values.push(yValue)

                // If exon not in CDS region
            } else {
                exonNumberOffset += 1
            };
        };

        return ({ "custData": customData, "shapeData": shapes, "xValues": x_values, "yValues": y_values, "symbolArray": text_symbol })
    };


    //Get Non-CDS Exons for a transcript
    function getNonCdsExons(gene, symbol, transcript, exon_2d_array, exonOffset, yValue, customData, shapes, x_values, y_values, text_symbol, plotRef) {

        //Add exon shapes
        for (index in exon_2d_array) {

            let exonPos = exon_2d_array[index]

            //Skip any full CDS exon
            if (exonPos[0] >= transcript.cdsstart & exonPos[1] <= transcript.cdsend) { continue };

            //Get Exon Positions
            let exonStart = exonPos[0]
            let exonEnd = exonPos[1]
            let partialCds = false

            //Adjust start and end position if they overlap CDS region
            /// IF exon start is before cds start and exon end is after, adjust end position
            if (exonPos[0] < transcript.cdsstart & exonPos[1] > transcript.cdsstart) {
                exonEnd = transcript.cdsstart
                partialCds = true
            };

            /// IF exon start is before cds end and exon end is after, adjust start position
            if (exonPos[0] < transcript.cdsend & exonPos[1] > transcript.cdsend) {
                exonStart = transcript.cdsend
                partialCds = true
            };

            //Add Shapes
            shapes.push({
                type: "rect", x0: exonStart, x1: exonEnd, y0: yValue - exonOffset,
                y1: yValue + exonOffset, opacity: 0.4, line_color: "black", fillcolor: "grey", yref: plotRef
            })


            //Add exon info to customData
            /// Add the full exon start and end position, even if the exon is partially in non-CDS region
            if (!(partialCds)) {
                customData.push(get_custom_data(transcript, exonStart, exonEnd, "NA", "False", "Non-Coding Exon",
                    gene.plot_coords.x, gene.plot_coords.g))
                text_symbol.push(symbol)

                x_values.push((exonStart + ((exonEnd - exonStart) / 2)))
                y_values.push(yValue)
            };

        };

        return ({ "custData": customData, "shapeData": shapes, "xValues": x_values, "yValues": y_values, "symbolArray": text_symbol })
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

        gene_layout = get_gene_plot_layout(gene)

        var depth_traces = get_depth_trace(gene)

        var unioned_transcript_traces = get_transcript_traces(gene, [gene.unioned_transcript], "y2")
        unioned_transcript_traces.forEach(t => depth_traces.push(t))

        if (showTranscripts) {
          var transcript_traces = get_transcript_traces(gene, gene.transcripts, "y3")
          transcript_traces.forEach(t => depth_traces.push(t))
          gene_layout.yaxis3.range = [0.5, -(1 + transcript_traces.length / 3)]

        };

        Plotly.newPlot("gene_plot", depth_traces, gene_layout)

        var d = document.getElementById("gene_plot")
        var hoverInfo = document.getElementById("hoverInfo")
        d.on("plotly_hover", data => {
            handle_hover(data, depth_traces, gene)
        }).on("plotly_unhover", data => hoverInfo.innerHTML = "");

    };

    function handle_hover(data, depth_traces, gene) {
        // don't handle hover in depth plot
        let ax = data.yaxes[0]._attr;
        if(ax == "yaxis") {return; }
        var x = Math.round(data.xvals[0])
        var y = Math.round(Math.abs(data.yvals[0])) // can now use this as an index to get the transcript.
        var transcript = ax == "yaxis2" ? gene.unioned_transcript : gene.transcripts[y]
        if(transcript == undefined) {return}
        console.log(x, transcript)

        let overlaps = transcript.parts().filter(p => {
            //console.log(x, p, p.start <= x, p.stop >= x)
            return p.start <= x && p.stop >= x
        })
        if(overlaps.length > 1) {
            overlaps = overlaps.filter(p => p.type == FeatureType.CDS)
        }
        if(overlaps.length == 0){
            overlaps.push(new Feature(transcript.txstart, transcript.txstop))
        }

        var start_idx = binary_search(gene.plot_coords.x, overlaps[0].start)
        var stop_idx = binary_search(gene.plot_coords.x, overlaps[0].stop)

        var means = {}
        hoverInfo.innerHTML = ""
        for(var sample in gene.plot_coords.depths) {
            let depths = gene.plot_coords.depths[sample];
            means[sample] = math.mean(depths.slice(start_idx, stop_idx))
            hoverInfo.innerHTML += `<b>${sample}</b> mean depth for ${overlaps[0].FeatureType}: ${means[sample]}<br>`
        }


    }


    /*
    ------------------------------------------------------------------------------------------------------------------
                                                     Event Handling
    ------------------------------------------------------------------------------------------------------------------
    */

    // Controller function for generating gene/region specific plots
    function generate_plots(selected_region) {

        var g = plot_data[selected_region]
        if(g.unioned_transcript.constructor.name == "Object") {
            g.unioned_transcript = new Transcript(g.unioned_transcript)
        }
        g.transcripts = g.transcripts.map(t => new Transcript(t))

        //Plot per base depths
        plot_per_base_depth(g)

    };


    //Load first region
    jQuery(document).ready(function () {
        let regions = []
        for (const [i, pd] of plot_data.entries()) {
            console.log('%d: %s', i, pd)
            let region = `${pd.symbol} ${pd.unioned_transcript.chr}:${pd.plot_coords.g[0]}-${pd.plot_coords.g[pd.plot_coords.g.length - 1]}`
            regions.push({ "id": i, "text": region })
        }
        generate_plots(0)
        // n = gene symbol, v = the index
        $("#region-select").select2({
            data: regions,
            selectOnClose: true,
            width: 'resolve',
            theme: 'bootstrap4',
            width: $(this).data('width') ? $(this).data('width') : $(this).hasClass('w-100') ? '100%' : 'style'
        })

        // register handler
        $('#region-select').on('change', function () {
            // untested (https://select2.org/programmatic-control/retrieving-selections)
            generate_plots(this.value)
        })

    })

    //Turn on or off UTR regions in the transcript plot
    $("#utrs").on('change', function (e) {

        utrsChecked = this.checked

        //Redraw per base depth plot
        //Index of plot in plot_data and gene data
        var plotIndex = parseInt($('#region-select').find(':selected')[0].value)
        var selectedGene = plot_data[plotIndex];

        //Plot per base depths
        plot_per_base_depth(selectedGene)

    });


    //Turn on or off Non-CDS Exons in the transcript plot
    $("#nonCdsExons").on('change', function (e) {

        nonCdsChecked = this.checked

        //Redraw per base depth plot
        //Index of plot in plot_data and gene data
        var plotIndex = parseInt($('#region-select').find(':selected')[0].value)
        var selectedGene = plot_data[plotIndex];

        //Plot per base depths
        plot_per_base_depth(selectedGene)

    });

    //Turn on and off the transcript plot
    $("#showTranscripts").click(function () {

        //Update bool value for transcripts
        showTranscripts = !(showTranscripts)

        // Update button text
        // Update div class
        if (showTranscripts) {
            $(this).text("Hide Transcripts")
            document.getElementById("gene_plot").className = "big_div";
        } else {
            $(this).text("Show Transcripts")
            document.getElementById("gene_plot").className = "med_div";
        };

        //Redraw per base depth plot
        //Index of plot in plot_data and gene data
        var plotIndex = parseInt($('#region-select').find(':selected')[0].value)
        var selectedGene = plot_data[plotIndex];

        //Plot per base depths
        plot_per_base_depth(selectedGene)

    });
