    var nan = NaN; // hack to support json dumped with NaN values.
    const gene_plot_height = 500
    const plot_colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf"]

    var initalIndex = 0
    var currentRegion = "";
    var nonCdsChecked = false;
    var utrsChecked = false;
    var showTranscripts = false;

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
                range: [0, 2],
                showlegend: false,
                zeroline: false,
                showticklabels: false,
                ticktext: gene.unioned_transcript.transcript,
                domain: [0.0, 0.3],

            },

            hovermode: 'closest',
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

            //Update yaxis domains
            layout.yaxis.domain = [0.55, 0.90]
            layout.yaxis2.domain = [0.4, 0.55]

            //Add 3rd axis
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
    function get_by_position_depth_trace(gene) {

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
    function get_transcript_shapes(gene, transcripts, plotRef) {

        let symbol = gene.symbol
        let top = 2;
        let bottom = 0;
        let exonOffset = 0.25;
        let utrOffset = 0.1;
        let shapes = [];
        let tick_labels = [];
        let tick_vals = [];
        let strandArray = [];
        let x_values = []
        let y_values = []
        let customData = []
        let text_symbol = []

        for (transcript of transcripts) {

            //Sort exon positions
            exon_2d_array = transcript.position.sort(function (a, b) { return a[0] - b[0] })

            //Count cds exons
            let cdsExonCount = 0
            for (index in exon_2d_array) {

                let exonPos = transcript.position[index]

                if (exonPos[1] > transcript.cdsstart & exonPos[0] < transcript.cdsend) { cdsExonCount += 1 };

            };

            //If neither UTRs or non-CDS Exons are being shown, and the CDS Exon Count is 0, skip this transcript.
            if (utrsChecked | nonCdsChecked | cdsExonCount > 0) {

                // Transcript
                shapes.push({
                    type: "line", x0: transcript.txstart, x1: transcript.txend,
                    y0: bottom + 1, y1: bottom + 1, line_color: "black", yref: plotRef
                })

                //Get strand annotation info
                let max_index = transcript.position.length - 1
                let strandMarker = transcript.strand > 0 ? ">" : "<";
                let markerOffset = showTranscripts ? 0.01 : 0.02;
                for (index in transcript.position) {

                    //Get strand annotation text
                    if (parseInt(index) == max_index) { break };

                    //Get start and end exon positions
                    let startPos = transcript.position[parseInt(index)][0]
                    let endPos = transcript.position[parseInt(index)][1]
                    let nextStartPos = transcript.position[(parseInt(index) + 1)][0]
                    let xPos = endPos + ((nextStartPos - endPos) / 2)

                    //Only add strand marker annotation to CDS region
                    if (endPos > transcript.cdsstart & nextStartPos < transcript.cdsend) {
                        //Add marker to strandArray
                        strandArray.push({ x: xPos, y: bottom + (1 + markerOffset), text: strandMarker, font: { color: "black", size: 16 }, opacity: 0.7, showarrow: false, yref: plotRef })
                    }
                };

            } else { continue };

            //Get the CDS Exons
            let cdsExonMap = getCdsExons(gene, symbol, transcript, exon_2d_array, cdsExonCount, exonOffset, (bottom + 1), customData, shapes, x_values, y_values, text_symbol, plotRef)

            //Update data objects with CDS Exon info
            customData = cdsExonMap.custData
            shapes = cdsExonMap.shapeData
            x_values = cdsExonMap.xValues
            y_values = cdsExonMap.yValues
            text_symbol = cdsExonMap.symbolArray

            //If UTR checkbox is true, add UTRs
            if (utrsChecked) {

                //Get UTR shapes
                let utrMap = getUTRInfo(gene, symbol, transcript, utrOffset, (bottom + 1), customData, shapes, x_values, y_values, text_symbol, plotRef)

                //Update data objects with UTR info
                customData = utrMap.custData
                shapes = utrMap.shapeData
                x_values = utrMap.xValues
                y_values = utrMap.yValues
                text_symbol = utrMap.symbolArray

            };

            //Add non CDS exons if checked
            if (nonCdsChecked) {

                //Get Non-CDS Exon Shapes
                let nonCdsExonMap = getNonCdsExons(gene, symbol, transcript, exon_2d_array, exonOffset, (bottom + 1), customData, shapes, x_values, y_values, text_symbol, plotRef)

                //Update data objects with CDS Exon info
                customData = nonCdsExonMap.custData
                shapes = nonCdsExonMap.shapeData
                x_values = nonCdsExonMap.xValues
                y_values = nonCdsExonMap.yValues
                text_symbol = nonCdsExonMap.symbolArray
            };

            //Get tick labels
            tick_labels.push(transcript.transcript)
            tick_vals.push(bottom + 1)

            //update bottom
            bottom = bottom - 1
        };

        //Add a trace for the hover info
        let trace = {
            x: x_values, y: y_values, text: text_symbol,
            type: 'scatter', mode: 'markers', customdata: customData, name: "Transcripts",
            opacity: 0.0, showlegend: false, yaxis: plotRef, marker: { color: "darkcyan", size: 15 },
            hovertemplate: '<b>Gene</b>: %{text}<br><b>Transcript ID</b>: %{customdata.tx_id}<br><b>Strand</b>: %{customdata.strand}<br><b>Feature</b>: %{customdata.Feature}<br><b>Start Position</b>: %{customdata.exon_start}<br><b>End Position</b>: %{customdata.exon_end}<br><b>Coding Region</b>: %{customdata.CDS}<br><b>CDS Exon Number</b>: %{customdata.exon_count}<br>',
        };

        return ({ "shapes": shapes, "top": top, "bottom": bottom, "tick_labels": tick_labels, "tick_vals": tick_vals, "strand_annotation": strandArray, "trace": trace })
    };


    //Get descriptive stats info for stats table
    function get_stats_table(gene, low_depth = 10) {

        var dataRows = []
        var columnNames = [
            { "title": "Sample" },
            { "title": "Mean Depth" },
            { "title": "Median Depth" },
            { "title": "STD" },
            { "title": "Min Depth" },
            { "title": "Max Depth" },
            { "title": "Low Coverage Sites (<=" + low_depth + ")" }
        ]

        for (sample in gene.plot_coords.depths) {

            // Depths
            let dp = gene.plot_coords.depths[sample]
            dp = dp.filter(value => value > -1000);

            //Get depth descriptive statistics
            let dp_min = Math.min(...dp)
            let dp_max = Math.max(...dp)
            let dp_mean = math.round(math.mean(dp), 3)
            let dp_median = math.median(dp)
            let dp_std = math.round(math.std(dp), 3)
            let dp_low = math.sum(dp.map(function (v) { return v <= low_depth ? 1 : 0 }))

            dataRows.push([sample, dp_mean, dp_median, dp_std, dp_min, dp_max, dp_low])
        };
        return ({ rowData: dataRows, colNames: columnNames })
    };


    //Get the traces for each sample depth distribution
    function get_depth_distribution_traces(gene) {

        let traces = []
        let group_min_dp = 1000
        let group_max_dp = -1

        // Iterate through each samples depths
        for (sample in gene.plot_coords.depths) {

            // Sample depths
            let dp = gene.plot_coords.depths[sample]

            //Get local min and max depths
            let dp_min = 1000
            let dp_max = -1
            for (depth of dp) {

                //Skip null depths
                if (depth < 0) { continue };

                //Update local min and max
                if (depth < dp_min) { dp_min = depth }
                if (depth > dp_max) { dp_max = depth }

            };

            //Update group min and max
            if (dp_min < group_min_dp) { group_min_dp = dp_min }
            if (dp_max > group_max_dp) { group_max_dp = dp_max }

            // X and Y values
            let x_values = pv.range(dp_min, dp_max + 1)
            let y_values = new Array(x_values.length).fill(0)

            // Add y values based on depth value
            for (depth of dp) {

                //Skip null depths
                if (depth < 0) { continue };

                //Increment depth at depth index
                y_values[(depth - dp_min)] += 1
            };

            //Create trace
            let trace = {
                x: x_values, y: y_values, mode: "lines",
                name: sample, line: { width: 1 },
                hovertemplate: "<b>Depth</b>:%{x}<br><b>Count</b>:%{y}"
            }

            traces.push(trace)
        };

        return ({ "traces": traces, "min": group_min_dp, "max": group_max_dp })

    };


    //Get the depth distribution layout
    function get_depth_distribution_layout(xaxis_min = 0, xaxis_max = 100) {

        var layout = {
            autosize: true,
            title: "Depth Distribution per Sample",
            xaxis: {
                title: "Depth",
                range: [xaxis_min, xaxis_max]
            },
            yaxis: {
                title: "Count",
                domain: [0.0, 0.9],
            },
            hovermode: 'closest',
            showlegend: true,
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
        return (layout)
    };


    /*
    ------------------------------------------------------------------------------------------------------------------
                                                     Plot Controllers
    ------------------------------------------------------------------------------------------------------------------
    */

    //Per base depth plot with transcript annotations
    function plot_per_base_depth(gene) {

        //Get gene specific layout
        plot_gene_layout = get_gene_plot_layout(gene)

        //Get gene specific depth traces
        by_position_depth_traces = get_by_position_depth_trace(gene)

        // get unionized transcript shapes
        //unioned_transcript_shapes = get_unioned_transcript_shapes(gene)
        unioned_transcript_shapes = get_transcript_shapes(gene, [gene.unioned_transcript], "y2")
        by_position_depth_traces.push(unioned_transcript_shapes.trace)

        //If show transcripts is turned on, show the transcript plot
        if (showTranscripts) {
            // Get per transcript annotation shapes
            transcript_map = get_transcript_shapes(gene, gene.transcripts, "y3")
            by_position_depth_traces.push(transcript_map.trace)

            //Add shapes to layout
            plot_gene_layout.shapes = unioned_transcript_shapes.shapes.concat(transcript_map.shapes)

            //Adjuste subplot 3 yaxis ranges based on the number of annotated  transcripts
            plot_gene_layout.yaxis3.range = [transcript_map.bottom, transcript_map.top]

            //Add subplot 3 yaxis tick labels (transcript ids)
            plot_gene_layout.yaxis3.tickvals = transcript_map.tick_vals
            plot_gene_layout.yaxis3.ticktext = transcript_map.tick_labels

            //Add subplot 3 strand annotations
            plot_gene_layout.annotations = transcript_map.strand_annotation.concat(unioned_transcript_shapes.strand_annotation)


        } else {

            //Add shapes to layout
            plot_gene_layout.shapes = unioned_transcript_shapes.shapes
            plot_gene_layout.annotations = unioned_transcript_shapes.strand_annotation

        };

        //Plot by base depth
        Plotly.newPlot("gene_plot", by_position_depth_traces, plot_gene_layout)

    };


    //Add a jQuery data table with descriptive stats
    function plot_stats_table(gene) {

        //Remove previous data table if it exists
        if ($.fn.DataTable.isDataTable('#stats_table')) {
            let table = $('#stats_table').DataTable()
            table.destroy();
        };

        //Get row data and column names
        rows_and_columns = get_stats_table(gene)

        //Create data table
        jQuery("#stats_table").DataTable({
            destory: true,
            data: rows_and_columns["rowData"],
            columns: rows_and_columns["colNames"],
            paging: true,
            scroller: true,
            scrollY: '400',
            scrollX: true,
            scrollCollapse: true,
        });

    };


    //Depth Distribution plot using per depth count
    function plot_depth_distribution(gene) {

        //Get depth distribution plot trace per sample
        var depth_map = get_depth_distribution_traces(gene)

        //Get the plot layout for the depth plot
        var depth_dist_layout = get_depth_distribution_layout(depth_map.min - 5, depth_map.max + 5)

        Plotly.newPlot("depth_dist_plot", depth_map.traces, depth_dist_layout)

    };



    /*
    ------------------------------------------------------------------------------------------------------------------
                                                     Event Handling
    ------------------------------------------------------------------------------------------------------------------
    */

    // Controller function for generating gene/region specific plots
    function generate_plots(selected_region) {

        let selected_gene = plot_data[selected_region]

        //Plot per base depths
        plot_per_base_depth(selected_gene)

        //Add a descriptive statistics table
        plot_stats_table(selected_gene)

        //Plot depth distribution distrubtion
        plot_depth_distribution(selected_gene)

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
