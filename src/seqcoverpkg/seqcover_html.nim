import strutils
import ./utils
import json

const TOKEN = "SEQCOVER_JAVASCRIPT"

const html = staticRead("./seqcover.html")
const seqcover_js = staticRead("./seqcover.js")
const transcript_js = staticRead("./transcript.js")


proc write_html*(path:string, plot_data:seq[GenePlotData]) =


  var pieces = html.split(TOKEN)

  pieces = @[pieces[0], "<script>", transcript_js, seqcover_js, "XX", "</script>", pieces[1]]

  var fh:File
  if not fh.open(path, fmWrite):
    quit "could not open file:" & path
  discard

  #pieces[3] = $(%(plot_data))
  for p in pieces:
    if p == "XX":
      fh.write("var plot_data = ")
      fh.write_line($(%(plot_data)))
    else:
      fh.write_line(p)

  fh.close

