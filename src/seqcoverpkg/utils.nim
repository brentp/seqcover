import os
import d4
import sets
import httpclient
import strformat
import sequtils
import hts/fai
import algorithm
import json
import strutils
import ./typeenum
import ./transcript
export transcript

proc get_gene(gene:string, species:string="human") =
  var C = newHttpClient()
  var r = C.getContent(&"http://mygene.info/v3/query?q={gene}&species={species}")
  let js = parseJson(r)
  var gene_id = $js["hits"][0]["_id"]
  gene_id = gene_id.strip(chars={'"'})

  echo &"http://mygene.info/v3/gene/{gene_id}?fields=name,symbol,exons"
  var exons = C.getContent(&"http://mygene.info/v3/gene/{gene_id}?fields=name,symbol,exons")
  echo exons

proc `$`(j:JsonNode): string =
  return json.`$`(j).strip(chars={'"'})

proc drop_alt_chroms(g:var Gene) =
  # prefer, e.g. 22 over 22_KI270879v1_alt
  var chromset = initHashSet[string]()
  for t in g.transcripts:
    chromset.incl($(t.`chr`))

  if chromset.len <= 1: return
  var chroms: seq[string] = toSeq(chromset)
  chroms.sort do (a:string, b:string) -> int:
    result = len(a) - len(b)
  var drop: seq[int]
  for i, t in g.transcripts:
    if t.`chr` != chroms[0]: drop.add(i)

  for i in reversed(drop):
    g.transcripts.delete(i)


proc get_genes*(genes:seq[string], species:string="human", hg19:bool=false): seq[Gene] =
  var C = newHttpClient()
  var data = newMultipartData()
  data["q"] = genes.join(",")
  data["species"] = species
  data["scopes"] = "symbol,entrezgene,ensemblgene,retired"

  var r = C.postContent("http://mygene.info/v3/query", multipart=data)
  var js = parseJson(r)
  var ids = newSeq[string]()
  for res in js:
    if "_id" in res:
      ids.add($res["_id"])
    elif "notfound" in res:
      stderr.write_line &"""[seqcover] {res["query"]} not found, skipping"""

  data = newMultipartData()
  data["species"] = species
  data["ids"] = ids.join(",")
  if hg19:
    data["fields"] = "name,symbol,exons_hg19"
  else:
    data["fields"] = "name,symbol,exons"


  for res in C.postContent("http://mygene.info/v3/gene", multipart=data).parseJson:
    var gene = Gene(symbol: $res["symbol"], description: $res["name"])

    if hg19:
      gene.transcripts = to(res["exons_hg19"], seq[Transcript])
    else:
      gene.transcripts = to(res["exons"], seq[Transcript])

    gene.drop_alt_chroms
    result.add(gene)


proc get_glob_samples*(paths: seq[string]): seq[string] =
  result = newSeqOfCap[string](paths.len)
  for i, p in paths:
    var n = 0
    for w in p.walkFiles:
      n.inc
      result.add(w)
    if n == 0:
      raise newException(OSError, "[seqcover]: file:" & p & " not found")


proc read_dps*(paths: seq[string]): seq[Cover] =
  result = newSeq[Cover](paths.len)
  for i, p in paths:
    result[i] = open_dp(p)

proc read_d4s_to_table*(paths: seq[string]): TableRef[string, Cover] =
  result = newTable[string, Cover]()
  for pat in paths:
    for p in walkPattern(pat):
      var name = splitFile(p).name
      if name in result:
        raise newException(OSError, "[seqcover] repeated d4 name:" & p)
      if name.endsWith(".bed"):
        name = name[0..<name.len - 4]
      if name.endsWith(".per-base"):
        name = name[0..<name.len - 9]
      var d:Cover = p.open_dp
      result[name] = d

proc check_same_lengths*(dps: var seq[Cover], chrom: string, length: uint32, fai:Fai) =
  var lens = dps[0].chromosomes(fai)
  for d in dps.mitems:
    var clens = d.chromosomes(fai)

    doAssert lens[chrom] == clens[chrom], "[seqcover] differing chromosome lengths among samples"

when isMainModule:
  import json
  import algorithm
  import sequtils

  when true:
     #var ge = Gene(symbol: "MUC5B", description: "mucin 5B, oligomeric mucus/gel-forming", transcripts: @[Transcript(cdsstart: 1223123, cdsend: 1261608, chr: "11", position: @[[1223065, 1223193], [1225680, 1225737], [1226204, 1226276], [1226614, 1226876], [1227030, 1227145], [1227307, 1227398], [1227674, 1227781], [1228563, 1228765], [1229169, 1229295], [1229689, 1229807], [1230004, 1230143], [1230489, 1230600], [1230935, 1231005], [1231422, 1231560], [1231995, 1232160], [1232449, 1232544], [1232643, 1232770], [1233012, 1233268], [1233792, 1233848], [1234204, 1234305], [1234528, 1234680], [1235084, 1235223], [1235302, 1235413], [1236385, 1236562], [1236924, 1237164], [1238870, 1239027], [1239437, 1239566], [1239798, 1239943], [1240044, 1240088], [1240177, 1240375], [1240850, 1251743], [1252342, 1252524], [1252808, 1252980], [1254091, 1254351], [1254693, 1254880], [1255040, 1255266], [1255382, 1255558], [1256155, 1256225], [1256670, 1256771], [1257239, 1257271], [1257529, 1257710], [1258098, 1258203], [1258329, 1258367], [1258941, 1259061], [1259755, 1259842], [1259962, 1260085], [1260350, 1260393], [1260625, 1260728], [1261388, 1262172]], strand: 1, transcript: "NM_002458", txstart: 1223065, txend: 1262172)])
    var d4s = read_d4s_to_table(toSeq(walkPattern("d4s/*.d4")))
    var fa:Fai

    doAssert fa.open("/data/human/Homo_sapiens_assembly38.fasta")

    var backgrounds = read_d4s_to_table(toSeq(walkPattern("seqcover-backgrounds/*.d4")))
    var gpt: seq[GenePlotData]
    for gene in get_genes(@["PIGA", "KCNQ2", "MUC5B", "ARX", "DNM1", "SLC25A22", "CDKL5", "GABRA1", "ITPA", "GABRB1", "FRRS1L", "SLC12A5", "SIK1", "GRIN2D", "FGF12", "CAD", "DENND5A", "MDH2", "SCN1B", "AP3B2", "DENND5A", "MDH2", "GUF1", "YWHAG", "CNPY3", "CPLX1", "RHOBTB2"]):
    #for gene in get_genes(@["PIGA"]):
      echo gene

      #[
      let u = gene.transcripts.union
      stderr.write_line ">>> u:", u
      #for t in gene.transcripts:
      #  if t.transcript == "NM_001324238":
      #    stderr.write_line ">>> t:", $t
      stderr.write_line "######################"

      #var t = gene.transcripts[1]
      #var u = gene.transcripts.union

      #var last = t.position[^1]
      #t.position = @[last]
      #stderr.write_line "after:", $t

      #stderr.write_line u.translate(t, extend=10)
      #stderr.write_line "################"

      var pd = gene.plot_data(d4s, backgrounds, extend=10, fai=fa, max_gap=50)
      gpt.add(pd)

      for i, x in pd.plot_coords.x:
        if i == 0: continue
        doAssert x >= pd.plot_coords.x[i-1]
        doAssert pd.plot_coords.g[i] >= pd.plot_coords.g[i-1]


      var dlen:int
      for s, d in pd.plot_coords.depths.mpairs:
        dlen = d.len
        break
      for s, d in pd.plot_coords.background_depths.mpairs:
        doAssert d.len == dlen

      for i, p in pd.unioned_transcript.position:
        var idx = pd.plot_coords.x.lowerBound(p[0].uint32)
        if u.position[i][0].int != pd.plot_coords.g[idx].int:
          stderr.write_line &"i[{i}][0]:"
          stderr.write_line &"ERR LEFT : idx:{idx}, pd.plot_coords.g[idx]: {pd.plot_coords.g[idx]}, u.position[i][0]: {u.position[i][0]} l..r: {pd.plot_coords.g[idx-1..idx+1]} diff: {u.position[i][0].int - pd.plot_coords.g[idx].int}"

        idx = pd.plot_coords.x.lowerBound(p[1].uint32)
        if u.position[i][1].int != pd.plot_coords.g[idx].int:
          stderr.write_line &"i[{i}][1]:"
          stderr.write_line &"ERR RIGHT idx:{idx}, pd.plot_coords.g[idx]: {pd.plot_coords.g[idx]}, u.position[i][1]: {u.position[i][1]} l..r: {pd.plot_coords.g[idx-1..idx+1]} diff {u.position[i][1].int - pd.plot_coords.g[idx].int}"
          stderr.write_line "OK"


      #for t in pd.transcripts:
      #  if t.transcript == "NM_001324238":
      #    stderr.write_line ">>> t:", $t

      stderr.write_line " unioned:", pd.unioned_transcript

    gpt.sort(proc(a, b: GenePlotData): int = cmp(a.symbol, b.symbol))
    echo "plot_data = ", (%(gpt))
    ]#



