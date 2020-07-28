import os
import d4
import httpclient
import strformat
import json
import strutils
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

proc get_genes*(genes:seq[string], species:string="human"): seq[Gene] =
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
  data["fields"] = "name,symbol,exons"

  for res in C.postContent("http://mygene.info/v3/gene", multipart=data).parseJson:
    var gene = Gene(symbol: $res["symbol"], description: $res["name"])

    gene.transcripts = to(res["exons"], seq[Transcript])
    result.add(gene)



  #var gene_id = $js["hits"][0]["_id"]
  #gene_id = gene_id.strip(chars={'"'})

  #echo &"http://mygene.info/v3/gene/{gene_id}?fields=name,symbol,exons"
  #var exons = C.getContent(&"http://mygene.info/v3/gene/{gene_id}?fields=name,symbol,exons")
  #echo exons



proc get_glob_samples*(paths: seq[string]): seq[string] =
  result = newSeqOfCap[string](paths.len)
  for i, p in paths:
    var n = 0
    for w in p.walkFiles:
      n.inc
      result.add(w)
    if n == 0:
      raise newException(OSError, "[seqcover]: file:" & p & " not found")


proc read_d4s*(paths: seq[string]): seq[D4] =
  result = newSeq[D4](paths.len)
  for i, p in paths:
    if not result[i].open(p):
      raise newException(OSError, "[seqcover] couldn't open d4 file:" & p)

proc read_d4s_to_table*(paths: seq[string]): TableRef[string, D4] =
  result = newTable[string, D4]()
  for p in paths:
    let name = splitFile(p).name
    if name in result:
      raise newException(OSError, "[seqcover] repeated d4 name:" & p)
    var d:D4
    if not d.open(p):
      raise newException(OSError, "[seqcover] couldn't open d4 file:" & p)
    result[name] = d

proc check_same_lengths*(d4s: seq[D4], chrom: string, length: uint32) =
  for d in d4s:
    doAssert d.chromosomes[chrom] == length, "[seqcover] differing chromosome lengths among samples"

when isMainModule:
  import json

  when true:
     #var ge = Gene(symbol: "MUC5B", description: "mucin 5B, oligomeric mucus/gel-forming", transcripts: @[Transcript(cdsstart: 1223123, cdsend: 1261608, chr: "11", position: @[[1223065, 1223193], [1225680, 1225737], [1226204, 1226276], [1226614, 1226876], [1227030, 1227145], [1227307, 1227398], [1227674, 1227781], [1228563, 1228765], [1229169, 1229295], [1229689, 1229807], [1230004, 1230143], [1230489, 1230600], [1230935, 1231005], [1231422, 1231560], [1231995, 1232160], [1232449, 1232544], [1232643, 1232770], [1233012, 1233268], [1233792, 1233848], [1234204, 1234305], [1234528, 1234680], [1235084, 1235223], [1235302, 1235413], [1236385, 1236562], [1236924, 1237164], [1238870, 1239027], [1239437, 1239566], [1239798, 1239943], [1240044, 1240088], [1240177, 1240375], [1240850, 1251743], [1252342, 1252524], [1252808, 1252980], [1254091, 1254351], [1254693, 1254880], [1255040, 1255266], [1255382, 1255558], [1256155, 1256225], [1256670, 1256771], [1257239, 1257271], [1257529, 1257710], [1258098, 1258203], [1258329, 1258367], [1258941, 1259061], [1259755, 1259842], [1259962, 1260085], [1260350, 1260393], [1260625, 1260728], [1261388, 1262172]], strand: 1, transcript: "NM_002458", txstart: 1223065, txend: 1262172)])
    var d4s = read_d4s_to_table(@["d4s/HG00096.final.d4", "d4s/HG00097.final.d4", "d4s/HG00099.final.d4"])
    var gpt: seq[GenePlotData]
    for gene in get_genes(@["PIGA", "KCNQ2", "MUC5B", "ARX", "DNM1", "SLC25A22", "CDLK5", "GABRA1", "ITPA", "GABRB1"]):

      stderr.write_line ">>> u:", gene.transcripts.union
      for t in gene.transcripts:
        if t.transcript notin ["NM_172108", "NM_172109"]: continue
        stderr.write_line ">>> t:", $t
      stderr.write_line "######################"

      #var t = gene.transcripts[1]
      #var u = gene.transcripts.union

      #var last = t.position[^1]
      #t.position = @[last]
      #stderr.write_line "after:", $t

      #stderr.write_line u.translate(t, extend=10)
      #stderr.write_line "################"

      var pd = gene.plot_data(d4s, extend=10)
      gpt.add(pd)

    echo "plot_data = ", (%(gpt))



