import os
import d4
import httpclient
import strformat
import json
import strutils


type Transcript* = object
  cdsstart*: int
  cdsend*: int
  `chr`: string
  position*: seq[array[2, int]]
  strand*: int
  transcript*: string
  txstart*:int
  txend*:int

#type Transcripts* = seq[Transcript]

type Gene* = object
  symbol*: string
  description*: string
  transcripts*: seq[Transcript]


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

proc get_genes(genes:seq[string], species:string="human"): seq[Gene] =
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
    echo gene
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

proc check_same_lengths*(d4s: seq[D4], chrom: string, length: uint32) =
  for d in d4s:
    doAssert d.chromosomes[chrom] == length, "[seqcover] differing chromosome lengths among samples"

when isMainModule:

  for gene in get_genes(@["KCNQ2", "MUC5B", "XXXXX"]):
    echo "gene"
    echo gene

