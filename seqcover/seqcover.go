package main

import (
	"bytes"
	"encoding/json"
	"io"
	"log"
	"strconv"
	"strings"

	arg "github.com/alexflint/go-arg"
	"github.com/brentp/bix"
	"github.com/brentp/irelate/interfaces"
	"github.com/brentp/irelate/parsers"
	"github.com/brentp/vcfgo"
	"github.com/brentp/xopen"
)

var populationFields = []string{"mean", "median", "5", "10", "15", "20", "30", "100"}
var excludedPopulationFields = []string{"5", "10", "15", "30", "100"}

type Args struct {
	VariantsOfInterest string `arg:"-v,help:vcf of variants of interest"`
	PopulationCoverage string `arg:"-p,help:path to population coverage files in gnomAD/ExAC format as a template such as: coverage/gnomad.exomes.r2.0.2.chr$chrom.coverage.txt.gz "`
	OutputPrefix       string `arg:"-o,help:prefix for output json files"`
	GeneBed            string `arg:"positional,required,help:path to bed file of exons with gene names"`
	PerBase            string `arg:"positional,required,help:path to per-base.bed.gz from mosdepth"`
	// GeneList...
}

type SiteOfInterest struct {
	Start int    `json:"start"`
	Stop  int    `json:"stop"`
	Text  string `json:"text"`
}

type Exon struct {
	Start int `json:"start"`
	Stop  int `json:"stop"`
	// we can get data before and after the exon so DataStart stores how far before
	// and how far after is implicit in len(Depths)
	DataStart int      `json:"datastart"`
	Depths    []uint16 `json:"depths"`
	// keyed by metric ,e.g. mean/median/percent at 1x,5x, etc.
	PopulationDepths map[string][]float32 `json:"populationdepths,omitempty"`
	SitesOfInterest  []SiteOfInterest     `json:"sitesofinterest,omitempty"`
}

type Stats struct {
	PercentBasesUnder30X                        float32
	PercentBasesUnder7X                         float32
	PercentBasesUnder7XWithGnomad99PercentAt20X float32
	PercentClinvarVariantsUnder7X               float32
}

// Gene contains the per-gene information including the full per-base exon info.
type Gene struct {
	Name            string `json:"gene"`
	Chrom           string `json:"chrom"`
	Exons           []Exon `json:"exons"`
	SampleStats     *Stats `json:"stats"`
	PopulationStats *Stats `json:"popstats"`
}

// SmallGene does not contain the per-base info.
// TODO: don't need popstats for every sample.
type SmallGene struct {
	Name            string `json:"gene"`
	Chrom           string `json:"chrom"`
	SampleStats     *Stats `json:"stats"`
	PopulationStats *Stats `json:"popstats"`
}

func (s *Stats) CheckDepth(exon *Exon, d uint16, i int) {
	// only use actual exonic bases
	if d < 30 {
		s.PercentBasesUnder30X++
		if d < 7 {
			s.PercentBasesUnder7X++
			if exon.PopulationDepths["20"][i] > 0.9 {
				s.PercentBasesUnder7XWithGnomad99PercentAt20X++
			}
		}
	}
}

func (g *Gene) CalcStats() {
	g.SampleStats = &Stats{}
	g.PopulationStats = &Stats{}
	exonicBases := 0
	exonicSitesOfInterest := 0
	for _, exon := range g.Exons {
		exonicBases += exon.Length()
		//log.Println(exon.PopulationDepths["20"])
		if len(exon.Depths) != len(exon.PopulationDepths["20"]) {
			log.Println(len(exon.Depths), len(exon.PopulationDepths["20"]))
			panic("expected equal-length depths")
		}
		off := exon.Start - exon.DataStart
		for i, d := range exon.Depths {
			if i < off {
				continue
			}
			if i > exon.Stop-exon.DataStart {
				break
			}
			g.SampleStats.CheckDepth(&exon, d, i)
			popd := exon.PopulationDepths["mean"][i]
			g.PopulationStats.CheckDepth(&exon, uint16(popd), i)

		}
		for _, soi := range exon.SitesOfInterest {
			var starti = soi.Start - exon.DataStart
			var stopi = soi.Stop - exon.DataStart
			if stopi > len(exon.Depths) {
				stopi = len(exon.Depths)
			}
			if starti < 0 {
				starti = 0
			}
			exonicSitesOfInterest++
			under := false
			for i := starti; i < stopi; i++ {
				if exon.Depths[i] < 7 {
					g.SampleStats.PercentClinvarVariantsUnder7X++
					under = true
				}
				if exon.PopulationDepths["mean"][i] < 7 {
					g.PopulationStats.PercentClinvarVariantsUnder7X++
				}
				if under {
					break
				}

			}

		}

	}
	g.SampleStats.Adjust(float32(exonicBases), exonicSitesOfInterest)
	g.PopulationStats.Adjust(float32(exonicBases), exonicSitesOfInterest)
}

func (s *Stats) Adjust(exonicBases float32, exonicSitesOfInterest int) {
	s.PercentBasesUnder30X /= float32(exonicBases)
	s.PercentBasesUnder7X /= float32(exonicBases)
	s.PercentBasesUnder7XWithGnomad99PercentAt20X /= float32(exonicBases)
	s.PercentClinvarVariantsUnder7X /= float32(max(1, exonicSitesOfInterest))

	s.PercentBasesUnder30X *= 100
	s.PercentBasesUnder7X *= 100
	s.PercentBasesUnder7XWithGnomad99PercentAt20X *= 100
	s.PercentClinvarVariantsUnder7X *= 100
}

func check(e error) {
	if e != nil {
		panic(e)
	}
}

func mustAtoi(a string) int {
	if v, err := strconv.Atoi(a); err == nil {
		return v
	} else {
		panic(err)
	}
}

func contains(needle string, haystack []string) bool {
	for _, h := range haystack {
		if h == needle {
			return true
		}
	}
	return false
}

func writeBigGenes(output_prefix string, big map[string]*Gene, small map[string]*SmallGene) {
	for g, v := range big {
		v.CalcStats()
		small[g] = &SmallGene{Name: v.Name, Chrom: v.Chrom, SampleStats: v.SampleStats, PopulationStats: v.PopulationStats}

		f, err := xopen.Wopen(output_prefix + v.Name + ".json")
		check(err)
		j, err := json.Marshal(v)
		check(err)
		f.Write(j)
		f.Close()
	}

}

func (e Exon) Length() int {
	return e.Stop - e.Start
}

func (e Exon) Overlaps(chrom string, start, end int) bool {
	if e.Start > end || e.Stop < start {
		return false
	}
	return true
}

// meets IPosition interface
type ip struct {
	chrom string
	start int
	stop  int
}

func (i ip) Start() uint32 {
	return uint32(i.start)
}
func (i ip) End() uint32 {
	return uint32(i.stop)
}
func (i ip) Chrom() string {
	return i.chrom
}

func readDepthLines(b *bix.Bix, chrom string, exon *Exon) {
	iter, err := b.FastQuery(ip{chrom: chrom, start: exon.Start, stop: exon.Stop})
	check(err)
	hasColon := true // are we parsing from quantized e.g. 7:10, or per-base. e.g. 7
	for {
		r, err := iter.Next()
		if r != nil {
			var iv = r.(*parsers.Interval)
			if bytes.Compare(iv.Fields[3], []byte{'0'}) == 0 || bytes.Compare(iv.Fields[3], []byte{'0', ':', '1'}) == 0 {
				continue
			}
			//toks := []string{string(iv.Fields[0]), string(iv.Fields[1]), string(iv.Fields[2]), string(iv.Fields[3])}
			c := string(iv.Fields[0])
			s := mustAtoi(string(iv.Fields[1]))
			e := mustAtoi(string(iv.Fields[2]))
			var v int
			var t3 = string(iv.Fields[3])

			if hasColon {
				ab := strings.SplitN(t3, ":", 2)
				if len(ab) == 1 {
					hasColon = false
				}
				v = mustAtoi(ab[0])
			} else {
				v = mustAtoi(t3)
			}

			if exon.Overlaps(c, s, e) {
				exon.Update(c, s, e, v)
			}
		}
		if err == io.EOF {
			break
		}
		check(err)
	}
}

func readVariants(b *bix.Bix, chrom string, exon *Exon) {

	iter, err := b.FastQuery(ip{chrom: chrom, start: exon.Start, stop: exon.Stop})
	check(err)

	for {
		r, err := iter.Next()
		if err == io.EOF {
			break
		}
		check(err)
		iv := r.(interfaces.VarWrap).IVariant
		v := iv.(*vcfgo.Variant)

		dbn, err := v.Info_.Get("CLNDN")
		if err != nil {
			dbn, err = v.Info_.Get("CLNDBN")
			if err != nil {
				continue
			}
		}

		var dbns string
		var ok bool
		if dbns, ok = dbn.(string); !ok {

			if dbnsa, ok2 := dbn.([]string); ok2 {
				dbns = strings.Join(dbnsa, ",")
			} else {
				continue
			}

		}

		text := "rsid:" + v.Id() + "<br>clinvar:" + dbns + "<br>position:" + strconv.Itoa(int(r.Start()))
		exon.SitesOfInterest = append(exon.SitesOfInterest, SiteOfInterest{Start: int(r.Start()), Stop: int(r.End()), Text: text})

	}

}

func readGenes(f string, perbase string, soi_path string, pop_coverage_tmpl string, output_prefix string) map[string]*SmallGene {
	fh, err := xopen.Ropen(f)
	check(err)

	bixperbase, err := bix.New(perbase)
	check(err)

	var soi, last_pop *bix.Bix
	var last_pop_chrom string
	if soi_path != "" {
		soi, err = bix.New(soi_path)
	}
	smallgenes := make(map[string]*SmallGene, 20)
	// biggenes is per gene, per-base info.
	biggenes := make(map[string]*Gene, 20)
	var last_chrom = ""
	var k = 0
	for {
		line, err := fh.ReadString('\n')
		if len(line) != 0 {
			k += 1
			if k%5000 == 0 {
				log.Println(k, line)
			}
			toks := strings.SplitN(strings.TrimSpace(line), "\t", 6)
			if toks[0] != last_chrom && last_chrom != "" {
				log.Println(last_chrom)
				writeBigGenes(output_prefix, biggenes, smallgenes)
				biggenes = make(map[string]*Gene, 20)
			}
			last_chrom = toks[0]
			gene, ok := biggenes[toks[3]]
			if !ok {
				gene = &Gene{Chrom: toks[0], Name: toks[3]}
			} else {
				if gene.Chrom != toks[0] {
					log.Printf("warning genes from multiple chromosomes found for %s", toks[3])
					continue
				}
			}
			e := Exon{Start: mustAtoi(toks[1]), Stop: mustAtoi(toks[2])}
			e.DataStart = e.Start
			readDepthLines(bixperbase, gene.Chrom, &e)
			if soi != nil {
				readVariants(soi, gene.Chrom, &e)
			}
			if (last_pop == nil || gene.Chrom != last_pop_chrom) && pop_coverage_tmpl != "" {
				last_pop, err = bix.New(strings.Replace(pop_coverage_tmpl, "$chrom", gene.Chrom, 1))
				check(err)
				last_pop_chrom = gene.Chrom
			}
			readPopDepths(last_pop, gene.Chrom, &e)
			gene.Exons = append(gene.Exons, e)
			biggenes[toks[3]] = gene

		}
		if err == io.EOF {
			break
		}
		check(err)
	}

	fh.Close()
	return smallgenes
}
func (exon *Exon) initDepths() {
	exon.PopulationDepths = make(map[string][]float32)
	for _, k := range populationFields {
		if contains(k, excludedPopulationFields) {
			continue
		}
		exon.PopulationDepths[k] = make([]float32, len(exon.Depths))
	}

}

func readPopDepths(b *bix.Bix, chrom string, exon *Exon) {
	if b == nil {
		return
	}
	iter, err := b.FastQuery(ip{chrom: chrom, start: exon.DataStart, stop: exon.DataStart + len(exon.Depths) - 1})
	check(err)
	exon.initDepths()

	for {
		r, err := iter.Next()
		if err == io.EOF {
			break
		}
		check(err)
		// TODO: figure out how to add text
		if pi, ok := r.(*parsers.Interval); ok {
			pos, err := strconv.Atoi(string(pi.Fields[1]))
			if err != nil {
				panic(err)
			}
			//pos--
			if pos == exon.Stop {
				break
			}

			for i, pf := range populationFields {
				if contains(pf, excludedPopulationFields) {
					continue
				}
				fv, err := strconv.ParseFloat(string(pi.Fields[i+2]), 32)
				if err != nil {
					panic(err)
				}
				exon.PopulationDepths[pf][pos-exon.DataStart] = float32(fv)
			}

		}

	}

}

func min(a, b int) int {
	if a < b {
		return a
	}
	return b
}

func max(a, b int) int {
	if a > b {
		return a
	}
	return b
}

func (l *Exon) Update(chrom string, start, end, value int) {
	if l.Start-start > 25 {
		start = l.Start - 25
	}
	if end-l.Stop > 25 {
		end = l.Stop + 25
	}

	l.DataStart = min(l.DataStart, start)
	start -= l.DataStart
	end -= l.DataStart

	if end > len(l.Depths) {
		//log.Printf("got end: %d for len %d with start: %d.", end, len(l.Depths), start)
		//log.Printf("len: %d. calc: %d", len(l.Depths), l.Stop-l.Start+1)
		var ext []uint16
		ext = make([]uint16, max(l.Length(), end)-len(l.Depths))
		l.Depths = append(l.Depths, ext...)
	}
	for i := start; i < end; i++ {
		if d := uint32(l.Depths[i]) + uint32(value); d > 65535 {
			l.Depths[i] = 65535
		} else {
			l.Depths[i] = uint16(d)
		}
	}
}

type Genes struct {
	// keyed by gene-name
	Genes map[string][]Gene `json:"exons,omitempty"`
}

func main() {
	cli := &Args{OutputPrefix: "seqcover-"}
	arg.MustParse(cli)

	genes := readGenes(cli.GeneBed, cli.PerBase, cli.VariantsOfInterest, cli.PopulationCoverage, cli.OutputPrefix)
	v, err := json.Marshal(genes)
	if err != nil {
		panic(err)
	}
	if fh, err := xopen.Wopen(cli.OutputPrefix + "summary.json"); err != nil {
		panic(err)
	} else {
		fh.Write(v)
	}

}
