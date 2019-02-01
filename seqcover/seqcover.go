package main

import (
	"bufio"
	"encoding/json"
	"fmt"
	"io"
	"log"
	"os"
	"os/exec"
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

type Gene struct {
	Name            string `json:"gene"`
	Chrom           string `json:"chrom"`
	Exons           []Exon `json:"exons"`
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

func readDepthLines(path string, chrom string, exon *Exon) {
	args := []string{path}
	args = append(args, fmt.Sprintf("%s:%d-%d", chrom, exon.Start+1, exon.Stop))
	cmd := exec.Command("tabix", args...)
	cmd.Stderr = os.Stderr
	stdout, err := cmd.StdoutPipe()
	check(err)
	check(cmd.Start())

	b := bufio.NewReader(stdout)
	hasColon := true // are we parsing from quantized e.g. 7:10, or per-base. e.g. 7

	for {
		line, err := b.ReadString('\n')
		if len(line) != 0 {
			toks := strings.Split(strings.TrimSpace(line), "\t")
			if toks[3] == "0" || toks[3] == "0:1" {
				continue
			}
			c := toks[0]
			s := mustAtoi(toks[1])
			e := mustAtoi(toks[2])
			var v int

			if hasColon {
				ab := strings.SplitN(toks[3], ":", 2)
				if len(ab) == 1 {
					hasColon = false
				}
				v = mustAtoi(ab[0])
			} else {
				v = mustAtoi(toks[3])
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
	check(cmd.Wait())
	if len(exon.Depths) == 0 {
		exon.Depths = make([]uint16, exon.Length())
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

func readGenes(f string, perbase string, soi_path string, pop_coverage_tmpl string, geneNames []string) map[string]*Gene {
	fh, err := xopen.Ropen(f)
	check(err)
	var soi, last_pop *bix.Bix
	var last_pop_chrom string
	if soi_path != "" {
		soi, err = bix.New(soi_path)
	}
	genes := make(map[string]*Gene, 20)
	for {
		line, err := fh.ReadString('\n')
		if len(line) != 0 {
			toks := strings.SplitN(strings.TrimSpace(line), "\t", 6)
			if !contains(toks[3], geneNames) {
				continue
			}
			gene, ok := genes[toks[3]]
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
			readDepthLines(perbase, gene.Chrom, &e)
			if soi != nil {
				readVariants(soi, gene.Chrom, &e)
			}
			if last_pop == nil || gene.Chrom != last_pop_chrom {
				last_pop, err = bix.New(strings.Replace(pop_coverage_tmpl, "$chrom", gene.Chrom, 1))
				check(err)
				last_pop_chrom = gene.Chrom
			}
			readPopDepths(last_pop, gene.Chrom, &e)
			gene.Exons = append(gene.Exons, e)
			genes[toks[3]] = gene

		}
		if err == io.EOF {
			break
		}
		check(err)
	}

	for _, g := range genes {
		g.CalcStats()
	}

	fh.Close()
	return genes
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

func readVariants(b *bix.Bix, chrom string, exon *Exon) {

	iter, err := b.Query(ip{chrom: chrom, start: exon.Start, stop: exon.Stop})
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
	iter, err := b.Query(ip{chrom: chrom, start: exon.DataStart, stop: exon.DataStart + len(exon.Depths) - 1})
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

var geneNames = []string{
	"MECP2",
}

type Genes struct {
	// keyed by gene-name
	Genes map[string][]Gene `json:"exons,omitempty"`
}

func main() {
	cli := &Args{}
	arg.MustParse(cli)

	genes := readGenes(cli.GeneBed, cli.PerBase, cli.VariantsOfInterest, cli.PopulationCoverage, geneNames)

	// TODO: cull variants of interest to those that overlap any exon.
	//meta.Coverage = readCoverage(cli.PopulationCoverage, meta.Exons)

	//v, err := json.MarshalIndent(genes, "", "  ")
	v, err := json.Marshal(genes)
	check(err)
	fmt.Println(string(v))

}
