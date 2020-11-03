[![Build Status](https://github.com/brentp/seqcover/workflows/tests/badge.svg?branch=master)](https://github.com/brentp/seqcover/actions)

seqcover is a tool for viewing and evaluating depth-of-coverage with the following aims. It should:

 - show a global view where it's easy to see problematic samples and genes
 - offer an interactive gene-wise view to explore coverage characteristics of individual samples within each gene
 - **not** require a server (single html page)
 - be responsive for up to 20 samples * 200 genes and be useful for a single-sample [see how we do this](#how-it-works)
 - highlight outlier samples based on any number of (summarized) background samples

It is available as a static linux binary.

## Example Output

https://brentp.github.io/seqcover/

### Usage

`seqcover` can accept per-base coverage files in [d4](https://github.com/38/d4-format) or bgzipped bedgaph format. Either of
these formats can be output by [mosdepth](https://github.com/brentp/mosdepth) but `d4` format will be much faster.

Generate a report:
```
# generate per base depth files if you don't have them. can also parallize this...
mkdir -p samples/
for b in *.bam; do
  n=$(basename $b .bam)
  mosdepth -x -t 4 samples/$n $b
done
seqcover report --genes PIGA,KCNQ2,ARX,DNM1,SLC25A22,CDKL5,GABRA1,CAD,MDH2,SCN1B,CNPY3,CPLX1,NEB,HNRNPA1,CCDC39,AIFM1,CHCHD10 \
		 --background seqcover/seqcover_p5.d4 \
		 --fasta $fasta samples/*.bed.gz \
		 -r my_genes_report.html
```

Generate a background level:
```
seqcover generate-background --percentile 5 -f $fasta -o seqcover/ d4s/HG00*.d4
```
Once generated, this can be sent to `seqcover report` to give a metric for each sample
of the number of bases below the 5th percentile of the backgrounds, which is a nice quality-control
value. 

These backgounds should be specific to the samples of interest, so if you are using exome data
the backgrounds should be generated from samples from the same exome capture kit (and prefereably
from the same sequencing center).

Generate a transcript file:
```
# generate a transcripts file that can be used as input to the report option with the --transcripts-file option
# this is useful if running on servers with restricted internet access

seqcover save-transcripts --genes PIGA,KCNQ2,ARX,DNM1,SLC25A22,CDKL5,GABRA1,CAD,MDH2,SCN1B,CNPY3,CPLX1,NEB,HNRNPA1,CCDC39,AIFM1,CHCHD10 \
		 --output-path transcripts.json \
		 --hg19

# use the file as input to report with the --transcripts-file option
seqcover report --genes PIGA,KCNQ2,ARX,DNM1,SLC25A22,CDKL5,GABRA1,CAD,MDH2,SCN1B,CNPY3,CPLX1,NEB,HNRNPA1,CCDC39,AIFM1,CHCHD10 \
		 --background seqcover/seqcover_p5.d4 \
		 --fasta $fasta samples/*.bed.gz \
		 -r my_genes_report.html \
		 --transcripts-file transcripts.json
```

## How It Works

### Performance

`seqcover` is a command-line tool that extracts depth information for requested genes and generates a terse report.
This is possible because we **excise introns** which are often the majority of bases in the gene. The user can specify to extend
into the intron beyond the default of 10 bases. As that extension increases, `seqcover` will show more and more of the intronic space.
But it's possible to display 20 samples * 100 genes in html because we remove introns and use webGL (via plotly) for rendering.

### Outliers

An outlier is data-dependent. Every exome will appear as an outlier given a set of genomes and every 30X genome will
appear as an outlier given a set of 60X genomes. Therefore, `seqcover` let's the user extract a depth percentile from a set
of chosen backgrounds. For example, the default is to extract the 5th percentile from a set of samples. This percentile can then
be shown in the report and used as a metric: **how many bases in each sample are below the 5th percentile**.

