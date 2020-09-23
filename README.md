[![Build Status](https://github.com/brentp/seqcover/workflows/tests/badge.svg?branch=master)](https://github.com/brentp/seqcover/actions)

seqcover is a tool for viewing and evaluating depth-of-coverage with the following aims. It should:

 - show a global view where it's easy to see problematic samples and genes
 - offer an interactive gene-wise view to explore coverage characteristics of individual samples within each gene
 - **not** require a server (single html page)
 - be responsive for up to 20 samples * 200 genes and be useful for a single-sample
 - highlight outlier samples based on any number of (summarized) background samples

It is available as a static linux binary.

### Usage

`seqcover` can accept per-base coverage files in [d4](https://github.com/38/d4-format) or bgzipped bedgaph format. Either of
these formats can be output by [mosdepth](https://github.com/brentp/mosdepth) but `d4` format will be much faster.

Generate a report:
```
seqcover report PIGA,KCNQ2,ARX,DNM1,SLC25A22,CDKL5,GABRA1,CAD,MDH2,SCN1B,CNPY3,CPLX1,RHOBTB --background seqcover/seqcover_p5.d4 \
		 --fasta $fasta samples/*.bed.gz
```

Generate a background level:
```
seqcover generate-background -f $fasta -o seqcover/ d4s/HG00*.d4
```




