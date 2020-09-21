[![Build Status](https://github.com/brentp/seqcover/workflows/tests/badge.svg?branch=master)](https://github.com/brentp/seqcover/actions)

seqcover is a tool for viewing depth-of-coverage with the following aims. It should:

 - show a global view where it's easy to see problematic samples and genes
 - offer an interactive gene-wise view to explore coverage characteristics of individual samples within each gene
 - **not** require a server (single html page)
 - be responsive for up to 20 samples * 200 genes and be useful for a single-sample
 - highlight outlier samples based on any number of (summarized) background samples
