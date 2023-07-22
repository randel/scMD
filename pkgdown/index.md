single cell Methylation Deconvolution (scMD)
=============================================
a cellular deconvolution framework to reliably estimate cell type fractions from tissue-level DNAm data with single-cell DNAm (scDNAm) references. 

Introduction
-------------------
`scMD` is designed to generate reliable references from scDNAm. Evidenced by rigorous benchmarking across multiple datasets, it excels in estimating cellular fractions from bulk DNAm data and in identifying Alzheimer's disease-associated cell types. This package provides processed signatures for instant use and also feasible code to conduct DNAm deconvolution.

Installation
-------------------

```r
devtools::install_github("randel/scMD")
```

Workflow
-----------------

<img src = "./man/figures/scMDFigure1.jpg">


Tutorial
-----------------

https://randel.github.io/scMD/




Reference
-----------------

Cai M, Zhou J, McKennan C, Wang J. (2023). scMD: cell type deconvolution using single-cell DNA methylation references.
