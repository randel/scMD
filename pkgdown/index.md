single cell Methylation Deconvolution (scMD)
=============================================
a cellular deconvolution framework to reliably estimate cell type fractions from tissue-level DNAm data with single-cell DNAm (scDNAm) references. 

Introduction
-------------------
`scMD` is designed to generate reliable references from scDNAm. Evidenced by rigorous benchmarking across multiple datasets, it excels in estimating cellular fractions from bulk DNAm data and in identifying Alzheimer's disease-associated cell types. This package provides processed signatures for instant use and also feasible code to conduct DNAm deconvolution.

Installation
-------------------
To run the complete tutorial, please install the scMD package.
```r
# install devtools if necessary
install.packages('devtools')

# install the scMD package
devtools::install_github('randel/scMD')

# load
library(scMD)
```

Data
-------------------
The scMD package offers a range of processed signatures obtained from the studies conducted by Lee et al. and Tian et al. These signatures are available in two versions: one comprising signatures under the 850k platform and the other under the 450k platform.

To access the processed signatures, please navigate to the following link: https://github.com/randel/scMD/tree/main/Processed_data_450k850k



Workflow
-------------------

<img src = "./man/figures/scMDFigure1.png">


Tutorial
-----------------

https://randel.github.io/scMD/




Reference
-----------------

1. Cai, Manqi, et al. "scMD: cell type deconvolution using single-cell DNA methylation references." bioRxiv (2023): 2023-08.

2. Lee, D. -S., Luo, C., Zhou, J., Chandran, S., Rivkin, A., Bartlett, A., Nery, J. R., Fitzpatrick, C., O’Connor, C., Dixon, J. R., et al. (2019). Simultaneous Profiling of 3D Genome Structure and DNA Methylation in Single Human Cells. Nature Methods, 16, 999–1006.

3. Tian, W., Zhou, J., Bartlett, A., Zeng, Q., Liu, H., Castanon, R. G., Kenworthy, M., Altshul, J., Valadon, C., Aldridge, A., et al. (2022). Epigenomic Complexity of the Human Brain Revealed by Single-Cell DNA Methylomes and 3D Genome Structures. bioRxiv, 2022–11.
