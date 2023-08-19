#' Lee et al. brain single-cell DNAm dataset
#' Processed scDNAm Data from Lee et al. (450k and 850k platforms)
#'
#' This dataset comprises the beta matrix of Lee's scDNAm data (`Lee_sig_all`) and
#' a data frame with details about differentially methylated CpG sites (`Lee_DF_450850`).
#'
#' @references Lee, D. -S., Luo, C., Zhou, J., Chandran, S., Rivkin, A., Bartlett, A., Nery,
#' J. R., Fitzpatrick, C., O’Connor, C., Dixon, J. R., et al. (2019).
#' Simultaneous Profiling of 3D Genome Structure and DNA Methylation in Single Human Cells. Nature Methods, 16, 999–1006.
#' @format Lee_sig_all: A matrix of dimensions 883,132 x 7. Columns represent cell types:
#' \describe{
#'   \item{"Astro"}{Astrocytes}
#'   \item{"Micro"}{Microglia}
#'   \item{"Endo"}{Endothelial cells}
#'   \item{"Oligo"}{Oligodendrocytes}
#'   \item{"OPC"}{Oligodendrocyte precursor cells}
#'   \item{"Inh"}{Inhibitory neurons}
#'   \item{"Exc"}{Excitatory neurons}
#' }
#' @format `Lee_DF_450850`:
#' A matrix with dimensions of 883,131 x 8. Columns provide details on:
#' \describe{
#'   \item{"TargetID"}{ID for the CpG site}
#'   \item{"Astro", "Micro", "Endo", "Oligo", "OPC", "Inh", "Exc"}{p-values indicating differential methylation for each cell type}
#' }
#' "Lee_7ct_450850"

#' Processed scDNAm Data from Lee et al. (WGBS platforms)
#'
#' This dataset comprises the beta matrix of Lee's scDNAm data (`Lee_sig_all_WGBS`) and
#' a data frame with details about differentially methylated CpG sites (`Lee_DF_WGBS`).
#' For computational efficiency and optimized storage, only the top 100,000 CpG
#' sites per cell type were retained, based on cell type-specific p-values.
#'
#' @references Lee, D. -S., Luo, C., Zhou, J., Chandran, S., Rivkin, A., Bartlett, A., Nery,
#' J. R., Fitzpatrick, C., O’Connor, C., Dixon, J. R., et al. (2019).
#' Simultaneous Profiling of 3D Genome Structure and DNA Methylation in Single Human Cells. Nature Methods, 16, 999–1006.
#'
#' @format `Lee_sig_all_WGBS`:
#' A matrix of dimensions 1,940,080 x 7. Columns represent cell types:
#' \describe{
#'   \item{"Astro"}{Astrocytes}
#'   \item{"Micro"}{Microglia}
#'   \item{"Endo"}{Endothelial cells}
#'   \item{"Oligo"}{Oligodendrocytes}
#'   \item{"OPC"}{Oligodendrocyte precursor cells}
#'   \item{"Inh"}{Inhibitory neurons}
#'   \item{"Exc"}{Excitatory neurons}
#' }
#'
#' @format `Lee_DF_WGBS`:
#' A matrix of dimensions 1,940,281 x 8. Columns detail:
#' \describe{
#'   \item{"TargetID"}{ID for the CpG site}
#'   \item{"Astro", "Micro", "Endo",  "Oligo", "OPC",   "Inh",   "Exc"}{p-values indicating differential methylation for each cell type}
#' }
#' "Lee_7ct_WGBS"

#' Processed scDNAm Data from Tian et al. (450k and 850k platforms)
#'
#' @references Tian, W., Zhou, J., Bartlett, A., Zeng, Q., Liu, H., Castanon, R. G., Kenworthy,
#' M., Altshul, J., Valadon, C., Aldridge, A., et al. (2022).
#' Epigenomic Complexity of the Human Brain Revealed by Single-Cell DNA Methylomes and 3D Genome Structures. bioRxiv, 2022–11.
#'
#' This dataset comprises the beta matrix of Tian's scDNAm data (`Tian_sig_all`) and
#' a data frame with details about differentially methylated CpG sites (`Tian_DF_450850`).
#' @format `Tian_sig_all`:
#' A matrix of dimensions 819,694 x 7. Columns represent cell types:
#' \describe{
#'   \item{"Astro"}{Astrocytes}
#'   \item{"Micro"}{Microglia}
#'   \item{"Endo"}{Endothelial cells}
#'   \item{"Oligo"}{Oligodendrocytes}
#'   \item{"OPC"}{Oligodendrocyte precursor cells}
#'   \item{"Inh"}{Inhibitory neurons}
#'   \item{"Exc"}{Excitatory neurons}
#' }
#'
#' @format `Tian_DF_450850`:
#' A matrix of dimensions 819,694 x 8. Columns detail:
#' \describe{
#'   \item{"TargetID"}{ID for the CpG site}
#'   \item{"Astro", "Micro", "Endo",  "Oligo", "OPC",   "Inh",   "Exc"}{p-values indicating differential methylation for each cell type}
#' }
#' "Tian_7ct_450850"

#' Processed scDNAm Data from Tian et al. (WGBS platforms)
#'
#' This dataset comprises the beta matrix of Tian's scDNAm data (`Tian_sig_all_WGBS`) and
#' a data frame with details about differentially methylated CpG sites (`Tian_DF_WGBS`).
#' For computational efficiency and optimized storage, only the top 100,000 CpG
#' sites per cell type were retained, based on cell type-specific p-values.
#'
#' @references Tian, W., Zhou, J., Bartlett, A., Zeng, Q., Liu, H., Castanon, R. G., Kenworthy,
#' M., Altshul, J., Valadon, C., Aldridge, A., et al. (2022).
#' Epigenomic Complexity of the Human Brain Revealed by Single-Cell DNA Methylomes and 3D Genome Structures. bioRxiv, 2022–11.
#'
#' @format `Tian_sig_all_WGBS`:
#' A matrix of dimensions 554,668 x 7. Columns represent cell types:
#' \describe{
#'   \item{"Astro"}{Astrocytes}
#'   \item{"Micro"}{Microglia}
#'   \item{"Endo"}{Endothelial cells}
#'   \item{"Oligo"}{Oligodendrocytes}
#'   \item{"OPC"}{Oligodendrocyte precursor cells}
#'   \item{"Inh"}{Inhibitory neurons}
#'   \item{"Exc"}{Excitatory neurons}
#' }
#'
#' @format `Tian_DF_WGBS`:
#' A matrix of dimensions 554,668 x 8. Columns detail:
#' \describe{
#'   \item{"TargetID"}{ID for the CpG site}
#'   \item{"Astro", "Micro", "Endo",  "Oligo", "OPC",   "Inh",   "Exc"}{p-values indicating differential methylation for each cell type}
#' }
#'
#' "Tian_7ct_WGBS"
