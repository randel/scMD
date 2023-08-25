# Author: Manqi Cai
###################################################
#' @title get_sig
#' @description  Generate Deconvolution Signatures
#'
#' This function generates deconvolution signatures based on the provided parameters.
#' It uses differentially methylated CpG sites from single-cell DNA methylation (scDNAm) data
#' to identify markers for various cell types.
#'
#' @param nmrk An integer specifying the number of marker CpG sites to consider for each cell type.
#' @param beta_mtx A matrix representing the scDNAm data, CpG sites by cell types. Default is NULL.
#' If NULL, the function will utilize the inbuilt scDNAm dataset.
#' @param DM_df A data frame detailing differentially methylated CpG sites from scDNAm data.
#' It should include columns for TargetID and celltype_ind, where columns corresponding to celltype_ind
#' should have p-values for each cell type.
#' @param ct_ind A character vector listing the cell types to include in the analysis.
#' Default is c("Astro","Micro", "Neuro","Oligo").
#'
#' @return A matrix of signatures that can be used for deconvolution, CpG sites by cell types.
#' @export
#'
get_sig <- function(nmrk, sig = sig_all, DM_df = Mrk_twosided,
                    ct_ind = c("Astro","Micro", "Neuro","Oligo")) {
  DM_df <- DM_df %>% filter(TargetID %in% rownames(sig))

  output_mrk <- lapply(ct_ind, function(i) {
    com <- DM_df %>%
      arrange(.data[[i]]) %>%
      dplyr::slice(1:nmrk)
    return(com$TargetID)
  })

  output_mrk <- unique(unlist(output_mrk))
  sig <- sig[output_mrk, ]

  return(sig)
}
