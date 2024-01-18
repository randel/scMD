# Author: Manqi Cai
###################################################
#' @title scMD framework
#'
#' @description This function performs signature generation and ensemble deconvolution using the scMD algorithm.
#' @param bulk A matrix of bulk data, with CpGs by sample.
#' @param bulk_type A character string indicating the technical platform for the bulk data. Available options include "450k_or_850k" or "WGBS",.
#' @param ncluster An integer specifying the number of cores for parallel computation. Default is 5.
#' @param dmet_list A character vector specifying the desired deconvolution methods to be used.
#' If NULL (the default), the following methods will be included:
#' c("CIBERSORT", "EPIC", "FARDEEP", "DCQ", "ICeDT", "Houseman", "NNLS", "RPC").
#' @param nmrk The number of markers per cell type used in deconvolution. Default is 100.
#' @param enableFileSaving Enable Saving of Intermediate Output
#'  (Optional) A boolean flag that controls the saving of intermediate outputs as separate files. 
#'          When set to TRUE, intermediate outputs of the analysis will be saved to files. 
#'          If not explicitly set, this parameter defaults to FALSE, meaning that intermediate 
#'          outputs will not be saved by default.
#' @param output_path A character specifying the path for saving intermediate results. Default is NULL.
#'
#' @return If scMD converges, return:
#' \itemize{
#'    \item {scMD_p: matrix, scMD estiamted CTS fractions.}
#'    \item {phat_all: list, estimated cell-type fractions.}
#'    }
#'
#' @examples
#' data("Guintivano")
#' Est.prop.scMD = scMD(bulk = Guintivano_bulk_sub, bulk_type = "450k_or_850k",output_path = "~/Desktop/")
#'
#' @importFrom EpiDISH epidish
#' @import dplyr
#' @import tidyverse
#' @import IlluminaHumanMethylation450kanno.ilmn12.hg19
#' @import EnsDeconv
#' @importFrom minfi makeGenomicRatioSetFromMatrix
#' @importFrom MIND est_frac
#' @importFrom quadprog solve.QP
#'
#' @export
scMD <- function(bulk, bulk_type = "450k_or_850k", ncluster =5,
                 dmet_list = NULL, nmrk = 100, enableFileSaving = FALSE,output_path = NULL){

  if(is.null(dmet_list)){
    dmet_list <- c("CIBERSORT", "EPIC", "FARDEEP", "DCQ", "ICeDT", "Houseman", "NNLS", "RPC")
  }

  if(bulk_type == "450k_or_850k"){
    data("Lee_7ct_450850")
    data("Tian_7ct_450850")

    Lee_res <- sc_MD_deconv(bulk,sc_mtx = Lee_sig_all,DM_df = Lee_DF_450850,ncluster,
                             dmet_list,nmrk,NULL,enableFileSaving,output_path)

    Tian_res <- sc_MD_deconv(bulk,sc_mtx = Tian_sig_all,DM_df = Tian_DF_450850,ncluster,
                             dmet_list,nmrk,NULL,enableFileSaving,output_path)
  }else{
    data("Lee_7ct_WGBS")
    data("Tian_7ct_WGBS")

    Lee_res <- sc_MD_deconv(bulk,sc_mtx = Lee_sig_all_WGBS,DM_df = Lee_DF_WGBS,ncluster,
                            dmet_list,nmrk,NULL,enableFileSaving,output_path)

    Tian_res <- sc_MD_deconv(bulk,sc_mtx = Tian_sig_all_WGBS,DM_df = Tian_DF_WGBS,ncluster,
                             dmet_list,nmrk,NULL,enableFileSaving,output_path)
  }

  # Call CTS_EnsDeconv_wrapper function
  scMD_p <- CTS_EnsDeconv_wrapper(phat_all = c(Lee_res, Tian_res))[["ensemble_p"]]


  return(list(scMD_p=scMD_p,phat_all = c(Lee_res, Tian_res)))
}

