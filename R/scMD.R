# Author: Manqi Cai
###################################################
#' @title scMD framework
#'
#' @description This function performs signature generation and ensemble deconvolution using the scMD algorithm.
#' @param bulk A matrix of bulk data, with CpGs by sample.
#' @param bulk_type A character string indicating the technical platform for the bulk data. Available options include "450k_or_850k", "WGBS", or "NULL".
#' When set to "NULL", users are required to provide both `sc_mtx` and `DM_df`.
#' @param use_sc A character specifying the inbuilt scDNAm to be used. Default is "Both". Options are "Luo", "Tian", "Both", or "NULL".
#' When set to NULL, users must provide `sc_mtx` and `DM_df`.
#' @param ncluster An integer specifying the number of cores for parallel computation. Default is 5.
#' @param dmet_list A character vector specifying the desired deconvolution methods to be used.
#' If NULL (the default), the following methods will be included:
#' c("CIBERSORT", "EPIC", "FARDEEP", "DCQ", "ICeDT", "Houseman", "NNLS", "RPC").
#' @param nmrk The number of markers per cell type used in deconvolution. Default is 100.
#' @param gen_sig A logical value indicating whether to generate signatures for deconvolution. Default is TRUE. When FALSE,
#' input data should already be preprocessed to only include markers for deconvolution.
#' @param celltype_ind A character vector specifying the cell types included in the scMD analysis. Default is NULL, which includes all cell types from the input data.
#' @param output_path A character specifying the path for saving intermediate results. Default is the current directory (".").
#' @param sc_mtx A beta matrix containing the scDNAm data. Default is NULL. When NULL, scMD uses the inbuilt scDNAm data.
#' @param DM_df A data frame with information about differentially methylated CpG sites from scDNAm data, including columns for TargetID and celltype_ind.
#' Columns for celltype_ind should contain p-values for each cell type. Default is NULL.
#' When NULL, scMD uses the inbuilt scDNAm data, or the provided `sc_mtx` is considered preprocessed.
#'
#' @return If scMD converges, return:
#' \itemize{
#'    \item {scMD_p: matrix, scMD estiamted CTS fractions.}
#'    \item {phat_all: list, estimated cell-type fractions.}
#'    \item {sig: matriix, signatures used in the deconvolution.}
#'    }
#'
#' @examples
#' # Assuming 'my_bulk' and 'my_bulk_type' are defined
#' scMD(bulk = my_bulk, bulk_type = my_bulk_type)
#'
#' @importFrom EpiDISH epidish
#'
#' @export
scMD <- function(bulk, bulk_type = "450k_or_850k",use_sc = "Both",ncluster =5,dmet_list = NULL,
                 nmrk = 100,gen_sig = T,celltype_ind =  NULL,output_path = ".",sc_mtx = NULL,DM_df = NULL,
                  ){
  if(bulk_type == "450k_or_850k"){
    load
  }
  overlap_features = intersect(rownames(sc_mtx),rownames(bulk))
  sc_mtx <- sc_mtx[overlap_features,]
  bulk <- bulk[overlap_features,]


  if(is.null(celltype_ind)){
    celltype_ind = colnames(sc_mtx)
  }

  if(gen_sig){
    sig <- get_sig(nmrk = nmrk,sc_mtx,DM_df,ct_ind = celltype_ind)
  }else{
    sig <- sc_mtx
  }

  bulk <- bulk[intersect(rownames(sig),rownames(bulk)),]


  frac_est <-  est_frac(sig, bulk)


  frac_rpc <-  EpiDISH::epidish(bulk, sig,method = 'RPC')$estF

  minfi_res <- DNAm_minfi_deconv_general(dat = sig,metaref = data.frame(SamplesName = colnames(sig), deconv_clust = colnames(sig)),
                                         bulk =bulk,true_frac = NULL,findMrks = F)
  phat_all <- list(NNLS = frac_est, Houseman_beta = minfi_res[[2]]$counts, Houseman_m = minfi_res[[1]]$counts,RPC = frac_rpc)

  saveRDS(phat_all, file = paste0(output_path,"phat_all.rds"))



    ref_list = list(beta = list(ref_matrix = sig, meta_ref = data.frame(SamplesName = colnames(sig), deconv_clust = colnames(sig)),
                                data_name = "beta"))

    suppressWarnings(res_Ens <- EnsDeconv(count_bulk = bulk,
                                          ref_list  = ref_list,
                                          params = get_params(data_type = "singlecell-rna", data_name = "beta",
                                                              n_markers = 50, Marker.Method = "none",
                                                              TNormalization = c("none"), CNormalization = c("none"),
                                                              Scale = c("linear","log"), dmeths = dmet_list),
                                          outpath = paste0(output_path,"beta"),parallel_comp = T,ncore = ncluster))



    phat_list1 <- unlist(lapply(res_Ens$allgene_res, getphat),recursive = FALSE)

    ref_list = list(Mval = list(ref_matrix = BetaToMvalue(sig), meta_ref = data.frame(SamplesName = colnames(sig), deconv_clust = colnames(sig)),
                                data_name = "Mval"))

    suppressWarnings(res_Ens_mval <- EnsDeconv(count_bulk = BetaToMvalue(bulk),
                                               ref_list  = ref_list,
                                               params = get_params(data_type = "singlecell-rna", data_name = "Mval",
                                                                   n_markers = 50, Marker.Method = "none",
                                                                   TNormalization = c("none"), CNormalization = c("none"),
                                                                   Scale = c("linear"), dmeths = dmet_list),
                                               outpath = output_path,parallel_comp = T,ncore = ncluster))

    phat_list2 <- unlist(lapply(res_Ens_mval$allgene_res, getphat),recursive = FALSE)
    phat_list1 <- append(phat_list1, phat_list2)

    phat_all <- append(phat_all, phat_list1)



  saveRDS(phat_all, file = paste0(output_path,"phat_all.rds"))
  scMD_p = CTS_EnsDeconv_wrapper2(phat_all = phat_all)[["ensemble_p"]]

  return(list(scMD_p=scMD_p,phat_all = phat_all, sig = sig))
}

get_sig <- function(nmrk,sig = sig_all,DM_info = Mrk_twosided,
                    ct_ind = c("Astro","Micro_Endo", "Neuro","Oligo")){
  DM_info <- DM_info %>% filter(TargetID%in% rownames(sig))

  output_mrk <- lapply(ct_ind, function(i){
    com <- DM_info%>% arrange(.data[[i]]) %>%
      dplyr::slice(1:nmrk)
    return(com$TargetID)
  })
  output_mrk <- unique(unlist(output_mrk))
  sig <- sig[output_mrk,]
  return(sig)
}
