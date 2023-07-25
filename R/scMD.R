#' scMD
#'
#' Perform signature generation and ensemble deconvolution using the scMD algorithm.
#'
#' @param dat A beta matrix containing the scDNAm data.
#' @param bulk A matrix of bulk data.
#' @param nc An integer specifying the number of cores for parallel computation. Default is 5.
#' @param DM_df A data frame containing information about differentially methylated CpG sites from scDNAm data.
#' It should include columns for TargetID and celltype_ind. The columns celltype_ind contain p-values for each cell type.
#' @param dmet_list A character vector specifying the list of deconvolution methods. Default is c("CIBERSORT", "EPIC", "FARDEEP", "DCQ",
#'                 "DeconRNASeq", "BayesPrism", "ICeDT", "dtangle", "hspe").
#' @param include_minfi A logical value indicating whether to include Minfi-based deconvolution. Default is TRUE.
#' @param nmrk Number of markers per cell type used in deconvolution. Default is 100.
#' @param celltype_ind A character vector specifying the cell types to be included in the scMD analysis. Default is NULL, which
#'                     includes all cell types from the input data.
#' @param includeMval A logical value indicating whether to include M values (beta values transformed to M values) for deconvolution.
#'                    Default is TRUE.
#' @param gen_sig A logical value indicating whether to generate signatures for deconvolution. Default is TRUE. When set to FALSE, it means
#'                the input data is already preprocessed to only include markers for deconvolution.
#' @param output_path A character specifying the path to save the output results. Default is the current directory (".").
#'
#' @return A list containing various results including phat_all (estimated cell-type fractions), scMD_p (scMD estiamted CTS fractions),
#'         nmrk (number of markers used), and sig (signatures used in the deconvolution).
#'
#' @examples
#' # Assuming 'dat',  'bulk', 'DM_df', and 'sig_all' are defined
#' scMD(dat = my_dat, bulk = my_bulk, DM_df = my_DM_df, sig_all = my_sig_all)
scMD <- function(dat,bulk,nc =5,
                 DM_df = NULL,dmet_list = c("CIBERSORT","EPIC","FARDEEP","DCQ","DeconRNASeq","BayesPrism","ICeDT","dtangle","hspe"),
                 include_minfi = TRUE,nmrk = 100,
                 celltype_ind =  NULL,
                 includeMval = T, gen_sig = T,
                 output_path = "."){

  overlap_features = intersect(rownames(dat),rownames(bulk))
  dat <- dat[overlap_features,]
  bulk <- bulk[overlap_features,]


  if(is.null(celltype_ind)){
    celltype_ind = colnames(dat)
  }

  if(gen_sig){
    sig <- get_sig(nmrk = nmrk,dat,DM_df,ct_ind = celltype_ind)
  }else{
    sig <- dat
  }

  bulk <- bulk[intersect(rownames(sig),rownames(bulk)),]

  if(include_minfi){
    frac_est <-  est_frac(sig, bulk)


    frac_rpc <-  EpiDISH::epidish(bulk, sig,method = 'RPC')$estF

    minfi_res <- DNAm_minfi_deconv_general(dat = sig,metaref = data.frame(SamplesName = colnames(sig), deconv_clust = colnames(sig)),
                                           bulk =bulk,true_frac = NULL,findMrks = F)
    if(includeMval){
      phat_all <- list(NNLS = frac_est, Houseman_beta = minfi_res[[2]]$counts, Houseman_m = minfi_res[[1]]$counts,RPC = frac_rpc)
    }else{
      phat_all <- list(NNLS = frac_est, Houseman_beta = minfi_res[[2]]$counts,RPC = frac_rpc)
    }
    saveRDS(phat_all, file = paste0(output_path,"phat_all.rds"))
  }


    ref_list = list(beta = list(ref_matrix = sig, meta_ref = data.frame(SamplesName = colnames(sig), deconv_clust = colnames(sig)),
                                data_name = "beta"))

    suppressWarnings(res_Ens <- EnsDeconv(count_bulk = bulk,
                                          ref_list  = ref_list,
                                          params = get_params(data_type = "singlecell-rna", data_name = "beta",
                                                              n_markers = 50, Marker.Method = "none",
                                                              TNormalization = c("none"), CNormalization = c("none"),
                                                              Scale = c("linear","log"), dmeths = dmet_list),
                                          outpath = paste0(output_path,"beta"),parallel_comp = T,ncore = nc))



    phat_list1 <- unlist(lapply(res_Ens$allgene_res, getphat),recursive = FALSE)

    if(includeMval){
      ref_list = list(Mval = list(ref_matrix = BetaToMvalue(sig), meta_ref = data.frame(SamplesName = colnames(sig), deconv_clust = colnames(sig)),
                                  data_name = "Mval"))

      suppressWarnings(res_Ens_mval <- EnsDeconv(count_bulk = BetaToMvalue(bulk),
                                                 ref_list  = ref_list,
                                                 params = get_params(data_type = "singlecell-rna", data_name = "Mval",
                                                                     n_markers = 50, Marker.Method = "none",
                                                                     TNormalization = c("none"), CNormalization = c("none"),
                                                                     Scale = c("linear"), dmeths = dmet_list),
                                                 outpath = output_path,parallel_comp = T,ncore = nc))

      phat_list2 <- unlist(lapply(res_Ens_mval$allgene_res, getphat),recursive = FALSE)
      phat_list1 <- append(phat_list1, phat_list2)
    }


    if(include_minfi){
      phat_all <- append(phat_all, phat_list1)
    }else{
      phat_all <- phat_list1
    }


  saveRDS(phat_all, file = paste0(output_path,"phat_all.rds"))
  scMD_p = CTS_EnsDeconv_wrapper2(phat_all = phat_all)[["ensemble_p"]]

  return(list(phat_all = phat_all, scMD_p=scMD_p,nmrk = nmrk,sig = sig))
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
