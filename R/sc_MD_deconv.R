sc_MD_deconv <- function(){
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

  if("Houseman" %in% dmet_list){
    frac_est <-  est_frac(sig, bulk)


    frac_rpc <-  EpiDISH::epidish(bulk, sig,method = 'RPC')$estF

    minfi_res <- DNAm_minfi_deconv_general(dat = sig,metaref = data.frame(SamplesName = colnames(sig), deconv_clust = colnames(sig)),
                                           bulk =bulk,true_frac = NULL,findMrks = F)

    phat_all <- list(NNLS = frac_est, Houseman_beta = minfi_res[[2]]$counts, Houseman_m = minfi_res[[1]]$counts,RPC = frac_rpc)

    saveRDS(phat_all, file = paste0(output_path,"phat_all.rds"))
  }


  ref_list = list(beta = list(ref_matrix = sig, meta_ref = data.frame(SamplesName = colnames(sig), deconv_clust = colnames(sig)),
                              data_name = "beta"))

  suppressWarnings(res_Ens <- gen_all_res_list(count_bulk = bulk,
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
                                             outpath = output_path,parallel_comp = T,ncore = nc))

  phat_list2 <- unlist(lapply(res_Ens_mval$allgene_res, getphat),recursive = FALSE)
  phat_list1 <- append(phat_list1, phat_list2)
  phat_all <- append(phat_all, phat_list1)






  saveRDS(phat_all, file = paste0(output_path,"phat_all.rds"))
  scMD_p = CTS_EnsDeconv_wrapper2(phat_all = phat_all)[["ensemble_p"]]
}
