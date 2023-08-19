sc_MD_deconv <- function(bulk, sc_mtx, DM_df,ncluster =5,dmet_list,
                         nmrk = 100,celltype_ind =  NULL,output_path = ".",gen_sig = T){
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

  phat_all <- list()

  if("Houseman" %in% dmet_list) {
    minfi_res <- DNAm_minfi_deconv_general(
      dat = sig,metaref = data.frame(SamplesName = colnames(sig), deconv_clust = colnames(sig)),
      bulk = bulk,true_frac = NULL,findMrks = F)
    phat_all <- append(phat_all, list(Houseman_beta = minfi_res[[2]]$counts, Houseman_m = minfi_res[[1]]$counts))

  } else if("NNLS" %in% dmet_list) {
    frac_est <- est_frac(sig, bulk)
    phat_all <- append(phat_all, list(NNLS = frac_est))

  } else if("RPC" %in% dmet_list) {
    frac_rpc <- epidish(bulk, sig, method = 'RPC')$estF
    phat_all <- append(phat_all, list(RPC = frac_rpc))
  }

  if(length(phat_all) > 1) {
    saveRDS(phat_all, file = paste0(output_path, "phat_all.rds"))
  }

  dmet_list <- setdiff(dmet_list, c("Houseman", "NNLS", "RPC"))

  # Create reference list for "beta"
  ref_list_beta <- list(
    beta = list(
      ref_matrix = sig,
      meta_ref = data.frame(SamplesName = colnames(sig), deconv_clust = colnames(sig)),
      data_name = "beta"
    )
  )

  # Generate results for "beta"
  suppressWarnings(
    res_Ens <- gen_all_res_list( count_bulk = bulk, ref_list = ref_list_beta,
                                 params = get_params(
                                   data_type = "singlecell-rna", data_name = "beta", n_markers = 50,  Marker.Method = "none",
                                   TNormalization = c("none"),  CNormalization = c("none"), Scale = c("linear","log"),
                                   dmeths = dmet_list ), outpath = paste0(output_path,"beta"),
                                 parallel_comp = T, ncore = ncluster ) )

  # Extract phat values for "beta"
  phat_list1 <- unlist(lapply(res_Ens$allgene_res, getphat), recursive = FALSE)

  # Create reference list for "Mval"
  ref_list_Mval <- list(
    Mval = list(
      ref_matrix = BetaToMvalue(sig),
      meta_ref = data.frame(SamplesName = colnames(sig), deconv_clust = colnames(sig)),
      data_name = "Mval"
    )
  )

  # Generate results for "Mval"
  suppressWarnings(
    res_Ens_mval <- gen_all_res_list(
      count_bulk = BetaToMvalue(bulk), ref_list = ref_list_Mval,
      params = get_params( data_type = "singlecell-rna",  data_name = "Mval",
                           n_markers = 50,  Marker.Method = "none", TNormalization = c("none"),
                           CNormalization = c("none"), Scale = c("linear"),  dmeths = dmet_list ),
      outpath = output_path, parallel_comp = T, ncore = nc ) )

  # Extract and combine phat values for "beta" and "Mval"
  phat_list2 <- unlist(lapply(res_Ens_mval$allgene_res, getphat), recursive = FALSE)
  phat_list1 <- append(phat_list1, phat_list2)

  phat_all <- append(phat_all, phat_list1)

  # Save the combined phat values
  saveRDS(phat_all, file = paste0(output_path,"phat_all.rds"))

  return(phat_all)

}


