#' Title
#'
#' @param nmrk
#' @param sig
#' @param DM_info
#' @param ct_ind
#'
#' @return
#' @export
#'
#' @examples
get_sig <- function(nmrk,sig = sig_all,DM_info = Mrk_twosided,
                    ct_ind = c("Astro","Micro", "Neuro","Oligo")){
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
