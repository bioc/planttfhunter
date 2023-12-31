

#' Identify TFs and classify them in families
#' 
#' @param domain_annotation A 2-column data frame with the gene ID in the first
#' column and the domain ID in the second column.
#' 
#' @return A 2-column data frame with the variables \strong{Gene} 
#' and \strong{Family} representing gene ID and TF family, respectively.
#' @export
#' @rdname classify_tfs
#' @examples 
#' data(gsu_annotation)
#' domain_annotation <- gsu_annotation
#' families <- classify_tfs(domain_annotation)
classify_tfs <- function(domain_annotation = NULL) {
    
    dom <- c(
        list_domains("dbd"), list_domains("auxiliary"), 
        list_domains("forbidden")
    )
    annot <- domain_annotation[domain_annotation$Domain %in% dom, ]
    annot_list <- split(annot, annot[, 1])
    
    fams <- Reduce(rbind, lapply(annot_list, function(x) {
        domains <- x$Domain
        
        # Apply all helper classification functions to object domains
        funcs <- list(
            check_ap2_erf, check_b3, check_c2c2, check_garp,
            check_hb, check_mads, check_myb, check_nf_y, check_smallfams
        )
        tffam <- unlist(lapply(funcs, function(f) f(domains)))
        if(length(tffam) == 0) {
            fam_df <- NULL
        } else {
            fam_df <- data.frame(
                Gene = unique(x[, 1]),
                Family = tffam
            )
            
        }
        return(fam_df)
    }))
    return(fams)
}



