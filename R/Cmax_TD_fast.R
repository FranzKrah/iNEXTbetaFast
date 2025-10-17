#' Fast computation of Cmax_joint for taxonomic diversity (TD) for abundance data
#' @param com A community data matrix (species x sites)
#' @param ncores Number of CPU cores to use for parallel processing
#' @param chunk_size Number of pairwise combinations to process per chunk (default = 5000)
#' @return The minimum Cmax_joint value across all pairwise combinations
#' @import iNEXT.3D iNEXT.beta3D
#' @importFrom parallel mclapply
#' @details This function computes the Cmax_joint for taxonomic diversity which is the minimal sample coverage between all pairs of the sites.
#' It splits the pairwise combinations of species into chunks and processes them in parallel across multiple CPU cores.
#' This approach significantly reduces computation time for large datasets. Note that this function only runs on UNIX based systems (Linux, MacOS) due to the use of `mclapply`. Runtime: On Ubuntu 24.04 with 192GB RAM on 16 cores, processing 200 sites with >20000 OTUs takes approximately 90 seconds. 1700 sites in 2280 seconds.
#' @references Chao, A., Thorn, S., Chiu, C.-H., Moyes, F., Hu, K.-H., Chazdon, R. L., Wu, J., Magnago, L. F. S., Dornelas, M., Zeleny, D., Colwell, R. K., and Magurran, A. E. (2023). Rarefaction and extrapolation with beta diversity under a framework of Hill numbers: the iNEXT.beta3D standardization. Ecological Monographs e1588.
#' @export

Cmax_TD_fast <- function(com, ncores = 16, chunk_size = 5000) {

  com <- as.data.frame(com)
  pairs <- combn(seq_len(ncol(com)), 2, simplify = FALSE)

  # Split all pairwise combinations into roughly equal chunks
  chunks <- split(pairs, ceiling(seq_along(pairs) / chunk_size))

  # Forked parallelization — each core processes one or several chunks
  res_min <- mclapply(
    chunks,
    function(chunk) {
      local_min <- Inf
      for (idx in chunk) {
        val <- tryCatch(
          iNEXT.3D:::TDinfo(as.vector(cbind(com[[idx[1]]], com[[idx[2]]])),
                            datatype = "abundance")$`SC(2n)`,
          error = function(e) NA_real_
        )
        if (!is.na(val) && val < local_min) local_min <- val
      }
      local_min
    },
    mc.cores = ncores,
    mc.preschedule = TRUE
  )

  # Global minimum across all chunks
  min(unlist(res_min), na.rm = TRUE)
}
