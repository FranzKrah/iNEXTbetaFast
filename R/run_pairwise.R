#' Run pairwise computations in parallel using forked processes
#' @param com A community matrix (samples x species)
#' @param fun A function to apply to each pair of species
#' @param level A numeric value indicating the level parameter for the function
#' @param ncores Number of cores to use for parallel processing
#' @param PDtree An optional phylogenetic tree object
#' @param chunk_size Number of pairs to process in each chunk
#' @return A data.table with the results of the pairwise computations
#' @references Chao, A., Thorn, S., Chiu, C.-H., Moyes, F., Hu, K.-H., Chazdon, R. L., Wu, J., Magnago, L. F. S., Dornelas, M., Zeleny, D., Colwell, R. K., and Magurran, A. E. (2023). Rarefaction and extrapolation with beta diversity under a framework of Hill numbers: the iNEXT.beta3D standardization. Ecological Monographs e1588.
#' @details Only use with ncores > 1.
#' @importFrom utils combn
#' @importFrom parallel mclapply
#' @importFrom data.table rbindlist
#' @export

run_pairwise <- function(com, fun, level = 0.8, ncores = 12, PDtree = NULL, chunk_size = 5000) {


  com_mat <- as.matrix(com)
  pairs <- combn(colnames(com_mat), 2, simplify = FALSE)

  # Split into chunks
  chunks <- split(pairs, ceiling(seq_along(pairs) / chunk_size))

  # Run computation in parallel using forked processes
  res_list <- mclapply(chunks, function(chunk) {
    lapply(chunk, function(pair_cols) {
      args <- list(pair_cols = pair_cols, com = com_mat, level = level)
      if (!is.null(PDtree)) {
        args$PDtree <- PDtree
        args$reft <- reft
      }
      do.call(fun, args)
    })
  }, mc.cores = ncores, mc.cleanup = TRUE, mc.preschedule = TRUE)

  # Flatten results incrementally into a data.table
  res_dt <- rbindlist(lapply(res_list, function(chunk_res) {
    rbindlist(lapply(chunk_res, function(x) as.list(x)))
  }))

  return(res_dt)
}
