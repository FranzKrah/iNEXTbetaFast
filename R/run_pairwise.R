#' Run pairwise similarity computations in parallel using forked processes
#' @param com A community matrix (samples x species)
#' @param fun A function to apply to each pair of species
#' @param level A numeric value indicating the level parameter for the function
#' @param ncores Number of cores to use for parallel processing
#' @param PDtree An optional phylogenetic tree object
#' @param chunk_size Number of pairs to process in each chunk, if NULL (default) then the size is computed automatically
#' @return A data.table with the results of the pairwise computations of similarities according to \code{fun}
#' @references Chao, A., Thorn, S., Chiu, C.-H., Moyes, F., Hu, K.-H., Chazdon, R. L., Wu, J., Magnago, L. F. S., Dornelas, M., Zeleny, D., Colwell, R. K., and Magurran, A. E. (2023). Rarefaction and extrapolation with beta diversity under a framework of Hill numbers: the iNEXT.beta3D standardization. Ecological Monographs e1588.
#' @details
#' Only use with ncores > 1.
#' The code was originally developed by Oliver Mitesser.
#'
#' @importFrom utils combn
#' @importFrom parallel mclapply
#' @importFrom data.table rbindlist
#' @examples
#' # library(microeco)
#' # data(otu_table_ITS)
#' # beta_td <- run_pairwise(otu_table_ITS, pair_c_a_td, ncores = 10)
#' @export

run_pairwise <- function(com, fun, level = 0.8, ncores = 12, PDtree = NULL, chunk_size = NULL) {

  com_mat <- as.matrix(com)
  pairs <- combn(colnames(com_mat), 2, simplify = FALSE)
  total_pairs <- length(pairs)

  # --- auto-tune chunk size ---
  if (is.null(chunk_size)) {
    # Aim for about 3 to 5 times as many chunks as cores for good load balance
    target_chunks <- ncores * 4
    chunk_size <- ceiling(total_pairs / target_chunks)
    # Keep it within reasonable bounds
    chunk_size <- max(200, min(5000, chunk_size))
    message(sprintf("Auto-selected chunk size = %d with approx %d total pairs and %d cores",
                    chunk_size, total_pairs, ncores))
  }

  # Split into chunks
  chunks <- split(pairs, ceiling(seq_along(pairs) / chunk_size))

  # Run computation in parallel using forked processes
  res_list <- mclapply(chunks, function(chunk) {
    lapply(chunk, function(pair_cols) {
      args <- list(pair_cols = pair_cols, com = com_mat, level = level)
      if (!is.null(PDtree)) {
        args$PDtree <- PDtree
      }
      do.call(fun, args)
    })
  }, mc.cores = ncores, mc.cleanup = TRUE, mc.preschedule = TRUE)

  # Flatten results incrementally into a data.table
  res_dt <- rbindlist(lapply(res_list, function(chunk_res) {
    rbindlist(lapply(chunk_res, function(x) as.list(x)))
  }))

  res_dt[, (3:5) := lapply(.SD, function(x) {
    x[x > 1] <- 1
    x[x < 0] <- 0
    x
  }), .SDcols = 3:5]

  return(res_dt)
}
