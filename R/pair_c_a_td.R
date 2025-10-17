#' Compute a single pairwise dissimilarity estimates based on TD diversity
#' @param pair_cols A vector of two column names representing the pair of sites.
#' @param com A community data matrix (data.frame or matrix) with sites as columns and species as rows.
#' @param level A numeric value indicating the coverage level for diversity estimation.
#' @return A list containing site names and dissimilarity estimates (Sørensen, Horn, Morisita-Horn).
#' @details This function computes pairwise dissimilarity estimates between two sites using TD diversity measures.
#' It calculates gamma and alpha diversities and derives dissimilarity indices based on these values.
#' @references Chao, A., Thorn, S., Chiu, C.-H., Moyes, F., Hu, K.-H., Chazdon, R. L., Wu, J., Magnago, L. F. S., Dornelas, M., Zeleny, D., Colwell, R. K., and Magurran, A. E. (2023). Rarefaction and extrapolation with beta diversity under a framework of Hill numbers: the iNEXT.beta3D standardization. Ecological Monographs e1588.
#' @import iNEXT.3D iNEXT.beta3D
#' @export

pair_c_a_td <- function(pair_cols, com, level) {

  # Slice the two rows only
  data_pair <- com[, pair_cols , drop = FALSE]

  # Compute gamma and alpha via estimate3D
  data_gamma <- rowSums(data_pair)
  data_gamma <- data_gamma[data_gamma > 0]
  data_alpha <- as.vector(as.matrix(data_pair))

  gamma <- estimate3D(data_gamma, diversity = 'TD',
                      datatype = 'abundance', base = "coverage", level = level, nboot = 0)$qTD
  alpha <- estimate3D(data_alpha, diversity = 'TD',
                      datatype = 'abundance', base = "coverage", level = level, nboot = 0)$qTD

  # Compute beta
  beta_vec <- gamma / (alpha / 2)

  # Compute dissimilarities
  C_q0 <- (beta_vec[1]^(1 - 0) - 1) / (2^(1 - 0) - 1)
  C_q1 <- log(beta_vec[2]) / log(2)
  C_q2 <- (beta_vec[3]^(1 - 2) - 1) / (2^(1 - 2) - 1)

  # Return results as data.frame row
  list(
    site1 = pair_cols[1],
    site2 = pair_cols[2],
    sor_est = as.numeric(1 - C_q0),
    hor_est = as.numeric(1 - C_q1),
    mor_hor_est = as.numeric(1 - C_q2)
  )
}
