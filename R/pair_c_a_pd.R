#' Compute a single phylogenetic dissimilarity estimates between two sites
#' @param pair_cols A vector of two column names representing the pair of sites.
#' @param com A community data matrix (data.frame or matrix) with sites as columns and species as rows.
#' @param level A numeric value indicating the coverage level for diversity estimation.
#' @param PDtree A phylogenetic tree object of class 'phylo'.
#' @return A list containing site names and dissimilarity estimates (Sørensen, Horn, Morisita-Horn).
#' @details PDtype A character string indicating the type of phylogenetic diversity to compute (currently hard coded as "meanPD"). Note that this computation is considerably slower than the taxonomic diversity.
#' Runtime: Note that this computation is considerably slower than the taxonomic diversity \code{pairs_c_a_td}. The same dataset may easily take by a factor 50 longer.
#' @references Chao, A., Thorn, S., Chiu, C.-H., Moyes, F., Hu, K.-H., Chazdon, R. L., Wu, J., Magnago, L. F. S., Dornelas, M., Zeleny, D., Colwell, R. K., and Magurran, A. E. (2023). Rarefaction and extrapolation with beta diversity under a framework of Hill numbers: the iNEXT.beta3D standardization. Ecological Monographs e1588.
#' @import tibble
#' @importFrom ape node.depth.edgelength drop.tip
#' @export

pair_c_a_pd <- function(pair_cols, com, level, PDtree, reft) {

  # --- load helpers once ---
  PhD.m.est <- safe_get("PhD.m.est", "iNEXT.3D")
  phyBranchAL_Abu <- safe_get("phyBranchAL_Abu", "iNEXT.3D")
  coverage_to_size <- safe_get("coverage_to_size", "iNEXT.beta3D")

  # --- subset columns for this pair ---
  data_pair <- as.matrix(com[, pair_cols, drop = FALSE])

  # --- pool (gamma) data ---
  data_gamma <- rowSums(data_pair)
  data_gamma <- data_gamma[data_gamma > 0]

  # --- basic counts ---
  n <- sum(data_pair)

  # --- coverage-based sizes ---
  m_gamma <- coverage_to_size(data_gamma, level, datatype = "abundance")
  m_alpha <- coverage_to_size(as.vector(data_pair), level, datatype = "abundance")

  # --- gamma ---
  aL_gamma <- phyBranchAL_Abu(PDtree, data_gamma, rootExtend = TRUE, refT = reft)
  aL_gamma$treeNabu$branch.length <- aL_gamma$BLbyT[, 1]
  aL_table_gamma <- aL_gamma$treeNabu[, c("branch.abun", "branch.length", "tgroup")]

  gamma <- PhD.m.est(
    ai = aL_table_gamma$branch.abun,
    Lis = as.matrix(aL_table_gamma$branch.length),
    m = m_gamma, nt = n, q = c(0, 1, 2),
    reft = reft, cal = "PD"
  )

  # --- alpha (loop over two samples only) ---
  aL_list <- vector("list", 2)
  for (i in 1:2) {
    x <- data_pair[, i]
    x <- x[x > 0]
    aL <- phyBranchAL_Abu(PDtree, x, rootExtend = TRUE, refT = reft)
    aL$treeNabu$branch.length <- aL$BLbyT[, 1]
    aL_list[[i]] <- aL$treeNabu[, c("branch.abun", "branch.length", "tgroup")]
  }

  aL_table_alpha <- rbind(aL_list[[1]], aL_list[[2]])

  alpha <- PhD.m.est(
    ai = aL_table_alpha$branch.abun,
    Lis = as.matrix(aL_table_alpha$branch.length),
    m = m_alpha, q = c(0, 1, 2),
    nt = n, reft = reft, cal = "PD"
  ) / 2

  # --- normalize (meanPD-like) ---
  gamma <- gamma / reft
  alpha <- alpha / reft

  # --- compute beta and dissimilarities ---
  beta <- gamma / alpha

  C_q0 <- (beta[1]^(1 - 0) - 1) / (2^(1 - 0) - 1)
  C_q1 <- log(beta[2]) / log(2)
  C_q2 <- (beta[3]^(1 - 2) - 1) / (2^(1 - 2) - 1)

  list(
    site1 = pair_cols[1],
    site2 = pair_cols[2],
    sor_est = as.numeric(1 - C_q0),
    hor_est = as.numeric(1 - C_q1),
    mor_hor_est = as.numeric(1 - C_q2)
  )
}

# pair_c_a_pd <- function(pair_cols, com, level, PDtree) {
#
#
#   PhD.m.est <- safe_get("PhD.m.est", "iNEXT.3D")
#   phyBranchAL_Abu <- safe_get("phyBranchAL_Abu", "iNEXT.3D")
#   coverage_to_size <- safe_get("coverage_to_size", "iNEXT.beta3D")
#
#   # Slice the two rows only
#   data_pair <- com[, pair_cols , drop = FALSE]
#
#   # for Hmax total
#   pool.data = rowSums(data_pair)
#
#   pool.name = names(pool.data[pool.data>0])
#   tip = PDtree$tip.label[-match(pool.name, PDtree$tip.label)]
#   mytree = drop.tip(PDtree, tip)
#
#   # H_max = get.rooted.tree.height(mytree)
#   reft = max(node.depth.edgelength(mytree))
#
#   n = sum(data_pair)
#   # Compute gamma and alpha via estimate3D
#   data_gamma <- rowSums(data_pair)
#   data_gamma <- data_gamma[data_gamma > 0]
#   data_alpha <- as.vector(as.matrix(data_pair))
#
#   m_gamma = coverage_to_size(data_gamma, level, datatype='abundance')
#   m_alpha = coverage_to_size(data_alpha, level, datatype='abundance')
#
#   pool.data = data_pair %>% data.frame %>% rownames_to_column()
#   # pool.data = full_join(pool.data, data_pair %>% data.frame %>% rownames_to_column(), 'rowname')
#   pool.data[is.na(pool.data)] = 0
#   pool.data = pool.data %>% column_to_rownames() %>% rowSums
#
#   pool.name = names(pool.data[pool.data>0])
#   tip = PDtree$tip.label[-match(pool.name, PDtree$tip.label)]
#   mytree = drop.tip(PDtree, tip)
#
#
#   aL <- phyBranchAL_Abu(phylo = PDtree, data = data_gamma, rootExtend = TRUE, refT = reft)
#   aL$treeNabu$branch.length = aL$BLbyT[,1]
#   suppressMessages({
#     aL_table_gamma = aL$treeNabu %>% select(branch.abun, branch.length, tgroup)
#   })
#
#   gamma = PhD.m.est(ai = aL_table_gamma$branch.abun, Lis = as.matrix(aL_table_gamma$branch.length), m = m_gamma, nt = n, q = c(0,1,2), reft = reft, cal = "PD")
#
#
#   aL_table_alpha = c()
#   N = 2
#   for (i in 1:N){
#
#     x = data_pair[data_pair[,i]>0,i]
#
#     aL = phyBranchAL_Abu(phylo = PDtree, data = x, rootExtend = T, refT = reft)
#     aL$treeNabu$branch.length = aL$BLbyT[,1]
#     suppressMessages({
#       aL_table = aL$treeNabu %>% select(branch.abun, branch.length, tgroup)
#     })
#
#     aL_table_alpha = rbind(aL_table_alpha, aL_table)
#
#   }
#
#
#   alpha = PhD.m.est(ai = aL_table_alpha$branch.abun, Lis = as.matrix(aL_table_alpha$branch.length), m = m_alpha, q = c(0,1,2), nt = n, reft = reft, cal = "PD")/N
#
#   # if (PDtype == 'meanPD') {
#     gamma = gamma/reft
#     alpha = alpha/reft
#   # }
#
#   beta = alpha
#   beta = gamma/alpha
#
#   # Compute dissimilarities
#   C_q0 <- (beta[1]^(1 - 0) - 1) / (2^(1 - 0) - 1)
#   C_q1 <- log(beta[2]) / log(2)
#   C_q2 <- (beta[3]^(1 - 2) - 1) / (2^(1 - 2) - 1)
#
#   # Return results as data.frame row
#   list(
#     site1 = pair_cols[1],
#     site2 = pair_cols[2],
#     sor_est = as.numeric(1 - C_q0),
#     hor_est = as.numeric(1 - C_q1),
#     mor_hor_est = as.numeric(1 - C_q2)
#   )
# }
