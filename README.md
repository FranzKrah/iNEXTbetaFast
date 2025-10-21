# iNEXTbetaFast

This package uses code from the R package iNEXT.beta3D to compute pairwise dissimilarity matrices along the Hill series.

# Note
The parallelization currently utilizes forking with the parallel and mclapply R packages and thus 
run on UNIX based systems only. The number of cores to use depends on the size of the task.
I recommend starting with half the cores available and checking closely how your memory usage 
behaves (Mac: Activity Monitor, Ubuntu: System Monitor). If memory usage saturates soon then you can probably increase the number of cores carefully.
If the memory used increases rapidly towards 100% you should stop the run and reduce the number of cores.
From my experience on mac N-2 is fine, on Ubuntu N-6 is ok but it all depends on your specific system so you need to try a bit to optimize speed.

## Installation

```r
install.packages("devtools")
devtools:install_github("FranzKrah/iNEXTbetaFast")
```

and then use it by typing ... 

```r
library(iNEXTbetafast)
```

## Example
```r
# Example dataset
library(microeco)
data(otu_table_ITS)

# Lowest pairwise observed sample coverage
Cmax_joint <- Cmax_TD_fast(otu_table_ITS, ncores = 2)

# Estimate
beta_td <- run_pairwise(otu_table_ITS, pair_c_a_td, level = Cmax_joint, ncores = 10)

# Distance matrix
d_q0 <- pair2dist(beta_td, colname = "sor_est")
d_q1 <- pair2dist(beta_td, colname = "hor_est")
d_q2 <- pair2dist(beta_td, colname = "mor_hor_est")

```

