library(copula)

MAX_MOI <- 9

# inclusion-exclusion principle vectors
inc_exclusion_mask <- lapply(1:MAX_MOI, function(n) t(expand.grid(rep(list(c(0,1)), n))))
inc_exclusion_sign <- lapply(inc_exclusion_mask, function(x) (-1)^colSums(x))

# coefficient look up table
# create a lookup table indexed MOI, locus_card, num_clones, overlap_card
coeff_lookup <- array(0, dim=rep(MAX_MOI, 4)+c(0,0,1,1))
for (MOI in 1:MAX_MOI) {
  for (locus_card in 1:MOI) {
    for (num_clones in 0:MOI) {
      min_index <- max(0, locus_card-MOI+num_clones)
      max_index <- min(num_clones, locus_card)
      for (overlap_card in min_index:max_index) {
        coeff_lookup[MOI, locus_card, num_clones+1, overlap_card+1] <-
          (choose(MOI, num_clones)/choose(locus_card, overlap_card))*
          (Stirling2(num_clones, overlap_card)/Stirling2(MOI, locus_card)*
             Stirling2(MOI-num_clones, locus_card-overlap_card))
      }
    }
  }
}

moment_matrix <- lapply(0:MAX_MOI, function(n) {
  lapply(1:MAX_MOI, function(m) {
    mm <- matrix(0, nrow=n+1, ncol=m)
    for (i_n in 0:n) {
      for (i_m in 1:m) {
        mm[i_n+1, i_m] <- i_n^i_m
      }
    }
    return(mm)})
})

usethis::use_data(MAX_MOI, inc_exclusion_mask, inc_exclusion_sign, coeff_lookup, 
                  moment_matrix, internal = TRUE, overwrite = TRUE)