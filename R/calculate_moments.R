#' Calculate moments for sample allele counts
#'
#' This function calculates moments of allele couns for the power set of a given
#' set.of alleles, summed over one or more baseline isolates
#'
#' @param recurrent_alleles Vector of alleles (strings); allele counts will be
#'  computed for the power set of these alleles
#' @param recurrent_MOI Maximum moment ordinal
#' @param baseline_genotypes A binary (0/1) matrix of genotypes for isolates
#' over which allele counts will be summed: named rows correspond to isolates;
#' named columns correspond to alleles; and matrix elements are set to 1 if the
#' relevant allele is called in a given isolate, and zero otherwise
#' @param baseline_MOI A named vector of MOIs for isolates over which allele
#' counts will be summed
#' @param error_matrix A row stochastic matrix (with row and column names
#' corresponding to alleles), where element (i,j) yields the probability
#' that an allele called as i in a baseline isolate matches an allele called
#' as j in the recurrent isolate
#' @return A matrix of moments of allele counts summed over baseline isolates for
#' subsets of recurrent_alleles: columns correspond to allelic subsets (governed by
#' inc_exclusion_mask); rows correspond to the moment ordinal (running from zero
#' to recurrent_MOI)
#'
locus_powerset_moments <- function(recurrent_alleles, recurrent_MOI, baseline_genotypes,
                                   baseline_MOI, error_matrix) {

  baseline_locus_card <- rowSums(baseline_genotypes)
  baseline_n_indiv <- length(baseline_locus_card)
  baseline_pop_MOI <- sum(baseline_MOI)

  # prob observed allele (rows) lies in subset of recurrent alleles (governed by
  # inc_exclusion_mask, each configuration encoded by a column)
  obs_allele_in_subset <-
    error_matrix[,recurrent_alleles] %*% inc_exclusion_mask[[length(recurrent_alleles)]]

  obs_allele_not_in_subset <-
    obs_allele_in_subset[, ncol(obs_allele_in_subset)]-obs_allele_in_subset

  # Poisson-binomial probability masses: list stratified by baseline individual,
  # each entry is a matrix: row_index+1 = h (prob mass at h)
  # column is allelic subset (governed by inc_exclusion_mask)
  pmf_allele_not_in_subset <-
    apply(baseline_genotypes, 1, function(x) {
      subset_probs <- obs_allele_not_in_subset[x>0,,drop=FALSE]
      apply(subset_probs, 2, function(pp) poisbinom::dpoisbinom(0:length(pp), pp))}, simplify = FALSE)

  # pmf for number of clones in each baseline individual with alleles in a given subset
  # list names: baseline individual, each entry is a matrix
  # column is allelic subset (governed by inc_exclusion mask), row is 0:MOI
  pmf_clones_not_in_subset <- lapply(1:baseline_n_indiv, function(i) {
    coeff_lookup[baseline_MOI[i], baseline_locus_card[i],
                 1:(baseline_MOI[i]+1), 1:(baseline_locus_card[i]+1)] %*%
      pmf_allele_not_in_subset[[i]]})

  # moments for number of clones in each baseline individual with alleles in a given subset
  # list names: baseline individual, each entry is a matrix
  # row is allelic subset (governed by inc_exclusion mask)
  # column is moment 1:recurrent_MOI
  moments_clones_not_in_subset <- lapply(1:baseline_n_indiv, function(i) {
    t(pmf_clones_not_in_subset[[i]]) %*% moment_matrix[[baseline_MOI[i]+1]][[recurrent_MOI]]
  })

  # cumulants for number of clones in each baseline individual with alleles in a given subset
  # list names: baseline individual, each entry is a matrix
  # column is allelic subset (governed by inc_exclusion mask)
  # row is cumulant 1:recurrent_MOI
  cumulants_clones_not_in_subset <- lapply(moments_clones_not_in_subset, function(x) {
    apply(x, 1, PDQutils::moment2cumulant)
  })

  # cumulants for number of clones with alleles not in subset summed over population
  # column is allelic subset (governed by inc_exclusion mask), row is cumulant
  pop_cumulants_not_in_subset <- Reduce(`+`, cumulants_clones_not_in_subset)

  if (recurrent_MOI==1) return(rbind(1, pop_cumulants_not_in_subset))

  # moments for number of clones with alleles not in subset summed over population
  # column is allelic subset (governed by inc_exclusion_mask), row is moment
  pop_moments_not_in_subset <- apply(pop_cumulants_not_in_subset, 2, PDQutils::cumulant2moment)

  return(rbind(1, pop_moments_not_in_subset))
}


#' Convert moments of sample allele counts to moments of population allele
#' frequencies derived under a multinomial-Dirichlet model (with a uniform
#' prior over the standard simplex)
#'
#' @param pop_moments_not_in_subset A matrix of moments of sample allele counts:
#' columns correspond to allelic subsets (governed by inclusion_exclusion_mask),
#' rows correspond to the moment ordinal (starting from zero)
#' @param recurrent_MOI Maximum moment ordinal
#' @param n_recurrent_alleles Number of distinct alleles over which allele counts
#' have been summed (the columns of pop_moments_not_in_subset iterate over the
#' relevant power set)
#' @param n_possible_alleles Number of possible alleles at the relevant locus
#' @param n_clones Number of clones over which allele counts have been computed
#' @return A matrix of moments of moments of sample allele frequencies
#'
dirichlet_transform <- function(pop_moments_not_in_subset, recurrent_MOI,
                                n_recurrent_alleles, n_possible_alleles, n_clones) {

  # correct previous error in derivation
  # subset_size <- n_possible_alleles - colSums(inc_exclusion_mask[[n_recurrent_alleles]])
  subset_size <- n_recurrent_alleles - colSums(inc_exclusion_mask[[n_recurrent_alleles]])

  calc_transform <- sapply(1:2^n_recurrent_alleles, function(i) {

    sapply(0:recurrent_MOI, function(m) {
      if (subset_size[i]==0) return(0)
      s1 <- sum(sapply(0:m, function(j) {
        sum(sapply(j:m, function(k) abs(copula::Stirling1(m, k))*choose(k, j)*subset_size[i]^(k-j)))*
          pop_moments_not_in_subset[j+1, i]}))
      s2 <- choose(n_clones+n_possible_alleles+m-1, n_clones+n_possible_alleles-1)*factorial(m)
      return(exp(log(s1)-log(s2)))
    })
  })

  return(calc_transform)
}


#' Convert moments of sample allele counts to moments of sample allele frequencies
#'
#' @param pop_moments_not_in_subset A matrix of moments of sample allele counts:
#' columns correspond to allelic subsets, rows correspond to the moment ordinal
#' (starting from zero)
#' @param recurrent_MOI Maximum moment ordinal
#' @param n_clones Number of clones over which allele counts have been computed
#' @return A matrix of moments of moments of sample allele frequencies
#'
sample_af_transform <- function(pop_moments_not_in_subset, recurrent_MOI, n_clones) {
  t(sapply(0:recurrent_MOI, function(i) pop_moments_not_in_subset[i+1,]/n_clones^i))
}
