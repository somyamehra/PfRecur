#' Prior distribution
#'
#' This function computes the prior distribution of recrudescent (vs
#' newly-inoculated) clones within a recurrent isolate (symmemetric
#' beta binomial)
#'
#' @param n_C Number of recrudescent clones
#' @param recurrent_MOI MOI of recurrent isolate
#' @param beta Shape parameter (in the range (0, 1]) for the symmetric beta binomial
#' prior (with support 0 to recurrent_MOI)
prior_mixture <- function(n_C, recurrent_MOI, beta) {
  VGAM::dbetabinom.ab(n_C, recurrent_MOI, shape1=beta, shape2=beta)
}

#' Posterior summary
#'
#' This function computes the posterior distribution of recrudescent vs
#' newly-inoculated clones within a recurrent isolate, handling isolates
#' with MOI 9 and fewer
#'
#' @param recurrent_isolate Name (string) of the recurrent isolate for which classification
#' is to be performed
#' @param ref_C Name (string) of the baseline isolate paired to the recurrent isolate
#' @param ref_I Vector of names (strings) of baseline isolates that are NOT
#' paired to the recurrent isolate, but used to estimate allele frequencies for
#' newly-inoculated clones
#' @param genotype_matrix Named marker-wise list where each element is a binary (0/1)
#' matrix of genotypes at a given marker: named rows correspond to isolates; named columns
#' correspond to alleles; and matrix elements are set to 1 if the relevant allele
#' is called in a given isolate, and zero otherwise. Every row must have sum at
#' least one; isolates with no alleles called at a given marker are dropped from the
#' genotype matrix for the relevant marker.
#' @param error_matrix Named marker-wise list where each element is a row stochastic
#' matrix (with row and column names corresponding to alleles at a given marker),
#' where element (i,j) yields the probability that an allele called as i in a
#' baseline isolate matches an allele called as j in the recurrent isolate at that marker
#' @param omega_vals Vector of values for omega (in the range (0, 1]), governing the locus-wise per-clone
#' probability of detection (applied to the paired baseline and recurrent isolates only)
#' @param beta Shape parameter (in the range (0, 1]) for the symmetric beta binomial prior(with support 0 to
#' recurrent_MOI) for the number of recrudescent clones in the recurrent isolate
#' @param keep_markers Set to "all" (by default) to calculate the posterior over
#' all markers; otherwise, a vector of markers (strings) over which to calculate
#' the posterior
#' @return A list of two elements: element "posterior" is a data frame for the
#' posterior distribution (posterior) of newly-inoculated (n_I) vs recrudescent (n_C)
#' clones in the recurrent isolate, stratified by the per-clone probability of
#' detection for the paired baseline/recurrent isolates (omega); element "metrics"
#' is a data frame detailing the posterior probability of at least one recrudescent
#' clone (M1) and the expected proprtion of recrudescent clones (M2) stratified by omega
#' @export
evaluate_posterior <- function(recurrent_isolate, ref_C, ref_I, genotype_matrix,
                               error_matrix, omega_vals=seq(0.75, 1, 0.05), beta, keep_markers="all") {

  if (beta<=0 | beta>1) stop("beta must lie in the range (0, 1]")

  if (any(omega_vals<=0) | any(omega_vals>1)) stop("omega must lie in the range (0, 1]")


  if (is.null(names(genotype_matrix)) | !is.list(genotype_matrix) |
      !all(sapply(genotype_matrix, is.matrix)) |
      !all(sapply(genotype_matrix, function(x) !is.null(rownames(x)))) |
      !all(sapply(genotype_matrix, function(x) !is.null(colnames(x)))) |
      !all(sapply(genotype_matrix, function(x) all(x==0 | x==1)))) {
    stop("genotype_matrix must be a named list of binary (0/1) matrices")
  }

  if (is.null(names(error_matrix)) | !is.list(error_matrix) |
      !all(sapply(error_matrix, is.matrix)) |
      !all(sapply(error_matrix, function(x) !is.null(rownames(x)))) |
      !all(sapply(error_matrix, function(x) !is.null(colnames(x)))) |
      !all(sapply(error_matrix, function(x) all(x>=0))) |
      !all(sapply(error_matrix, function(x) all(sort(rownames(x))==sort(colnames(x))))) |
      !all(round(unlist(sapply(error_matrix, rowSums, simplify = FALSE)), 5)==1)) {
    stop("error_matrix must be a named list of square row stochastic matrices (i.e. non-negative entries, row sums 1)")
  }

  if (!all(names(genotype_matrix) %in% names(error_matrix))) {
    stop("genotype_matrix contains markers than are not present in error_matrix")
  }

  if (!all(mapply(function(x, y) all(sort(x)==sort(y)),
                  lapply(genotype_matrix[keep_markers], colnames),
                  lapply(error_matrix[keep_markers], rownames)))) {
    stop("genotype_matrix contains alleles (at one or more marker) that are not present in error_matrix")
  }

  if (length(keep_markers)==1) {
    if (keep_markers=="all") keep_markers <- names(genotype_matrix)
  }

  if (!all(keep_markers %in% names(genotype_matrix))) {
    stop("keep_markers contains markers that are not present in genotype_matrix")
  }

  # ignore missing markers
  missing_markers <- sapply(keep_markers, function(m)
    !(recurrent_isolate %in% rownames(genotype_matrix[[m]])) |
      !(ref_C %in% rownames(genotype_matrix[[m]])) |
      !any(ref_I %in% rownames(genotype_matrix[[m]])))

  keep_markers <- keep_markers[!missing_markers]

  # incomparable due to missing genotypes
  if (length(keep_markers)==0) return(NULL)

  recurrent_alleles <- lapply(keep_markers, function(m)
    names(which(genotype_matrix[[m]][recurrent_isolate,]==1)))
  names(recurrent_alleles) <- keep_markers

  recurrent_MOI <- max(sapply(recurrent_alleles, length))

  ref_I_genotypes <- lapply(genotype_matrix[keep_markers], function(x)
    x[sort(intersect(rownames(x), ref_I)),,drop=FALSE])

  ref_I_MOI <- sapply(ref_I_genotypes, rowSums, simplify = FALSE) %>%
    bind_rows() %>% apply(2, max, na.rm=T)

  ref_C_genotypes <- lapply(genotype_matrix[keep_markers], function(x)
    x[sort(intersect(rownames(x), ref_C)),,drop=FALSE])

  ref_C_MOI <- max(sapply(ref_C_genotypes, rowSums))

  if (recurrent_MOI>9 | ref_C_MOI>9 | any(ref_I_MOI>9)) {
    stop("All isolate MOIs must be 9 or fewer")
  }

  if (min(sapply(ref_C_genotypes, function(x) min(rowSums(x))))==0 |
      min(sapply(ref_I_genotypes, function(x) min(rowSums(x))))==0 |
      min(sapply(recurrent_alleles, length))==0) {
    stop("Every row of the marker-wise genotyping matrix must have at least one non-zero entry")
  }

  posterior_distribution <- lapply(keep_markers, function(m)
    locus_likelihood(recurrent_alleles[[m]], recurrent_MOI, ref_I_genotypes[[m]],
                     ref_I_MOI, ref_C_genotypes[[m]], ref_C_MOI,
                     error_matrix[[m]], omega_vals)) %>%
    bind_rows(.id="marker") %>%
    group_by(omega, n_I, n_C) %>%
    summarise(lhood=prod(lhood)) %>%
    mutate(prior=prior_mixture(n_C, recurrent_MOI, beta), prob=lhood*prior) %>%
    split(f=.$omega) %>% lapply(function(x) {
      x$posterior=x$prob/sum(x$prob); return(x)}) %>%
    bind_rows() %>% select(omega, n_I, n_C, posterior) %>%
    mutate(recurrence=recurrent_isolate, posterior=round(posterior, 6))

  posterior_metrics <- posterior_distribution %>% group_by(recurrence, omega) %>%
    summarise(M1=sum((n_C>0)*posterior), M2=sum(n_C*posterior)/recurrent_MOI)

  return(list(posterior=posterior_distribution, metrics=posterior_metrics))
}
