#' Locuswise likelihood conditional on detected clones
#' 
#' This function evaluates the locus-wise likelihood of observations for the recurrent 
#' isolate conditional on the number of successfully genotyped clones in the paired
#' baseline and recurrent isolates
#'
#' @param recurrent_alleles Vector of alleles (strings) observed in the recurrent isolate
#' @param recurrent_MOI MOI of the recurrent isolate
#' @param ref_I_genotypes A binary (0/1) matrix of genotypes for baseline isolates 
#' over which population I allele frequencies are derived: named rows correspond to isolates; 
#' named columns correspond to alleles; and matrix elements are set to 1 if the 
#' relevant allele is called in a given isolate, and zero otherwise
#' @param ref_I_MOI A named vector of MOIs for baseline isolates over which 
#' population I allele frequencies are derived
#' @param ref_C_genotypes A binary (0/1) single-row matrix of genotypes for the
#' paired baseline isolates: named rows correspond to isolates: named columns 
#' correspond to alleles; and matrix elements are set to 1 if the relevant allele 
#' is called in a given isolate, and zero otherwise
#' @param ref_C_MOI MOI for the paired baseline isolate
#' @param error_matrix A row stochastic matrix (with row and column names 
#' corresponding to alleles), where element (i,j) yields the probability
#' that an allele called as i in a baseline isolate matches an allele called
#' as j in the recurrent isolate
#' @return A data frame for the conditional locus-wise likelihood (p) of the observed
#' genotypes for the recurrent isolate, given the number of newly-inoculated (n_I)
#' and recrudescent (n_C) clones comprising the recurrent isolate; and the number
#' of succesfully genotyped clones in the paired baseline (genotyped_C) and
#' recurrent (genotyped_recurrent) isolates respectively
#'           
locus_likelihood_given_detected_clones <- 
  function(recurrent_alleles, recurrent_MOI, ref_I_genotypes, ref_I_MOI, 
           ref_C_genotypes, ref_C_MOI, error_matrix) {
  
  # calculate moments for population I allele frequencies
  ref_I_freq <- locus_powerset_moments(recurrent_alleles, recurrent_MOI, ref_I_genotypes,
                                       ref_I_MOI, error_matrix)
  
  ref_I_freq <- dirichlet_transform(ref_I_freq, recurrent_MOI, length(recurrent_alleles),
                                    nrow(error_matrix), sum(ref_I_MOI))
  
  # calculate moments for allele frequencies over successfully genotyped
  # clones in the paired baseline islolate
  ref_C_freq <- lapply(sum(ref_C_genotypes):ref_C_MOI, function(nonlatent) {
    ref_C_mom <- locus_powerset_moments(recurrent_alleles, recurrent_MOI, ref_C_genotypes,
                                        nonlatent, error_matrix)
    return(sample_af_transform(ref_C_mom, recurrent_MOI, nonlatent))})
  
  # locuswise likelihood conditional on the number of successfully genotyped
  # clones in the paired baseline isolate (ref_C) and recurrent isolate
  mixture_prob <- mapply(function(genotyped_C, genotyped_C_freq) {
    y <- expand.grid(n_I=0:recurrent_MOI, n_C=recurrent_MOI:0, genotyped_C=genotyped_C) %>% 
      mutate(genotyped_recurrent=n_I+n_C) %>% 
      subset(genotyped_recurrent>=length(recurrent_alleles) & 
               genotyped_recurrent<=recurrent_MOI & n_C<=genotyped_C)
    y$p <- apply(y, 1, function(x) sum(ref_I_freq[x[1]+1,]*genotyped_C_freq[x[2]+1,]*
                                         inc_exclusion_sign[[length(recurrent_alleles)]]))
    return(y)},
    sum(ref_C_genotypes):ref_C_MOI, ref_C_freq, SIMPLIFY = FALSE) %>% bind_rows()
  
  return(mixture_prob)
}


#' Convolve the conditional likelihood over the distribution of detected clones
#' 
#' This function evaluates convolves the conditional locus-wise likelihood of 
#' observations for the recurrent isolate, given the number of successfully genotyped 
#' clones in the paired baseline and recurrent isolates, over the distribution
#' of detected clones
#'
#' @param cond_lhood Data frame of conditional likelihoods generated with the 
#' function locus_likelihood_given_detected_clones: this encodes the conditional 
#' locus-wise likelihood (p) of the observed genotypes for the recurrent isolate, 
#' given the number of newly-inoculated (n_I) and recrudescent (n_C) clones 
#' comprising the recurrent isolate; and the number of succesfully genotyped clones 
#' in the paired baseline (genotyped_C) and recurrent (genotyped_recurrent) isolates 
#' respectively
#' @param recurrent_MOI MOI of the recurrent isolate
#' @param recurrent_locus_card cardinality of the recurrent locus at the given marker
#' @param ref_C_MOI MOI for the paired baseline isolate
#' @param omega Locus-wise per-clone probability of detection (applied to the
#' paired baseline and recurrent isolates only)
#' @return A data frame for the locus-wise likelihood (lhood) of the observed
#' genotypes for the recurrent isolate, stratified by the number of newly-inoculated (n_I)
#' vs recrudescent (n_C) clones
#'           
convolve_over_detected_clones <- function(omega, cond_lhood, recurrent_locus_card, recurrent_MOI, ref_C_MOI) {
  
  # convolve over the distribution of detected clones in the paired
  # baseline isolate (locus-wise truncated binomial model)
  convolve_over_paired_baseline <- cond_lhood %>% split(f=.$genotyped_recurrent) %>%
    lapply(function(x) {
      min_ref_C_MOI <- min(x$genotyped_C)
      genotyped_recurrent <- x[1, "genotyped_recurrent"]
      
      mix_prob <- 
        data.frame(genotyped_recurrent=genotyped_recurrent,
                   n_C=0:ref_C_MOI, n_I=genotyped_recurrent:(genotyped_recurrent-ref_C_MOI)) %>%
        subset(n_I>=0)
      
      mix_prob$lhood <- sapply(mix_prob$n_C, function(n_recrudescent) {
        mix_lhood <- x %>% subset(n_C<=n_recrudescent) %>% 
          mutate(p=p*dbinom(n_C, n_recrudescent, genotyped_C/ref_C_MOI)*
                   dbinom(genotyped_C, ref_C_MOI, omega)/
                   (1-pbinom(min_ref_C_MOI-1, ref_C_MOI, omega)))
        return(sum(mix_lhood$p))})
      
      return(mix_prob)
    }) %>% bind_rows()
  
  
  # convolve over the distribution of detected clones in the recurrent
  # isolate (locus-wise truncated multinomial model)
  recurrent_prob <- data.frame(n_C=0:min(ref_C_MOI, recurrent_MOI), omega=omega) %>%
    mutate(n_I=recurrent_MOI-n_C)
  recurrent_prob$lhood <- sapply(recurrent_prob$n_C, function(true_n_C) 
    sum(dbinom(convolve_over_paired_baseline$n_C, true_n_C, omega)*
          dbinom(convolve_over_paired_baseline$n_I, recurrent_MOI-true_n_C, omega)*
          convolve_over_paired_baseline$lhood))/
    (1-pbinom(recurrent_locus_card-1, recurrent_MOI, omega))
  
  
  return(recurrent_prob)
}

#' Locuswise likelihood
#' 
#' This function evaluates the locus-wise likelihood of observations for the 
#' recurrent isolate
#'
#' @param recurrent_alleles Vector of alleles (strings) observed in the recurrent isolate
#' @param recurrent_MOI MOI of the recurrent isolate
#' @param ref_I_genotypes A binary (0/1) matrix of genotypes for baseline isolates 
#' over which population I allele frequencies are derived: named rows correspond to isolates; 
#' named columns correspond to alleles; and matrix elements are set to 1 if the 
#' relevant allele is called in a given isolate, and zero otherwise
#' @param ref_I_MOI A named vector of MOIs for baseline isolates over which 
#' population I allele frequencies are derived
#' @param ref_C_genotypes A binary (0/1) single-row matrix of genotypes for the
#' paired baseline isolates: named rows correspond to isolates: named columns 
#' correspond to alleles; and matrix elements are set to 1 if the relevant allele 
#' is called in a given isolate, and zero otherwise
#' @param ref_C_MOI MOI for the paired baseline isolate
#' @param error_matrix A row stochastic matrix (with row and column names 
#' corresponding to alleles), where element (i,j) yields the probability
#' that an allele called as i in a baseline isolate matches an allele called
#' as j in the recurrent isolate
#' @param omega_vals Vector of values for omega, governing the locus-wise per-clone 
#' probability of detection (applied to the paired baseline and recurrent isolates only)
#' @return A data frame for the locus-wise likelihood (lhood) of the observed
#' genotypes for the recurrent isolate, given the number of newly-inoculated (n_I)
#' and recrudescent (n_C) clones comprising the recurrent isolate; and the per-clone
#' probability of detection for the paired baseline/recurrent isolates (omega)
#'  
locus_likelihood <- function(recurrent_alleles, recurrent_MOI, ref_I_genotypes, ref_I_MOI, 
                             ref_C_genotypes, ref_C_MOI, error_matrix, omega_vals) {
  
  cond_lhood_detected_clones <- 
    locus_likelihood_given_detected_clones(recurrent_alleles, recurrent_MOI, ref_I_genotypes, 
                                           ref_I_MOI, ref_C_genotypes, ref_C_MOI, error_matrix)
  
  lhood <- sapply(omega_vals, convolve_over_detected_clones, cond_lhood_detected_clones, 
                  length(recurrent_alleles), recurrent_MOI, ref_C_MOI, simplify = FALSE)
  
  return(bind_rows(lhood))
  
}