#' Genotypic data from a TES conducted in Angola by Dimbu et al (2024)
#' 
#' This contains genotypes at a panel of 7 neutral microsatellite markers
#' (M313, M383, TA1, POLYA, PFPK2, M2490, TA109) for 182 isolates sampled
#' across 3 study sites (Benguela, Lunda Sul and Zaire) in a TES conducted
#' in Angola by Dimbu et al (2024)
#'
#' A parsed version of Supplemental Table S3 from the paper:
#' Dimbu PR, Labuda S, Ferreira CM, Caquece F, André K, Pembele G, Pode D, João MF, 
#' Pelenda VM, Nieto Andrade B, Horton B. Therapeutic response to four artemisinin-based 
#' combination therapies in Angola, 2021. Antimicrobial Agents and Chemotherapy. 
#' 2024 Apr 3;68(4):e01525-23.
#'
#' @format A named list of binary (0/1) genotype matrices: list names correspond
#' to markers; matrix row names correspond to isolates; matrix column names
#' correspond to alleles; and matrix elements are set to 1 if the relevant allele
#' is called in a given isolate at that marker, and zero otherwise.
#' 
#' @source <https://doi.org/10.1128/aac.01525-23>
"Angola_TES_genotype_matrix"

#' Isolate groupings for a TES conducted in Angola by Dimbu et al (2024)
#' 
#' For each recurrent isolate genotyped by Dimbu et al (2024) 
#' in a TES conducted in Angola this specifies the paired baseline isolate 
#' (i.e., the day 0 isolate from the same study participant); and the set of 
#' unpaired baseline isolates (i.e., day 0 isolates from different study participants
#' but the same study site) that are used to derive allele frequencies for
#' newly-inoculated clones
#'
#' Derived from Supplemental Table S3 of the paper:
#' Dimbu PR, Labuda S, Ferreira CM, Caquece F, André K, Pembele G, Pode D, João MF, 
#' Pelenda VM, Nieto Andrade B, Horton B. Therapeutic response to four artemisinin-based 
#' combination therapies in Angola, 2021. Antimicrobial Agents and Chemotherapy. 
#' 2024 Apr 3;68(4):e01525-23.
#' 
#' @format A list of 70 lists (each corresponding to a recurrent isolate), where for
#' each list:
#' \describe{
#'   \item{recurrent}{Identifier (string) for the recurrent isolate}
#'   \item{ref_C}{Identifier (string) for the baseline isolate paired to the recurrent isolate}
#'   \item{ref_I}{Vector of identifiers (strings) for baseline isolates that are not paired
#'   to the recurrent isolate, but from the same study site}
#' }
#' 
#' @source <https://doi.org/10.1128/aac.01525-23>
"Angola_TES_isolates"

