#' Normalised geometric genotyping error model
#'
#' Length-dependent genotyping error, adapted to length polymorphic markers with 
#' integer lengths: genotyping error rate EPSILON, with a geometric dependence
#' on the allele length difference
#'
#' @param marker_set A named marker-wise list of vectors, listing all possible
#' alleles at each marker (must be strings, but encoding integers)
#' @param EPSILON Genotyping error rate 
#' @return A named marker-wise list of row stochastic matrices (with row and column names 
#' corresponding to alleles), where element (i,j) yields the probability
#' that an allele called as i in a baseline isolate matches an allele called
#' as j in the recurrent isolate
#' @export
geometric_error_matrix <- function(marker_set, EPSILON) {
  
  if (!all(sapply(marker_set, function(x) all(!grepl("[^0-9]", x))))) {
    print("This error matrix can only be generated for length-polymorphic markers with integer allele lengths")
    return(NULL)
  } 
  
  delta <- lapply(marker_set, function(markers) {
    
    n_markers <- length(markers)
    repeat_length <- as.numeric(markers)
    
    error_matrix <- matrix(1, nrow=n_markers, ncol=n_markers,
                           dimnames = list(markers, markers))
    
    if (n_markers==1) return(error_matrix)
    
    diag(error_matrix) <- 1-EPSILON
    
    if (EPSILON>0) {
      for (i in 1:n_markers) {
        non_self <- setdiff(repeat_length, repeat_length[i])
        p <- dgeom(abs(non_self-repeat_length[i]), 1-EPSILON)
        error_matrix[i, as.character(non_self)] <- p*EPSILON/sum(p)
      }
    } else {
      error_matrix <- 0
      diag(error_matrix) <- 1
    }
    return(error_matrix) }
  )
  
  return(delta)
}

#' Uniform error model
#'
#' Adapted to nominal marker data: genotyping error rate EPSILON; distributed
#' uniformly over all possible alleles at the marker
#'
#' @param marker_set A named marker-wise list of vectors, listing all possible
#' alleles (strings) at each marker
#' @param EPSILON Genotyping error rate 
#' @return A named marker-wise list of row stochastic matrices (with row and column names 
#' corresponding to alleles), where element (i,j) yields the probability
#' that an allele called as i in a baseline isolate matches an allele called
#' as j in the recurrent isolate
#' @export
uniform_error_matrix <- function(marker_set, EPSILON) {
  
  if (!all(sapply(marker_set, function(x) all(!grepl("[^0-9]", x))))) {
    print("This error matrix can only be generated for length-polymorphic markers with integer allele lengths")
    return(NULL)
  } 
  
  delta <- lapply(marker_set, function(markers) {
    
    n_markers <- length(markers)

    error_matrix <- matrix(1, nrow=n_markers, ncol=n_markers,
                           dimnames = list(markers, markers))
    
    if (n_markers==1) return(error_matrix)
    
    error_matrix <- EPSILON/(n_markers-1)
    diag(error_matrix) <- 1-EPSILON
    
    return(error_matrix) 
  })
  
  return(delta)
}
