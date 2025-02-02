% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{Angola_TES_isolates}
\alias{Angola_TES_isolates}
\title{Isolate groupings for a TES conducted in Angola by Dimbu et al (2024)}
\format{
A list of 70 lists (each corresponding to a recurrent isolate), where for
each list:
\describe{
\item{recurrent}{Identifier (string) for the recurrent isolate}
\item{ref_C}{Identifier (string) for the baseline isolate paired to the recurrent isolate}
\item{ref_I}{Vector of identifiers (strings) for baseline isolates that are not paired
to the recurrent isolate, but from the same study site}
}
}
\source{
\url{https://github.com/MateuszPlucinski/AngolaTES2021}
}
\usage{
Angola_TES_isolates
}
\description{
For each recurrent isolate genotyped by Dimbu et al (2024)
in a TES conducted in Angola this specifies the paired baseline isolate
(i.e., the day 0 isolate from the same study participant); and the set of
unpaired baseline isolates (i.e., day 0 isolates from different study participants
but the same study site) that are used to derive allele frequencies for
newly-inoculated clones
}
\details{
Dimbu PR, Labuda S, Ferreira CM, Caquece F, André K, Pembele G, Pode D, João MF,
Pelenda VM, Nieto Andrade B, Horton B. Therapeutic response to four artemisinin-based
combination therapies in Angola, 2021. Antimicrobial Agents and Chemotherapy.
2024 Apr 3;68(4):e01525-23.
}
\keyword{datasets}
