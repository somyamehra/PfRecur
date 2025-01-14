# PfRecur

This repository contains an R package to probabilistically classify late treatment failure for uncomplicated *Plasmodium falciparum* malaria. A detailed description of the underlying Bayesian model is provided in an accompanying preprint:

Mehra S, Taylor AR, Imwong M, White NJ, Watson JA. Probabilistic classification of late treatment failure in uncomplicated malaria. 2025.

The main function `PfRecur:: evaluate_posterior` can be used to generate a posterior summary for the number of newly-inoculated vs recrudescent clones within a recurrent *P. falciparum* infection based on bulk genotypic data, with user-defined parameters governing genotyping error and the detectability of individual clones.

Please refer to the vignette for an application to genotypic data from a recent study conducted by [Dimbu et al (2024)](https://doi.org/10.1128/aac.01525-23).
