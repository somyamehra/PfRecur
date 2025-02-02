# PfRecur

This repository contains an R package to probabilistically classify late treatment failure for uncomplicated *Plasmodium falciparum* malaria. A detailed description of the underlying Bayesian model is provided in an accompanying preprint:

Mehra S, Taylor AR, Imwong M, White NJ and Watson JA, 2025. Probabilistic classification of late treatment failure in uncomplicated malaria. *medRxiv:2025.01.21.25320790*. Available at: https://www.medrxiv.org/content/10.1101/2025.01.21.25320790v1.

The main function `PfRecur::evaluate_posterior` can be used to generate a posterior summary for the number of newly-inoculated vs recrudescent clones within a recurrent *P. falciparum* infection based on bulk genotypic data, with user-defined parameters governing genotyping error and the detectability of individual clones.

The package largely relies on base R functionality, with external dependencies on the functions `copula::Stirling2`, `copula::Stirling1`, `PDQutils::cumulant2moment`, `PDQutils::moment2cumulant`, `poisbinom::dpoisbinom` and `VGAM::dbetabinom.ab`.

Please refer to the vignette for an application to genotypic data from a recent study conducted by [Dimbu et al (2024)](https://doi.org/10.1128/aac.01525-23).
