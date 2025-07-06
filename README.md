# PfRecur

This repository contains an R package to probabilistically classify late treatment failure for uncomplicated *Plasmodium falciparum* malaria. A detailed description of the underlying Bayesian model is provided in an accompanying preprint:

Mehra S, Taylor AR, Imwong M, White NJ and Watson JA, 2025. Probabilistic classification of late treatment failure in uncomplicated malaria. *medRxiv:2025.01.21.25320790*. Available at: https://www.medrxiv.org/content/10.1101/2025.01.21.25320790v1.

The main function `PfRecur::evaluate_posterior` can be used to generate a posterior summary for the number of newly-inoculated vs recrudescent clones within a recurrent *P. falciparum* infection based on bulk genotypic data, with user-defined parameters governing genotyping error and the detectability of individual clones.

## Example dataset

As an example, the package includes genotypic data from a recent study conducted by Dimbu et al (2024), and retrieved from an accompanying GitHub repository, in the lists `PfRecur::Angola_TES_genotype_matrix` and `PfRecur::Angola_TES_isolates`.

**References**

Dimbu PR, Labuda S, Ferreira CM, Caquece F, André K, Pembele G, Pode D, João MF, Pelenda VM, Nieto Andrade B, Horton B. Therapeutic response to four artemisinin-based combination therapies in Angola, 2021. Antimicrobial Agents and Chemotherapy. 2024 Apr 3;68(4):e01525-23. Available at: <https://journals.asm.org/doi/full/10.1128/aac.01525-23>

Plucinski M. *AngolaTES2021*. GitHub. 2024. <https://github.com/MateuszPlucinski/AngolaTES2021>

## Dependencies

The package largely relies on base R functionality, with external dependencies on the functions `copula::Stirling2`, `copula::Stirling1`, `PDQutils::cumulant2moment`, `PDQutils::moment2cumulant`, `poisbinom::dpoisbinom` and `VGAM::dbetabinom.ab`.

We have tested the package with copula V1.1-1, PDQutils V0.1.6, poisbinom V1.0.1 and VGAM V1.1-9.

## Installation

The package can be installed using the command:

`devtools::install_github("somyamehra/PfRecur")`

Installation should take under a minute on a standard laptop.

## Vignette

Please refer to the vignette for instructions to the run the package and an application to genotypic data from the study of Dimbu et al (2024).

Running the function `PfRecur::evaluate_posterior` for the 70 paired recurrences in the dataset of Dimbu et al (2024) should take under 10 seconds on a standard laptop for a given parameter set (i.e., a user-defined genotyping error matrix and per-clone detection probability). The runtime may be significantly higher for datasets with elevated multiplicities of infection.

To analyse larger datasets or perform sensitivity analyses for user-defined parameters, it may be useful to parallelise classification across paired recurrences. We provide an example of this in the vignette, using the functions `parallel::mclapply` and `purrr::pmap`.
