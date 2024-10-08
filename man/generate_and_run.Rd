% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_fits.R
\name{generate_and_run}
\alias{generate_and_run}
\title{Generate synthetic data and run the fit}
\usage{
generate_and_run(
  comb_matrix,
  py,
  fits_path = NULL,
  data_path = NULL,
  seeds = 1:30,
  catalogue_sbs,
  alpha_range = c(0.15, 0.2),
  alpha_sigma = 0.1,
  pi_conc = 1,
  frac_rare = 1,
  n_muts_range = 100:5000,
  shared_sbs = c("SBS1", "SBS5"),
  reference_catalogue = COSMIC_filt,
  subset_reference = c("SBS1", "SBS5"),
  keep_sigs = c("SBS1", "SBS5"),
  hyperparameters = NULL,
  lr = 0.005,
  n_steps = 3000,
  nonparametric = TRUE,
  enforce_sparsity = TRUE,
  reg_weight = 0,
  regularizer = "cosine",
  regul_denovo = TRUE,
  store_fits = FALSE,
  seed_list = c(4, 17, 22),
  CUDA = TRUE,
  do.fits = FALSE,
  cohort = "",
  check_present = TRUE,
  check_linear_comb = FALSE,
  inference_type = c("flat", "hier", "clust"),
  ...
)
}
\arguments{
\item{comb_matrix}{Combinations of parameters to test}

\item{py}{Either \code{NULL} or the \code{Python} package to use}

\item{fits_path}{Path to store the fits}

\item{data_path}{Path to store the data}

\item{seeds}{List of seeds to use}

\item{reference_catalogue}{Reference catalogue to use for the inference}

\item{reg_weight}{Regularization weight}

\item{regularizer}{Type of regularization}

\item{CUDA}{Logical as whether run it on GPU or not}

\item{do.fits}{Logical as whether to run the fits}

\item{shared}{Signatures shared across groups}

\item{private}{Signatures private to one group}

\item{catalogue}{Reference catalogue}

\item{mut_range}{Number of mutation range for each patient}

\item{input_catalogue}{Input catalogue to use for the fit}

\item{verbose}{Logical as whether show the progression}

\item{new_model}{Logical as whether to run the new model}
}
\value{
nothing
}
\description{
Generate synthetic data and run the fit
}
