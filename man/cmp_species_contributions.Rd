% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/singleSpec_contributions.R
\name{cmp_species_contributions}
\alias{cmp_species_contributions}
\title{Evaluate species contributors for a single metabolite with qPCR/species abundance data}
\usage{
cmp_species_contributions(j, cmps_sub_good, all_rxns, subjects, norm_kos,
  ko_net, spec_abunds, ref_kos, cor_with = T)
}
\arguments{
\item{j}{metabolite # (usually from lapply/sapply)}

\item{cmps_sub_good}{CMP scores for metabolites with abundance data}

\item{all_rxns}{list of relevant reactions for each metabolite}

\item{subjects}{vector of subjects}

\item{norm_kos}{data.table of gene abundances}

\item{ko_net}{network created by generate_genomic_network}

\item{spec_abunds}{original species abundances}

\item{ref_kos}{gene abundances for each species}

\item{cor_with}{whether to look at the correlation of CMP scores of each species by itself with the metabolite, or of the whole community with that species removed}
}
\value{
List of length two. Item 1 is a table of all relevant species and correlation between the scores alone and the true scores, while item 2 is a vector of taxa that were classified as potential key taxa.
}
\description{
Evaluate species contributors for a single metabolite with qPCR/species abundance data
}
\examples{
sapply(1:length(metabolites), cmp_species_contributions, cmp_scores, all_rxns,
all_subjects, ko_abunds, ko_net, spec_abunds, ref_kos)
}
