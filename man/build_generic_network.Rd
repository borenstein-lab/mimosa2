% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/speciesPicrustFunctions.R
\name{build_generic_network}
\alias{build_generic_network}
\title{Function called by run_pipeline to get species-specific reaction network}
\usage{
build_generic_network(contribution_table, kegg_paths)
}
\arguments{
\item{contribution_table}{OTU/seq table}

\item{kegg_paths}{paths to KEGG network files}
}
\value{
Table of species-specific KEGG reactions
}
\description{
Function called by run_pipeline to get species-specific reaction network
}
\examples{
build_generic_network(species_table, "Greengenes 13_5 or 13_8", picrust_paths, kegg_paths)
}
