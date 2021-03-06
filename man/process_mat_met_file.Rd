% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/taxa_agora.R
\name{process_mat_met_file}
\alias{process_mat_met_file}
\title{Returns the column from the picrust tables that corresponds to the genomic content of the indicated OTU}
\usage{
process_mat_met_file(mat_met_file = "BiKEGG-master/AllKEGG2BiGGmet.mat")
}
\arguments{
\item{mat_met_file}{Matlab file containing mapping of AGORA metabolite IDs to KEGG}
}
\value{
data.table of KEGG/AGORA metabolite mapping
}
\description{
Returns the column from the picrust tables that corresponds to the genomic content of the indicated OTU
}
\examples{
met_id_table = process_mat_met_file(mat_met_file)
}
