### PICRUSt-single species model functions - modified from BURRITO

#' Set column name of OTU table to "OTU"
#'
#' @import data.table
#' @param species_table
#' @return Species table with column name as OTU
#' @examples
#' spec_table_fix(spec_table)
#' @export
spec_table_fix = function(species_table){
  otu_col = names(species_table)[names(species_table) %in% c("OTU", "#OTU ID", "Sequence", "ASV")]
  if(length(otu_col) != 1) stop("Species data must have a valid OTU/sequence column name")
  setnames(species_table, otu_col, "OTU")
  return(species_table)
}




#' Function called by run_pipeline to get species-specific KOs and reactions
#'
#' @import data.table
#' @param species_table OTU/seq table
#' @param database ESV/GReengenes/SILVA
#' @param picrust_paths paths to PICRUSt info files
#' @param kegg_paths paths to KEGG network files
#' @return Table of species-specific KEGG reactions
#' @examples
#' build_generic_network(species_table, "Greengenes 13_5 or 13_8", picrust_paths, kegg_paths)
#' @export
build_generic_network = function(species_table, database, picrust_paths, kegg_paths, repSeqPath = NA){
  if(database=="Sequence variants (recommended for AGORA)"){
    seq_list = species_table[,OTU]
    if(!is.na(repSeqPath)) stop("Path to representative sequences required!")
    species_table = get_otus_from_seqvar(species_table[,OTU], ) #Run vsearch to get gg OTUs
  } else if(database != "Greengenes 13_5 or 13_8"){
    stop("Only Greengenes currently implemented")
  }
  contribution_table = generate_contribution_table_using_picrust(species_table, picrust_norm_file = picrust_paths[1], picrust_ko_table_directory = picrust_paths[2], picrust_ko_table_suffix = picrust_paths[3])
  contribution_table = contribution_table[contribution != 0]
  all_kegg = get_kegg_reaction_info(kegg_paths[2], reaction_info_file = kegg_paths[3], save_out = F)
  network_template = generate_network_template_kegg(kegg_paths[1], all_kegg = all_kegg, write_out = F) #We should really speed this thing up
  otu_list = contribution_table[,sort(unique(OTU))]
  spec_models = rbindlist(lapply(otu_list, function(x){
    spec_mod = generate_genomic_network(contribution_table[OTU==x, unique(Gene)], keggSource = "KeggTemplate", degree_filter = 40, rxn_table = network_template, normalize = F)[[3]]
    spec_mod[,OTU:=x]
    return(spec_mod)
  }))
  return(spec_models)
}


#' Returns the column from the picrust tables that corresponds to the genomic content of the indicated OTU
#'
#' @import data.table
#' @param otu OTU ID
#' @param picrust_ko_table_directory Directory of PICRUSt genome OTU predictions
#' @param picrust_ko_table_suffix File naming of PICRUSt genome OTU predictions
#' @return Genomic content for a single OTU
#' @examples
#' get_genomic_content_from_picrust_table(otu, picrust_dir, picrust_suffix)
#' @export
get_genomic_content_from_picrust_table = function(otu, picrust_ko_table_directory, picrust_ko_table_suffix){

  # Read in the file for this otu's genomic content
  otu_genomic_content = fread(paste(picrust_ko_table_directory, otu, picrust_ko_table_suffix, sep=""), header=T)
  colnames(otu_genomic_content) = c("OTU", "Gene", "copy_number")
  otu_genomic_content[,OTU:= as.character(OTU)]


  return(otu_genomic_content)
}

#' Returns the melted picrust ko table corresponding to the given set of OTUs
#'
#' @import data.table
#' @param otus List of OTUs to find genomic content for
#' @param picrust_ko_table_directory Directory of PICRUSt genome OTU predictions
#' @param picrust_ko_table_suffix File naming of PICRUSt genome OTU predictions
#' @return Genomic content for list of OTUs
#' @examples
#' get_subset_picrust_ko_table(otus, picrust_dir, picrust_suffix)
#' @export
get_subset_picrust_ko_table = function(otus, picrust_ko_table_directory, picrust_ko_table_suffix){

  # Get all columns corresponding to the otus and transpose
  subset_picrust_ko_table = rbindlist(lapply(otus, get_genomic_content_from_picrust_table, picrust_ko_table_directory = picrust_ko_table_directory, picrust_ko_table_suffix = picrust_ko_table_suffix), use.names=T)

  return(subset_picrust_ko_table)
}

#' Uses the provided OTU table to generate a contribution table based on the PICRUSt 16S normalization and genomic content tables
#'
#' @import data.table
#' @param otu_table OTU/seq table
#' @param picrust_norm_file File path to PICRUSt 16S normalization reference data
#' @param picrust_ko_table_directory Directory of PICRUSt genome OTU predictions
#' @param picrust_ko_table_suffix File naming of PICRUSt genome OTU predictions
#' @return Table of PICRUSt-based contribution abundances for all OTUs
#' @examples
#' generate_contribution_table_using_picrust(otu_table, picrust_norm_file, picrust_dir, picrust_suffix)
#' @export
generate_contribution_table_using_picrust = function(otu_table, picrust_norm_file, picrust_ko_table_directory, picrust_ko_table_suffix){

  #Melt table and convert to relative abundances
  otu_table = melt(otu_table, id.var = "OTU", value.name = "abundance", variable.name = "Sample")
  otu_table[,abundance:=abundance/sum(abundance), by=Sample]
  otu_table[,OTU:=as.character(OTU)]
  #otu_table = data.table(OTU = otu_table[,as.character(OTU)], otu_table[,lapply(.SD, function(x){ return(x/sum(x))}), .SDcols = names(otu_table)[names(otu_table) != "OTU"]])

  # Read the normalization table and standardize column names
  picrust_normalization_table = fread(paste("gunzip -c ", picrust_norm_file, sep=""), header=T)
  colnames(picrust_normalization_table) = c("OTU", "norm_factor")
  picrust_normalization_table[,OTU:= as.character(OTU)]

  if(!all(otu_table[,OTU] %in% picrust_normalization_table[,OTU])) stop("Not all OTUs found in PICRUSt table")


  # Get the subset of the ko table mapping present OTUs to their genomic content
  subset_picrust_ko_table = get_subset_picrust_ko_table(otu_table[,unique(OTU)], picrust_ko_table_directory = picrust_ko_table_directory,
                                                        picrust_ko_table_suffix = picrust_ko_table_suffix)

  # Merge with the table of 16S normalization factors
  contribution_table = merge(otu_table, picrust_normalization_table, by = "OTU", allow.cartesian = TRUE, sort = FALSE)

  # Merge with the PICRUSt genomic content table
  contribution_table = merge(contribution_table, subset_picrust_ko_table, by = "OTU", allow.cartesian = TRUE, sort = FALSE)

  # Make a contribution column by dividing the OTU abundances by the normalization factors and multiplying with the KO copy number
  contribution_table[,contribution := (abundance * copy_number) / norm_factor]

  contribution_table = contribution_table[,list(Sample, OTU, Gene, contribution)]

  # Convert sample, OTU, and function names to character type
  contribution_table[,Sample:= as.character(Sample)]
  contribution_table[,OTU:= as.character(OTU)]
  contribution_table[,Gene:= as.character(Gene)]

  return(contribution_table)
}
