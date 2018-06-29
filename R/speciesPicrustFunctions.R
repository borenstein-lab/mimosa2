### PICRUSt-single species model functions - modified from BURRITO

#Function called by run_pipeline to get species-specific KOs and reactions
build_generic_network = function(species_table, database, picrust_paths, kegg_paths){
  if(database=="Greengenes 13_5 or 13_8"){
    contribution_table = generate_contribution_table_using_picrust(species_table, picrust_norm_file = picrust_paths[1], picrust_ko_table_directory = picrust_paths[2], picrust_ko_table_suffix = picrust_paths[3])
    contribution_table = contribution_table[contribution != 0]
    all_kegg = get_kegg_reaction_info(kegg_paths[2], reaction_info_file = kegg_paths[3], save_out = F)
    network_template = generate_network_template_kegg(kegg_paths[1], all_kegg = all_kegg, write_out = F)
    otu_list = contribution_table[,sort(unique(OTU))]
    spec_models = rbindlist(lapply(otu_list, function(x){
      spec_mod = network_template[KO %in% contribution_table[OTU==x, unique(KO)]]
      spec_mod[,OTU:=x]
      return(spec_mod)
    }))
  } else {
    stop("Only Greengenes currently implemented")
  }
  return(spec_models)
}


# get_genomic_content_from_picrust_table(otu)
#
# Returns the column from the picrust tables that corresponds to the genomic content of the indicated OTU
get_genomic_content_from_picrust_table = function(otu, picrust_ko_table_directory, picrust_ko_table_suffix){

  # Read in the file for this otu's genomic content
  otu_genomic_content = fread(paste(picrust_ko_table_directory, otu, picrust_ko_table_suffix, sep=""), header=T)
  colnames(otu_genomic_content) = c("OTU", "Gene", "copy_number")
  otu_genomic_content[,OTU:= as.character(OTU)]


  return(otu_genomic_content)
}

# get_subset_picrust_ko_table(otus)
#
# Returns the melted picrust ko table corresponding to the given set of OTUs
get_subset_picrust_ko_table = function(otus, picrust_ko_table_directory, picrust_ko_table_suffix){

  # Get all columns corresponding to the otus and transpose
  subset_picrust_ko_table = rbindlist(lapply(otus, get_genomic_content_from_picrust_table, picrust_ko_table_directory = picrust_ko_table_directory, picrust_ko_table_suffix = picrust_ko_table_suffix), use.names=T)

  return(subset_picrust_ko_table)
}

# generate_contribution_table_using_picrust(otu_table)
#
# Uses the provided OTU table to generate a contribution table based on the PICRUSt 16S normalization and genomic content tables
generate_contribution_table_using_picrust = function(otu_table, picrust_norm_file, picrust_ko_table_directory, picrust_ko_table_suffix){

  #Convert table to relative abundances
  otu_table = data.table(OTU = otu_table[,OTU], otu_table[,lapply(.SD, x/sum(x)), .SDcols = names(otu_table)[names(otu_table) != "OTU"]])

  # Read the normalization table and standardize column names
  picrust_normalization_table = fread(paste("zcat ", picrust_norm_file, sep=""), header=T)
  colnames(picrust_normalization_table) = c("OTU", "norm_factor")
  picrust_normalization_table[,OTU:= as.character(OTU)]

  if(!all(otu_table[,OTU] %in% picrust_normalization_table[,OTU])) stop("Not all OTUs found in PICRUSt table")

  # Get the subset of the ko table mapping present OTUs to their genomic content
  subset_picrust_ko_table = get_subset_picrust_ko_table(levels(factor(otu_table[,OTU])))

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
