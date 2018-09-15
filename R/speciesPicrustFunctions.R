### PICRUSt-single species model functions - modified from BURRITO

#' Set column name of OTU table to "OTU", and remove taxa with no abundance
#'
#' @import data.table
#' @param species_table Species/OTU table
#' @return Species table with column name as OTU
#' @examples
#' spec_table_fix(spec_table)
#' @export
spec_table_fix = function(species_table){
  otu_col = names(species_table)[names(species_table) %in% c("OTU", "#OTU ID", "Sequence", "ASV")]
  if(length(otu_col) != 1) stop("Species data must have a valid OTU/sequence column name")
  setnames(species_table, otu_col, "OTU")
  species_table = species_table[rowSums(species_table[,names(species_table) != "OTU", with=F]) > 0]
  return(species_table)
}

#' Set column name to compound and filter metabolites with too many zeros or NAs
#'
#' @import data.table
#' @param met_table Metabolite table
#' @param nzero_filt Number of samples that must have a nonzero/non-NA value to include that metabolite
#' @return Filtered metabolite table
#' @examples
#' met_table_fix(met_table)
#' @export
met_table_fix = function(met_table, nzero_filt = 5){
  met_col_name = names(met_table)[names(met_table) %in% c("compound", "KEGG", "Compound", "metabolite", "Metabolite")]
  if(length(met_col_name) != 1) stop("Ambiguous metabolite ID column name, must be one of Compound/compound/KEGG/Metabolite/metabolite")
  setnames(met_table, met_col_name, "compound")
  #Set NAs to 0
  for(j in names(met_table)){
    set(met_table ,which(is.na(met_table[[j]])),j,0)
  }
  met_table = met_table[apply(met_table[,names(met_table) != "compound", with=F], 1, function(x){
    return(length(x[x != 0]) >= nzero_filt)
  })]
  return(met_table)
}

#' Filter species found in few samples and/or low abundance
#'
#' @import data.table
#' @param species_dat OTU table
#' @param filter_type Either "mean" or "fracNonzero"
#' @param minMeanAbund minimum mean relative abundance for inclusion of a taxon
#' @param minSampFrac minimum fraction of samples with nonzero abundance for inclusion of a taxon
#' @return Filtered species table
#' @examples
#' filter_species_abunds(species_dat)
#' @export
filter_species_abunds = function(species_dat, filter_type = "mean", minMeanAbund = 0.001, minSampFrac = 0.01){
  species2 = melt(species_dat, id.vars = "OTU")
  species2[,relAbund:=value/sum(value), by=variable]
  if(filter_type=="mean"){
    mean_abunds = species2[,mean(relAbund), by=OTU]
    good_species = mean_abunds[V1 >= minMeanAbund, OTU]
  } else if(filter_type=="fracNonzero"){
    abund_counts = species2[,sum(relAbund != 0)/length(relAbund), by=OTU]
    good_species = abund_counts[V1 >= minSampFrac, OTU]
  }
  return(species_dat[OTU %in% good_species])
}

#' Function called by run_pipeline to get species-specific reaction network
#'
#' @import data.table
#' @param contribution_table OTU/seq table
#' @param kegg_paths paths to KEGG network files
#' @return Table of species-specific KEGG reactions
#' @examples
#' build_generic_network(species_table, "Greengenes 13_5 or 13_8", picrust_paths, kegg_paths)
#' @export
build_generic_network = function(contribution_table, kegg_paths){
  if(length(list.files(path = gsub("/reaction.*", "", kegg_paths[3]), pattern = "network_template.txt")) > 0){
    network_template = fread(paste0(gsub("/reaction.*", "", kegg_paths[3]), "/network_template.txt"))
  } else {
    all_kegg = get_kegg_reaction_info(kegg_paths[2], reaction_info_file = kegg_paths[3], save_out = F)
    network_template = generate_network_template_kegg(kegg_paths[1], all_kegg = all_kegg, write_out = F) #We should really speed this thing up
  }
  otu_list = contribution_table[,sort(unique(OTU))]
  spec_models = rbindlist(lapply(otu_list, function(x){
    spec_mod = generate_genomic_network(contribution_table[OTU==x, unique(Gene)], keggSource = "KeggTemplate", degree_filter = 0, rxn_table = network_template, return_mats = F)
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
  otu_table[,abundance:=as.numeric(abundance)]
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
  contribution_table = contribution_table[contribution != 0]
  return(contribution_table)
}

#' Imports pre-generated KEGG network models for each species
#'
#' @import data.table
#' @param species_list list of OTUs (typically Greengenes)
#' @param net_path Path to metabolic network models for each species
#' @param net_suffix File suffix for species-level network models
#' @return Table of KEGG reactions for each species
#' @examples
#' get_kegg_network(species_table)
#' @export
get_kegg_network = function(species_list, net_path = "data/picrustGenomeData/indivModels/", net_suffix = "_rxns.txt"){
  all_net = rbindlist(lapply(species_list, function(x){
    net_file = paste0(net_path, x, net_suffix)
    if(file.exists(net_file)) return(fread(net_file)) else return(NULL) #Silently skip if no network file. This is a little dangerous.
  }))
  if(nrow(all_net)==0) stop("Network files not found for this set of species")
  return(all_net)
}
