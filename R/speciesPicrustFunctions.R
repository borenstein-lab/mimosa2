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
  if(length(otu_col) != 1) stop("Species data must have a valid OTU/sequence column name, one of: OTU/#OTU ID/Sequence/ASV")
  setnames(species_table, otu_col, "OTU")
  species_table = species_table[rowSums(species_table[,names(species_table) != "OTU", with=F], na.rm = T) > 0]
  if(any(duplicated(species_table[,OTU]))){ #Deal with duplicated features
    dup_otus = species_table[duplicated(OTU), OTU]
    for(k in 1:length(dup_otus)){
      dup_rows = species_table[,which(OTU==dup_otus[k])]
      species_table = species_table[-dup_rows[2:length(dup_rows)]] #Remove rows
      warning(paste0(dup_otus[k], " is specified in ", length(dup_rows), " rows of the microbiome table. Duplicate features are not allowed; all rows with this feature after the first one will be removed"))
    }
  }
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
  if(length(met_col_name) != 1) stop("Ambiguous metabolite ID column name, must be one of: Compound/compound/KEGG/Metabolite/metabolite")
  setnames(met_table, met_col_name, "compound")
  #Set NAs to 0
  for(j in names(met_table)){
    set(met_table ,which(is.na(met_table[[j]])),j,0)
  }
  met_table = met_table[apply(met_table[,names(met_table) != "compound", with=F], 1, function(x){
    return(length(x[x != 0]) >= nzero_filt)
  })]
  if(any(duplicated(met_table[,compound]))){ #Deal with duplicated features
    dup_mets = met_table[duplicated(compound), compound]
    for(k in 1:length(dup_mets)){
      dup_rows = met_table[,which(compound==dup_mets[k])]
      met_table = met_table[-dup_rows[2:length(dup_rows)]] #Remove rows
      warning(paste0(dup_mets[k], " is specified in ", length(dup_rows), " rows of the metabolite table. Duplicate features are not allowed; all rows with this feature after the first one will be removed"))
    }
  }
  
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
filter_species_abunds = function(species_dat, filter_type = "mean", minMeanAbund = 0, minSampFrac = 0){ #0.001, 0.01 previously
  species_dat[,OTU:=as.character(OTU)]
  species2 = melt(species_dat, id.vars = "OTU")
  species2[,relAbund:=value/sum(value, na.rm = T), by=variable]
  if(filter_type=="mean"){
    mean_abunds = species2[,mean(relAbund, na.rm = T), by=OTU]
    good_species = mean_abunds[V1 >= minMeanAbund, OTU]
  } else if(filter_type=="fracNonzero"){
    abund_counts = species2[,sum(relAbund != 0, na.rm = T)/length(relAbund), by=OTU]
    good_species = abund_counts[V1 >= minSampFrac, OTU]
  }
  return(species_dat[OTU %in% good_species])
}

#' Convert stoichiometric matrix to edge list of reactions
#'
#' @import data.table
#' @param emm Stoichiometric matrix (rows are metabolites, columns are reactions)
#' @return A data.table where each row is a pair of compounds being exchanged in a reaction
#' @examples
#' emm_to_edge_list(s_mat)
#' @export
#' 
emm_to_edge_list = function(emm){
  net_melted = melt(emm, id.var = "Compound")[value != 0]
  net_melted[,Prod:=ifelse(value > 0, Compound,0)]
  net_melted[,Reac:=ifelse(value < 0, Compound,0)]
  all_rxn_ids = net_melted[,unique(as.character(variable))]
  edge_list = data.table()
  
  for(k in 1:length(all_rxn_ids)){
    rxn_sub = net_melted[variable==all_rxn_ids[k]]
    if(nrow(rxn_sub[Prod !=0 ]) > 0 & nrow(rxn_sub[Reac != 0]) > 0){
      edge_list_sub = data.table(expand.grid(rxn_sub[,unique(Reac[Reac !=0])], rxn_sub[,unique(Prod[Prod != 0])]))
      setnames(edge_list_sub, c("Reac", "Prod"))
    } else {
      edge_list_sub = rxn_sub[,list(Reac, Prod)]
      edge_list_sub[Reac==0, Reac:=NA]
      edge_list_sub[Prod==0, Prod:=NA]
    }
    edge_list_sub[,KO:=all_rxn_ids[k]]
    edge_list_sub[,stoichReac:=sapply(Reac, function(x){ return(rxn_sub[Compound==x,abs(value)])})]
    edge_list_sub[,stoichProd:=sapply(Prod, function(x){ return(rxn_sub[Compound==x,value])})]
    edge_list = rbind(edge_list, edge_list_sub)
  }
  return(edge_list)
}

# emm_to_edge_list = function(emm){
#   all_rxn_ids = names(emm)[2:ncol(emm)]
#   net_melted = melt(emm, id.var = "Compound")[value != 0]
#   net_melted[,Prod:=ifelse(value > 0, Compound,0)]
#   net_melted[,Reac:=ifelse(value < 0, Compound,0)]
#   edge_list = rbindlist(lapply(all_rxn_ids, function(x){
#     rxn_sub = net_melted[variable==x]
#     if(nrow(rxn_sub[Prod !=0 ]) > 0 & nrow(rxn_sub[Reac != 0]) > 0){
#       edge_list_sub = data.table(expand.grid(rxn_sub[,unique(Reac[Reac !=0])], rxn_sub[,unique(Prod[Prod != 0])]))
#       setnames(edge_list_sub, c("Reac", "Prod"))
#     } else {
#       edge_list_sub = rxn_sub[,list(Reac, Prod)]
#       edge_list_sub[Reac==0, Reac:=NA]
#       edge_list_sub[Prod==0, Prod:=NA]
#     }
#     edge_list_sub[,KO:=x]
#     reac_info = rxn_sub[Reac != 0, list(Reac, value)]
#     setnames(reac_info, "value", "stoichReac")
#     prod_info = rxn_sub[Prod != 0, list(Prod, value)]
#     setnames(prod_info, "value", "stoichProd")
#     if(nrow(rxn_sub[Reac != 0]) > 0) edge_list_sub = merge(edge_list_sub, reac_info, by = "Reac", all.x = T)
#     if(nrow(rxn_sub[Prod != 0]) > 0) edge_list_sub = merge(edge_list_sub, prod_info, by = "Prod", all.x = T)
#   }), fill = T)
#   edge_list[,stoichReac:=-1*stoichReac]
#   return(edge_list)
# }

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
  file_id = paste(picrust_ko_table_directory, otu, picrust_ko_table_suffix, sep="")
  if(file.exists(file_id)){
    otu_genomic_content = fread(file_id, header=T)
    colnames(otu_genomic_content) = c("OTU", "Gene", "copy_number")
    otu_genomic_content[,OTU:= as.character(OTU)]
    return(otu_genomic_content)
  } else {
    warning(paste0("OTU file ", otu, " not found"))
    return(NULL)
  }
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
#' @param copynum_column Whether to include a copy number column in the contribution table
#' @param otu_rel_abund Whether to convert OTU table to relative abundances first
#' @return Table of PICRUSt-based contribution abundances for all OTUs
#' @examples
#' generate_contribution_table_using_picrust(otu_table, picrust_norm_file, picrust_dir, picrust_suffix)
#' @export
generate_contribution_table_using_picrust = function(otu_table, picrust_norm_file, picrust_ko_table_directory, picrust_ko_table_suffix, copynum_column = F, otu_rel_abund = T){

  #Melt table and convert to relative abundances
  otu_table = melt(otu_table, id.var = "OTU", value.name = "abundance", variable.name = "Sample")
  otu_table[,abundance:=as.numeric(abundance)]
  if(otu_rel_abund){
    otu_table[,abundance:=abundance/sum(abundance), by=Sample]
  }
  otu_table[,OTU:=as.character(OTU)]
  #otu_table = data.table(OTU = otu_table[,as.character(OTU)], otu_table[,lapply(.SD, function(x){ return(x/sum(x))}), .SDcols = names(otu_table)[names(otu_table) != "OTU"]])

  # Read the normalization table and standardize column names
  picrust_normalization_table = fread(picrust_norm_file, header=T) #cmd = paste("gunzip -c ", picrust_norm_file, sep="")
  colnames(picrust_normalization_table) = c("OTU", "norm_factor")
  picrust_normalization_table[,OTU:= as.character(OTU)]

  if(!all(otu_table[,OTU] %in% picrust_normalization_table[,OTU])) warning("Not all OTUs found in PICRUSt table")

  # Get the subset of the ko table mapping present OTUs to their genomic content
  subset_picrust_ko_table = get_subset_picrust_ko_table(otu_table[OTU %in% picrust_normalization_table[,OTU],unique(OTU)], picrust_ko_table_directory = picrust_ko_table_directory,
                                                        picrust_ko_table_suffix = picrust_ko_table_suffix)

  # Merge with the table of 16S normalization factors
  contribution_table = merge(otu_table, picrust_normalization_table, by = "OTU", allow.cartesian = TRUE, sort = FALSE)

  # Merge with the PICRUSt genomic content table
  contribution_table = merge(contribution_table, subset_picrust_ko_table, by = "OTU", allow.cartesian = TRUE, sort = FALSE)

  # Make a contribution column by dividing the OTU abundances by the normalization factors and multiplying with the KO copy number
  contribution_table[,contribution := (abundance * copy_number) / norm_factor]

  if(copynum_column){
    contribution_table = contribution_table[,list(Sample, OTU, Gene, contribution, copy_number/norm_factor)]
    setnames(contribution_table, "V5", "copy_number")
  } else {
    contribution_table = contribution_table[,list(Sample, OTU, Gene, contribution)]
  } 

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

#' Returns representative 16S sequences for a set of closed-reference OTU IDs
#'
#' @import data.table
#' @import Biostrings
#' @param otus List of OTU IDs
#' @param database Database/method used for generating populations in the table - currently Greengenes or SILVA
#' @param rep_seq_path File path to OTU representative sequences
#' @examples
#' get_rep_seqs_from_otus(otu_list, "Greengenes 13_5 or 13_8", "data/rep_seqs/")
#' @return Biostrings object containing representative sequences for each OTU in order
#' @export
get_rep_seqs_from_otus = function(otus, database = database_choices[2], rep_seq_path = "data/rep_seqs/"){
  if(grepl("Greengenes", database)){
    dat_file = paste0(rep_seq_path, "gg_13_5.fasta.gz")
  } else {
    # dat_file = paste0(rep_seq_path, "silva_132_99_16S.fna")
  }
  dat_index = data.table(fasta.index(dat_file))
  desired_otus = dat_index[desc %in% otus]
  if(nrow(desired_otus)==0) stop("No matching sequences found, is your database format correct?")
  seq_dat = Biostrings::readDNAStringSet(dat_file)
  seq_dat = seq_dat[which(dat_index[,desc] %in% otus)]
  return(seq_dat)
}



#' Generate preprocessed reference databases for MIMOSA2
#'
#' @param database One of "KEGG", "picrustGG", "AGORA", or "embl_gems"
#' @param model_table_file Optionally, a table listing models from a version of AGORA or embl_gems other than what is included in package data
#' @param picrust_ko_path Path to unzipped precalculated PICRUSt 1 files
#' @param kegg_paths Vector of file paths to the 3 KEGG files required for MIMOSA network template construction (1) reaction_mapformula.lst, 2) reaction_ko.list, and 3) reaction)
#' @param dat_path For the AGORA or embl_gems options, path to directory of .mat model files
#' @param out_path Path to write processed reaction network outputs for each taxon
#'
#' @return No return, writes processed reaction network to file
#' @export
#'
#' @examples
#' generate_preprocessed_networks("AGORA", dat_path = "data/AGORA/", )
generate_preprocessed_networks = function(database, model_table_file = NULL, picrust_ko_path = "data/picrustGenomeData/", 
                                                  kegg_paths = c("data/KEGGfiles/reaction_mapformula.lst", "data/KEGGfiles/reaction_ko.list", "data/KEGGfiles/reaction"),
                                                  dat_path = paste0("data/", database, "/"), out_path = paste0("data/", database, "/RxnNetworks/")
                                                  ){
  if(database %in% c("KEGG", "picrustGG")){ #anything kegg related
    #Get KEGG network template (same as MIMOSA1)
    if(length(list.files(path = out_path, pattern = "network_template.txt")) > 0){
      network_template = fread(paste0(out_path, "/network_template.txt"))
    } else {
      all_kegg = get_kegg_reaction_info(kegg_paths[2], reaction_info_file = kegg_paths[3], save_out = F)
      network_template = generate_network_template_kegg(kegg_paths[1], all_kegg = all_kegg, write_out = F) 
    }
    if(database == "picrustGG"){
      ## Get list of all OTUs
      picrust_ko_table_directory = picrust_ko_path
      picrust_ko_table_suffix = "_genomic_content.tab"
      all_otus = gsub(picrust_ko_table_suffix, "", list.files(picrust_ko_table_directory))
      
      picrust_norm_file = paste0(picrust_ko_path, "16S_13_5_precalculated.tab")
      #Get normalization data
      picrust_normalization_table = fread(picrust_norm_file, header = T)#fread(paste("gunzip -c ", picrust_norm_file, sep=""), header=T)
      colnames(picrust_normalization_table) = c("OTU", "norm_factor")
      picrust_normalization_table[,OTU:= as.character(OTU)]
    
      for(x in all_otus){
        genomic_content = get_genomic_content_from_picrust_table(x, picrust_ko_table_directory, picrust_ko_table_suffix)
        spec_mod = generate_genomic_network(genomic_content[,unique(Gene)], keggSource = "KeggTemplate", degree_filter = 0, rxn_table = network_template, return_mats = F) #Not going to filter anymore
        #Also get copy number info
        spec_mod = merge(spec_mod, genomic_content, by.x = "KO", by.y = "Gene", all.x=T)
        spec_mod = merge(spec_mod, picrust_normalization_table, by = "OTU", all.x=T, all.y=F)
        spec_mod[,normalized_copy_number:=copy_number/norm_factor]
        spec_mod = spec_mod[,list(OTU, KO, Reac, Prod, stoichReac, stoichProd, normalized_copy_number)]
        write.table(spec_mod, file = paste0(out_path, x, "_rxns.txt"), quote=F, row.names=F, sep = "\t")
      }
    } else {
      write.table(network_template, file = paste0(out_path, "network_template.txt"), quote = F, row.names = F, sep = "\t")
    }
  }else if(grepl("embl", database)){ #Embl gems
    genome_info = get_text("embl_gems_model_data")
    if(is.null(genome_info)){ stop("List of models is required") }
    otu_list = list.files(path = dat_path, pattern = ".mat$")
    otu_list = gsub(".mat$", "", otu_list)
    bad_genomes = genome_info[duplicated(ModelID), ModelID] #One double genome
    genome_info = genome_info[!(ModelID %in% bad_genomes & CopyNum16S==0)]
    for(spec in otu_list){
      mod1 = load_agora_models(spec, agora_path = dat_path)
      mod1 = get_S_mats(mod1, spec, edge_list = T)
      ##Read back in and add copy number info????
      print(spec)
      mod1[,copy_number:=1]
      mod1[,Reac:=gsub("[C_c]", "[c]", Reac, fixed = T)]
      mod1[,Prod:=gsub("[C_c]", "[c]", Prod, fixed = T)]
      mod1[,Reac:=gsub("[C_e]", "[e]", Reac, fixed = T)]
      mod1[,Prod:=gsub("[C_e]", "[e]", Prod, fixed = T)]
      
      mod1 = merge(mod1, genome_info[,list(ModelID, CopyNum16S)], all.x = T, by.x = "Species", by.y = "ModelID")
      mod1[,normalized_copy_number:=ifelse(CopyNum16S==0|is.na(CopyNum16S), 1, copy_number/CopyNum16S)]
      mod1 = mod1[,list(Species, KO, Reac, Prod, stoichReac, stoichProd, normalized_copy_number, LB, UB, Rev)]
      write.table(mod1, file = paste0(out_path, spec, "_rxns.txt"), quote=F, row.names=F, sep = "\t")
    }
    
  } else { #AGORA - type
    genome_info = get_text("agora_model_data")
    if(is.null(genome_info)){ stop("List of models is required")}
    otu_list = list.files(path = dat_path, pattern = ".mat$")
    otu_list = gsub(".mat$", "", otu_list)
    print(otu_list)
    genome_info = fread("data/blastDB/AGORA_full_genome_info.txt")
    for(spec in otu_list){
      ##Read back in and add copy number info
      print(spec)
      mod1 = load_agora_models(spec, agora_path = dat_path)
      mod1 = get_S_mats(mod1, spec, edge_list = T)
      mod1[,copy_number:=1]
      mod1 = merge(mod1, genome_info[,list(ModelAGORA, CopyNum)], all.x = T, by.x = "Species", by.y = "ModelAGORA")
      mod1[,normalized_copy_number:=ifelse(CopyNum==0|is.na(CopyNum), 1, copy_number/CopyNum)]
      mod1 = mod1[,list(Species, KO, Reac, Prod, stoichReac, stoichProd, normalized_copy_number, LB, UB, Rev)]
      write.table(mod1, file = paste0(dat_path, spec, "_rxns.txt"), quote=F, row.names=F, sep = "\t")
    }
  }
}

#' Uses the provided OTU table to generate a contribution table based on the PICRUSt 16S normalization and genomic content tables
#'
#' @import biomartr
#' @param seq_db Microbiome data source type: must be one of "Sequence variants (ASVs)", "Greengenes 13_5 or 13_8 OTUs", or "SILVA 132 OTUs"  
#' @param target_db Metabolic network source: must be one of "AGORA genomes and models" or "RefSeq/EMBL_GEMs genomes and models" (for KEGG models, see documentation)
#' @param save_to File path to save the downloaded data. Default is a folder named "data/" in the current working directory
#' @return None
#' @examples
#' download_reference_data("Sequence variants (ASVs)", "AGORA genomes and models")
#' @export
download_reference_data = function(seq_db = get_text("database_choices")[1], target_db = get_text("source_choices")[2], save_to = "data/"){
  if(grepl("AGORA", target_db, ignore.case = T)){
    target1 = "AGORA"
  } else if(grepl("EMBL", target_db, ignore.case = T)){
    target1 = "EMBL"
  }
  if(grepl("ASV", seq_db, ignore.case = T)){
    source1 = "ASV"
  } else if(grepl("OTU", seq_db, ignore.case = T)){
    source1 = "OTU"
  }
  if(grepl("KEGG", target_db)){
    stop("KEGG resources cannot be downloaded")
  }
  file_id = paste0(source1, "_", target1, ".tar.gz")
  download1 = download.file(paste0("http://elbo-spice.gs.washington.edu/shiny/MIMOSA2shiny/refData/", file_id), destfile = file_id)
  if(download1 != 0) stop("Download failed, please download manually") else {
    cat("Data downloaded to ", getwd())
  }
  untar(file_id)
  file.rename(from = paste0(source1, "_", target1, "/data/"), to = save_to)
  cat("Set up the data directory, ready for mimosa2 analysis\n")
}


#' Uses the provided OTU table to generate a contribution table based on the PICRUSt 16S normalization and genomic content tables
#'
#' @import biomartr
#' @param db Either AGORA or embl_gems
#' @param picrust_norm_file File path to PICRUSt 16S normalization reference data
#' @param picrust_ko_table_directory Directory of PICRUSt genome OTU predictions
#' @param picrust_ko_table_suffix File naming of PICRUSt genome OTU predictions
#' @param copynum_column Whether to include a copy number column in the contribution table
#' @param otu_rel_abund Whether to convert OTU table to relative abundances first
#' @return Table of PICRUSt-based contribution abundances for all OTUs
#' @examples
#' generate_contribution_table_using_picrust(otu_table, picrust_norm_file, picrust_dir, picrust_suffix)
#' @export
download_ribosomal_ref_seqs = function(db = "AGORA", out_path){
  
}
