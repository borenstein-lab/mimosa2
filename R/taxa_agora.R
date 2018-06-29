### Taxa mapping functions, MIMOSA2
# June 2018

#' Runs BLAST to find close AGORA models and then imports those for each species. Builds PICRUSt-based network for remaining species
#'
#' @import data.table
#' @param species_dat OTU/sequence/species table of abundances
#' @param database Database/method used for generating populations in the table - currently Greengenes, SILVA, or sequence variants
#' @param closest Boolean for whether to pick the closest model or any within a particular sequence similarity threshold. If F simThreshold needs to be specified.
#' @param simThreshold Sequence similarity threshold for linking Greengenes to AGORA
#' @return Table of reactions for each species
#' @examples
#' build_species_networks_w_agora(species_table, "Greengenes 13_5 or 13_8", closest = F, simThreshold = 0.99)
#' @export
build_species_networks_w_agora = function(species_dat, database, closest, simThreshold = NA){
  if(database != "Sequence variants"){
    seqs = get_rep_seqs_from_otus(species_dat[,get(otu_col)], database = database)
  } else {
    seqs = species_dat[,get(otu_col)]
  }
  blast_results = blast_seqs(seqs)
  #We should change blast DB to include same naming as models
  #Now load AGORA models
  agora_mods = load_agora_models()
  #mapping_file =
  if(closest == T){

  }
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
get_rep_seqs_from_otus = function(otus, database = "Greengenes 13_5 or 13_8", rep_seq_path = "data/rep_seqs/"){
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

#' BLASTS a set of sequences against a database
#'
#' @import data.table
#' @param seqs Biostrings object of sequences
#' @param blast_path Path to BLAST directory
#' @return BLAST results
#' @examples
#' blast_seqs(seqs, blast_path = "data/blastDB")
#' @export
blast_seqs = function(seqs, blast_path = "data/blastDB/"){
  blastdb = rBLAST::blast(db = blast_path)
  blast_out = list()
  for(j in 1:length(seqs)){
    blast_out[[j]] = rBLAST::predict(blastdb, seqs[j,])
  }
  return(blast_out)
}

#' Returns the column from the picrust tables that corresponds to the genomic content of the indicated OTU
#'
#' @import data.table
#' @param out OTU/seq table
#' @param picrust_ko_table_directory
#' @param picrust_ko_table_suffix
#' @examples
#' get_genomic_content_from_picrust_table(otu, picrust_dir, picrust_suffix)
#' @export
getModelInfo = function(matFile){
  return(readMat(matFile)[[1]][,,1])
}

#' Reads AGORA models from Matlab files
#'
#' @import data.table
#' @import R.matlab
#' @param agora_path Path to AGORA files
#' @return List of model objects from files
#' @examples
#' load_agora_models("data/AGORA/")
#' @export
load_agora_models = function(agora_path = "data/AGORAmodels/"){
  all_mod_files = list.files(agora_path)
  all_mods = list()
  for(j in 1:length(all_mod_files)){
    all_mods[[j]] = getModelInfo(paste0(agora_path, all_mod_files[j]))
  }
  return(all_mods)
}

#' Returns the column from the picrust tables that corresponds to the genomic content of the indicated OTU
#'
#' @import data.table
#' @import R.matlab
#' @param mat_met_file Matlab file containing mapping of AGORA metabolite IDs to KEGG
#' @return data.table of KEGG/AGORA metabolite mapping
#' @examples
#' met_id_table = process_mat_met_file(mat_met_file)
#' @export
process_mat_met_file = function(mat_met_file = "BiKEGG-master/AllKEGG2BiGGmet.mat"){
  met_table = readMat(mat_met_file)[[1]][,,1]
  kegg_ids = sapply(met_table[[1]], function(x){
    return(x[[1]][[1]][[1]][1])
  })
  bigg_ids = sapply(met_table[[2]], function(x){
    return(x[[1]][[1]][[1]][1])
  })
  kegg_table = data.table(KEGG = kegg_ids, BiGG = bigg_ids)
  kegg_table[,genericBiGG:=gsub("_[a-z]$", "", BiGG)]
  return(kegg_table)
}

#' Finds KEGG IDs for a list of AGORA metabolites
#'
#' @import data.table
#' @param agora_ids List of AGORA metabolite IDs
#' @param kegg_table Dictionary of AGORA-KEGG IDs
#' @return vector of KEGG IDs corresponding to agora_ids
#' @examples
#' get_kegg_met_ids(agora_ids, kegg_table)
#' @export
get_kegg_met_ids = function(agora_ids, kegg_table){
  agora_ids_search = gsub("\\[.*", "", agora_ids)
  return(kegg_table[genericBiGG %in% agora_ids_search, KEGG])
}



