### Taxa mapping functions, MIMOSA2
# June 2018

build_species_networks_w_agora = function(species_dat, database, closest, simThreshold){
  otu_col = names(species_dat)[names(species_dat) %in% c("OTU", "#OTU ID", "Sequence", "ASV")]
  if(length(otu_col) != 1) stop("Species data must have a valid OTU/sequence column name")
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

get_rep_seqs_from_otus = function(otus, database = "Greengenes 13_5 or 13_8", rep_seq_path = "data/rep_seqs/"){
  if(grepl("Greengenes", database)){
    dat_file = paste0(rep_seq_path, "gg_13_5.fasta.gz")
  } else {
    dat_file = paste0(rep_seq_path, "silva_132_99_16S.fna")
  }
  dat_index = data.table(fasta.index(dat_file))
  desired_otus = dat_index[desc %in% otus]
  if(nrow(desired_otus)==0) stop("No matching sequences found, is your database format correct?")
  seq_dat = Biostrings::readDNAStringSet(dat_file)
  seq_dat = seq_dat[which(dat_index[,desc] %in% otus)]
  return(seq_dat)
}

blast_seqs = function(seqs, blast_path = "./data/blastDB/"){
  blastdb = rBLAST::blast(db = blast_path)
  blast_out = list()
  for(j in 1:length(seqs)){
    blast_out[[j]] = rBLAST::predict(blastdb, seqs[j,])
  }
  return(blast_out)
}


load_agora_models = function(agora_path = "data/AGORA/"){

}
