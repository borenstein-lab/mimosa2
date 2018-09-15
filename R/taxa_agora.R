### Taxa mapping functions, MIMOSA2
# June 2018

#' Finds close AGORA models and then imports those for each species. Builds PICRUSt-based network for remaining species
#'
#' @import data.table
#' @param species_list List of species from mapping
#' @param usePreprocessed Whether to process .mat files or use pre-computed versions
#' @param agora_path Path to AGORA models
#' @return A list with 2 entries - one a species abundance table in terms of the new species, the second a table of reactions for each species
#' @examples
#' build_species_networks_w_agora(species_table, "Greengenes 13_5 or 13_8", simThreshold = 0.99)
#' @export
build_species_networks_w_agora = function(mod_list, usePreprocessed = T, agora_path = "~/Documents/MIMOSA2shiny/data/AGORA/"){
  #Now load AGORA models
  if(usePreprocessed == F){
    agora_mods = load_agora_models(mod_list, agora_path = agora_path) #This takes a long time
    agora_mats = get_S_mats(agora_mods, mod_list, edge_list = T)
    setnames(agora_mats, "Species", "OTU")
  } else {
    agora_mats = rbindlist(lapply(mod_list, function(x){
      fread(paste0(agora_path, x, "_rxns.txt"))
    }))
  }
  setnames(agora_mats, "Species", "OTU")
  return(agora_mats)
}

#' Convert stoichiometric matrix to edge list of reactions
#'
#' @import data.table
#' @param emm Stoichiometric matrix (rows are metabolites, columns are reactions)
#' @return A data.table where each row is a pair of compounds being exchanged in a reaction
#' @examples
#' emm_to_edge_list(s_mat)
#' @export
emm_to_edge_list = function(emm){
  all_rxn_ids = names(emm)[2:ncol(emm)]
  net_melted = melt(emm, id.var = "Compound")[value != 0]
  net_melted[,Prod:=ifelse(value > 0, Compound,0)]
  net_melted[,Reac:=ifelse(value < 0, Compound,0)]
  edge_list = rbindlist(lapply(all_rxn_ids, function(x){
    rxn_sub = net_melted[variable==x]
    if(nrow(rxn_sub[Prod !=0 ]) > 0 & nrow(rxn_sub[Reac != 0]) > 0){
      edge_list_sub = data.table(expand.grid(rxn_sub[,unique(Reac[Reac !=0])], rxn_sub[,unique(Prod[Prod != 0])]))
      setnames(edge_list_sub, c("Reac", "Prod"))
    } else {
      edge_list_sub = rxn_sub[,list(Reac, Prod)]
      edge_list_sub[Reac==0, Reac:=NA]
      edge_list_sub[Prod==0, Prod:=NA]
    }
    edge_list_sub[,KO:=x]
    reac_info = rxn_sub[Reac != 0, list(Reac, value)]
    setnames(reac_info, "value", "stoichReac")
    prod_info = rxn_sub[Prod != 0, list(Prod, value)]
    setnames(prod_info, "value", "stoichProd")
    if(nrow(rxn_sub[Reac != 0]) > 0) edge_list_sub = merge(edge_list_sub, reac_info, by = "Reac", all.x = T)
    if(nrow(rxn_sub[Prod != 0]) > 0) edge_list_sub = merge(edge_list_sub, prod_info, by = "Prod", all.x = T)
  }), fill = T)
  edge_list[,stoichReac:=-1*stoichReac]
  return(edge_list)
}

#' Convert stoichiometric matrix to edge list of reactions
#'
#' @import data.table
#' @param all_mods List of Cobra-formatted metabolic models
#' @param species_names List of species names for each model
#' @param edge_list Whether to return the metabolic network as a stoichiometric matrix or edge list
#' @return Either a stoichiometric matrix or edge list of reactions
#' @examples
#' get_S_mats(all_mods, my_species, edge_list = T)
#' @export
get_S_mats = function(all_mods, species_names, edge_list = F){
  #get S matrices
  all_S_mat = list()
  for(j in 1:length(all_mods)){
    all_S_mat[[j]] = as.matrix(all_mods[[j]]$S)
    row.names(all_S_mat[[j]]) = unlist(all_mods[[j]]$mets)
    colnames(all_S_mat[[j]]) = unlist(all_mods[[j]]$rxns)
  }
  if(edge_list == T){
    all_S_mat = rbindlist(lapply(1:length(all_S_mat), function(x){
      foo = emm_to_edge_list(data.table(all_S_mat[[x]], Compound = row.names(all_S_mat[[x]])))
      foo[,Species:=species_names[x]]
      return(foo)
    }))
  }
  return(all_S_mat)
}

build_model_components = function(all_mods, species_names, remove_rev = T, missing_rxns = F, agora = F){
  #ID reversible reactions
  reversible_rxns = data.table()
  for(j in 1:length(all_mods)){
    reversible_rxns = rbind(reversible_rxns, data.table(Rxn = unlist(all_mods[[j]]$rxns), Rev = unlist(all_mods[[j]]$rev), Species = species_names[j]))
  }
  reversible_rxns = reversible_rxns[Rxn != "biomass0"]
  reversible_rxns = reversible_rxns[!grepl("EX_", Rxn)]

  #get S matrices
  all_S_mat = list()
  for(j in 1:length(all_mods)){
    all_S_mat[[j]] = as.matrix(all_mods[[j]]$S)
    row.names(all_S_mat[[j]]) = unlist(all_mods[[j]]$mets)
    colnames(all_S_mat[[j]]) = unlist(all_mods[[j]]$rxns)
  }

  all_comps = lapply(all_S_mat, function(x){ return(row.names(x))})
  #   length(Reduce(intersect, all_comps)) #572
  #   length(Reduce(union, all_comps)) #1448

  all_rxns = lapply(all_S_mat, function(x){ return(colnames(x))})
  #   length(Reduce(intersect, all_rxns)) #389
  #   length(Reduce(union, all_rxns)) #1930
  #Well, lots of differences

  all_S_mat = lapply(all_S_mat, function(x){
    foo = data.table(Compound = row.names(x), x)
    foo[,biomass0:=NULL] #don't need this
    foo = foo[,which(!grepl("EX_",names(foo))),with=F]
    return(foo)})

  all_S_mats = all_S_mat[[1]]
  for(j in 2:length(all_S_mat)){
    all_S_mats = merge(all_S_mats, all_S_mat[[j]], by=intersect(names(all_S_mats), names(all_S_mat[[j]])), all=T)
  }
  #Community EMM
  #Get rid of duplicates
  all_S_mats = all_S_mats[,lapply(.SD, function(x){ if(length(x[!is.na(x)]) > 0 ) return(unique(x[!is.na(x)])) else return(0) }), by=Compound]
  setkey(all_S_mats, NULL)

  if(agora == F){
    #Get counts of each Rxn for every species
    all_geneRxn_mats = lapply(1:length(all_mods), function(x){
      mat = data.table(t(as.matrix(all_mods[[x]]$rxnGeneMat)))
      setnames(mat, unlist(all_mods[[x]]$rxns))
      mat[,Gene:=unlist(all_mods[[x]]$genes)]
      mat[,Species:=species_names[x]]
      mat = mat[Gene != "kb"]
      mat[Gene=="Unknown", Gene:=paste0("Unknown",species_names[x])]
      return(mat)
    })
    if(missing_rxns){
      #For every reaction annotated in a species' S mat but without any gene in the gene mat, add unknown gene
      #Mostly transport reactions
      for(k in 1:length(all_mods)){
        missing_genes = names(all_S_mat[[k]])[which(unlist(all_geneRxn_mats[[k]][,lapply(.SD, function(x){ sum(x) }), .SDcols = names(all_S_mat[[k]])[which(names(all_S_mat[[k]]) != "Compound")]])==0)+1]
        for(i in 1:length(missing_genes)){
          set(all_geneRxn_mats[[k]], all_geneRxn_mats[[k]][,grep("Unknown",Gene)[1]], which(names(all_geneRxn_mats[[k]])==missing_genes[i]), 1)
        }
      }
    }
    all_geneRxn_mat = rbindlist(all_geneRxn_mats, fill = T)
    for (j in 1:ncol(all_geneRxn_mat)){
      set(all_geneRxn_mat,which(is.na(all_geneRxn_mat[[j]])),j,0)
    }
    specRxns = all_geneRxn_mat[,lapply(.SD, sum), by=Species, .SDcols=which(!names(all_geneRxn_mat)  %in% c("Gene","Species", "biomass0") & !grepl("EX_", names(all_geneRxn_mat)))] #
    specRxns = specRxns[,c(1, which(colSums(specRxns[,2:ncol(specRxns),with=F])!=0)+1),with=F]
    specRxns = data.table(t(specRxns[,2:ncol(specRxns),with=F]), RxnID = names(specRxns)[names(specRxns) != "Species"])
    setnames(specRxns, c(species_names, "RxnID"))
  } else {
    specRxns = rbindlist(lapply(1:length(all_mods), function(x){
      rxn_dat = data.table(Species = species_names, RxnID = unlist(all_mods[[x]]$rxns), value = 1)
    }))
    specRxns = dcast(specRxns, RxnID~Code, value.var="value", fill = 0)
  }

  #Separate reactions annotated differently in different species
  dups = unique(all_S_mats[duplicated(Compound),Compound])
  if(length(dups) > 0){
    conflict_rxns = list()
    for(i in 1:length(dups)){
      nonzero = all_S_mats[Compound==dups[i]][,names(which(sapply(.SD, uniqueN)!=1))]
      conflict_rxns[[i]] = names(all_S_mats[Compound==dups[i]][,nonzero,with=F])
    }
    all_conflict_rxns = sort(unique(unlist(conflict_rxns)))

    for(j in 1:length(all_conflict_rxns)){
      dup_vals = all_S_mats[Compound %in% dups, unique(get(all_conflict_rxns[j])),by=Compound][duplicated(Compound),Compound]
      vals = all_S_mats[Compound %in% dup_vals, list(Compound,get(all_conflict_rxns[j]))]
      #We need to figure out which compound coefficients go with which from the original S mat
      all_options = list()
      spec_list = list()
      count_opt = 1
      any_match = F
      for(m in 1:length(species_names)){ #check which version each species has
        if(all_conflict_rxns[j] %in% names(all_S_mat[[m]])){
          opt1 = all_S_mat[[m]][Compound %in% dup_vals,list(Compound,get(all_conflict_rxns[j]))][order(Compound)]
          any_match = F
          q = 1
          while(any_match ==F & q <= length(all_options)){
            if(all(opt1 == all_options[[q]])){
              spec_list[[q]] = c(spec_list[[q]], species_names[m])
              any_match = T
            }
            q = q+1
          }
          if(any_match == F){ #if no matching version for this yet
            all_options[[count_opt]] = opt1
            spec_list[[count_opt]] = species_names[m]
            count_opt = count_opt + 1
          }
        }
      }
      for(k in 1:length(all_options)){
        setnames(all_options[[k]], "V2", paste0(all_conflict_rxns[j], "_",k))
        all_S_mats = merge(all_S_mats, all_options[[k]], by="Compound", all = T)
        set(all_S_mats, i=which(is.na(all_S_mats[,get(paste0(all_conflict_rxns[j], "_",k))])), j=which(names(all_S_mats)==paste0(all_conflict_rxns[j], "_",k)), value = all_S_mats[is.na(get(paste0(all_conflict_rxns[j], "_",k))), get(all_conflict_rxns[j])])
        specRxns = rbind(specRxns, specRxns[RxnID==all_conflict_rxns[j]])
        specRxns[nrow(specRxns), RxnID:=paste0(RxnID, "_", k)]
        species_names_bad = species_names[(Species %in% unlist(spec_list)) & !(Species %in% spec_list[[k]]),Code]
        specRxns[nrow(specRxns), (spec_codes_bad):=0]
        reversible_rxns[Rxn==all_conflict_rxns[j] & Species %in% spec_list[[k]], Rxn:=paste0(all_conflict_rxns[j], "_",k)]
        #rename rxns in the original S_mats
        spec_ids = match(spec_list[[k]], spec_codes[,Species])
        for(m in 1:length(spec_ids)){
          setnames(all_S_mat[[spec_ids[m]]], all_conflict_rxns[j], paste0(all_conflict_rxns[j], "_", k))
        }
      }
      all_S_mats[,(all_conflict_rxns[j]):=NULL]
      specRxns = specRxns[RxnID != all_conflict_rxns[j]]
      setkey(all_S_mats, NULL)
    }
    all_S_mats = unique(all_S_mats)
  }

  #Merge specRxns with reversibility info
  specRxns_melt = melt(specRxns, id.var = "RxnID", variable.name = "Code", value.name = "CopyNum")
  specRxns_melt = merge(specRxns_melt, reversible_rxns, by.x = c("RxnID","Code"), by.y = c("Rxn", "Code"), all.x=T)
  specRxns_melt = specRxns_melt[CopyNum != 0]

  if(remove_rev){
    specRxns_melt_rev = copy(specRxns_melt)
    specRxns_rev = copy(specRxns) ##So that we can modify one and not the other later
    specRxns_melt = specRxns_melt[Rev.V1==0]
    rev_rxn_ids = specRxns_melt_rev[Rev.V1==1, unique(RxnID)] #Reactions which are reversible for any subset of species
    for(k in 1:length(rev_rxn_ids)){ #For species that have the reversible version, remove from specRxns
      spec_fix = specRxns_melt_rev[RxnID==rev_rxn_ids[k] & Rev.V1==1, Code]
      set(specRxns, i = specRxns[,which(RxnID==rev_rxn_ids[k])], j=which(names(specRxns) %in% spec_fix), 0)
    }
    specRxns = specRxns[rowSums(specRxns[,which(names(specRxns) != "RxnID"),with=F]) != 0]

    all_S_mats_rev = all_S_mats
    all_S_mats = all_S_mats[,names(all_S_mats) %in% c("Compound", specRxns[,unique(RxnID)]),with=F]
    all_S_mat_rev = all_S_mat
    for(j in 1:length(all_S_mat)){
      all_S_mat[[j]] = all_S_mat[[j]][,names(all_S_mat[[j]]) %in% c("Compound", specRxns_melt[Species==spec_codes[j,Species],RxnID]),with=F]
    }
  }
  return(list(all_S_mat, all_S_mats, specRxns, specRxns_melt, reversible_rxns, all_S_mat_rev, all_S_mats_rev, specRxns_rev, specRxns_melt_rev))
}



#' Returns AGORA IDs that have been previously mapped to a set of closed-reference OTU IDs
#'
#' @import data.table
#' @param otus List of OTU IDs
#' @param database Database/method used for generating populations in the table - currently Greengenes or SILVA
#' @param data_path File path to OTU-AGORA mapping
#' @examples
#' get_rep_seqs_from_otus(otu_list, "Greengenes 13_5 or 13_8", "data/rep_seqs/")
#' @return Biostrings object containing representative sequences for each OTU in order
#' @export
get_agora_from_otus = function(otus, database = database_choices[2], data_path = "data/rep_seqs/agora_mapping/"){
  if(grepl("Greengenes", database)){
    dat_file = paste0(data_path, "gg_agora_vsearch_matches.txt")
  } else {
    dat_file = paste0(data_path, "silva_agora_vsearch_matches.txt") # dat_file = paste0(rep_seq_path, "silva_132_99_16S.fna")
  }
  seq_mapping = fread(dat_file)
  return(seq_mapping[OTU %in% otus])
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

randomString <- function() { #5 random letters and 5 random numbers
  a <- do.call(paste0, replicate(5, sample(LETTERS, 1, TRUE), FALSE))
  paste0(a, sprintf("%05d", sample(9999, 1, TRUE)))
}

#' Returns AGORA species that a list of Greengenes OTUs were mapped to
#'
#' @import data.table
#' @param species GG OTU abundance table
#' @param map_file_path File path to output of gg->agora vsearch mappings
#' @examples
#' gg_to_agora(otu_list)
#' @return List of AGORA species to use for model building
#' @export
gg_to_agora = function(otus, map_file_path = "data/rep_seqs/gg_13_8_99_toAGORA_97_map.txt"){
  map_table = fread(map_file_path, colClasses = "character")
  species_melt = melt(species, id.var = "OTU", variable.name = "Sample")
  species_melt = merge(species_melt, map_table, all.x = T)
  new_species = dcast(species_melt[,sum(value), by=list(AGORA_ID, Sample)], AGORA_ID~Sample, fill = 0)
  setnames(new_species, "AGORA_ID", "OTU")
  return(new_species)
}

#' Assign seq vars to OTUs or AGORA models using vsearch
#'
#' @import devtools
#' @import data.table
#' @import Biostrings
#' @param seqs vector of sequence variants
#' @param repSeqDir File path to reference database
#' @param repSeqFile File name with reference sequences
#' @param method Only "vsearch" currently implemented
#' @param vsearch_path Path to vsearch executable
#' @param file_prefix File name prefix for output
#' @param seqID threshold for vsearch --usearch-global search
#' @param add_agora_names Whether to add AGORA IDs to the table
#' @param otu_tab Whether to return an OTU table, or just the matched sequences
#' @return Table of alignment results (original sequence, hit ID)
#' @examples
#' get_otus_from_seqvar(seqs)
#' @export
get_otus_from_seqvar = function(seqs, repSeqDir = "~/Documents/MIMOSA2shiny/data/blastDB/", repSeqFile = "agora_NCBI_16S.udb", method = "vsearch", vsearch_path = "vsearch",
                                file_prefix = "seqtemp", seqID = 0.99, add_agora_names = T, otu_tab = F){
  file_prefix = paste0(file_prefix, randomString())
  seqList = DNAStringSet(seqs)
  names(seqList) = paste0("seq", 1:length(seqList))
  if(method=="vsearch"){
    writeXStringSet(seqList, filepath = paste0(repSeqDir, file_prefix, ".fasta"))
    command_to_run = paste0(vsearch_path, " --usearch_global ", repSeqDir, file_prefix, ".fasta --db ", repSeqDir, repSeqFile, " --id ", seqID," --strand both --blast6out ", repSeqDir, file_prefix, "vsearch_results.txt")
    if(otu_tab) command_to_run = paste0(command_to_run, " --otutabout ", repSeqDir, file_prefix, "otu_tab.txt")
    if(add_agora_names) command_to_run = paste0(command_to_run, " --maxaccepts 20 --maxrejects 500") #More comprehensive search
    print(command_to_run)
    system(command_to_run)
    # results = readDNAStringSet(paste0(repSeqDir, file_prefix, "vsearch_results.fna"))
    # seq_matches = data.table(seqID = names(results)[seq(1,length(results), by = 2)], databaseID = names(results)[seq(2,length(results), by = 2)])
    # seq_matches[,OrigSeq:=as.vector(results)[seq(1,length(results), by = 2)]]
    results = fread(paste0(repSeqDir, file_prefix, "vsearch_results.txt"), header = F)
    setnames(results, paste0("V", 1:6), c("seqID", "dbID", "matchPerc", "alnlen", "mism", "gapopens"))
    if(results[,length(dbID), by=seqID][,any(V1 != 1)]){
      results[,max_ID:=matchPerc==max(matchPerc), by=seqID]
      results_keep = results[max_ID==T]
      results_keep[,longestAln:=abs(alnlen-max(alnlen)) < 5, by=seqID]
      results_keep = results_keep[longestAln==T]
      #Just keep the first one after this
      results_keep[,count:=order(dbID), by=seqID]
      seq_matches = results_keep[count==1, list(seqID, dbID, matchPerc, alnlen)]
    } else {
      seq_matches = results
    }
    if(add_agora_names){
      seq_data = fread(paste0(repSeqDir, "agora_NCBItax_processed_nodups.txt"))
      seq_matches = merge(seq_matches, seq_data, by.x = "dbID", by.y = "databaseID", all.x = T)
    }
    if(otu_tab){
      otu_table = fread(paste0(repSeqDir, file_prefix, "otu_tab.txt"))
      system(paste0("rm ", repSeqDir, file_prefix, "*"))
      return(otu_table)
    } else {
      system(paste0("rm ", repSeqDir, file_prefix, "*"))
      return(seq_matches)
    }
  }
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

#' #' Aligns a set of sequences against a database using vsearch
#' #'
#' #' @import data.table
#' #' @param seqs Biostrings object of sequences
#' #' @param database_path Path to database
#' #' @return Mapping results
#' #' @examples
#' #' vsearch_map_seqs(seqs, database_path = "path")
#' #' @export
#' vsearch_map_seqs = function(seqs, database_path = "data/blastDB/"){
#'   write(seqs, file = "") #Write temp fasta (Biostrings?)
#'   for(j in 1:length(seqs)){
#'     system(paste0("vsearch --usearch_global ", seqs[j], " --strand=both --top_hits_only --db=", database_path))
#'   }
#'   system("rm temp*") #rm temp fasta
#'   hits_list = fread("vsearch_out")
#'   system("rm vsearch_out*") #rm temp fasta
#'
#'   return(hits_list)
#' }




#' Reads metabolic model in Cobra/Matlab format
#'
#' @import R.matlab
#' @param matFile Path to matlab-formatted AGORA model
#' @examples
#' getModelInfo(otu_file)
#' @export
getModelInfo = function(matFile){
  return(readMat(matFile)[[1]][,,1])
}

#' Reads a set of AGORA models from Matlab files
#'
#' @import data.table
#' @import R.matlab
#' @param mod_list List of species/model IDs to load
#' @param agora_path Path to AGORA files
#' @return List of model objects from files
#' @examples
#' load_agora_models(species_list, "data/AGORA/")
#' @export
load_agora_models = function(mod_list, agora_path = "data/AGORA/"){
  all_mod_files = paste0(mod_list, ".mat")
  all_mods = vector("list", length(mod_list))
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
#' @return vector of KEGG IDs corresponding to agora_ids
#' @examples
#' agora_kegg_mets(agora_ids)
#' @export
agora_kegg_mets = function(agora_ids){
  setkey(kegg_mapping, met)
  return(kegg_mapping[agora_ids, KEGG])
}

#' Finds AGORA IDs for a list of KEGG metabolites (i.e. when adding to AGORA model)
#'
#' @import data.table
#' @param kegg_ids List of KEGG metabolite IDs
#' @return vector of AGORA/Cobra IDs corresponding to agora_ids
#' @examples
#' kegg_agora_mets(kegg_ids)
#' @export
kegg_agora_mets = function(kegg_ids){
  setkey(kegg_mapping, KEGG)
  return(kegg_mapping[kegg_ids, met])
}

#' Determines whether a list of compounds is in KEGG or AGORA format or neither
#'
#' @import data.table
#' @param comp_list List of compound IDs
#' @return Either "KEGG" or "AGORA"
#' @examples
#' get_compound_format(met_ids)
#' @export
get_compound_format = function(comp_list){
  num_match_kegg = kegg_mapping[,sum(comp_list %in% KEGG)]
  num_match_agora = kegg_mapping[,sum(comp_list %in% met)]
  if(num_match_kegg < 2 & num_match_agora < 2) stop("Metabolites do not appear to be in either KEGG or COBRA format")
  if(num_match_kegg > num_match_agora) return("KEGG") else return("AGORA")
}


