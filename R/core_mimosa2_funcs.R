### MIMOSA linear model w/residual, then get contributions

#' Fits model scaling total CMPs to metabolite concentrations
#'
#' @import data.table
#' @param species_cmps Table of species contribution abundances
#' @param mets_melt Table of metabolite concentrations
#' @return List of 2 data.tables - one with model summary results, one with model residuals
#' @examples
#' fit_cmp_mods(species_cmps, met_data)
#' @export
fit_cmp_mods = function(species_cmps, mets_melt){
  tot_cmps = species_cmps[,sum(CMP), by=list(compound, Sample)]
  tot_cmps = merge(tot_cmps, mets_melt[,list(compound, Sample, value)], by = c("compound", "Sample"))
  all_comps = tot_cmps[,unique(compound)]
  model_dat = data.table(compound = all_comps)
  resid_dat = data.table(expand.grid(compound = all_comps, Sample = tot_cmps[,unique(Sample)]))
  for(x in 1:length(all_comps)){
    scaling_mod = tot_cmps[compound==all_comps[x], lm(value~V1)]
    scaling_coefs = coef(scaling_mod)
    scaling_resids = resid(scaling_mod)
    model_dat[x,Intercept:=scaling_coefs[1]]
    model_dat[x,Slope:=scaling_coefs[2]]
    model_dat[x,Rsq:=summary(scaling_mod)$r.squared]
    if(length(scaling_resids) != nrow(resid_dat[compound==all_comps[x]])) stop("Missing residuals")
    resid_dat[compound==all_comps[x], Resid:=scaling_resids]
  }
  return(list(model_dat, resid_dat))
}

#
#' Returns a melted data table of species and their cmps, along with residual contributions
#'
#' @import data.table
#' @param species_cmps data table of species and their CMPs
#' @param model_dat Table of model results (from fit_cmp_mods)
#' @param resid_dat Table of model residuals (from fit_cmp_mods)
#' @return A melted data table of species and their CMPs, including "Residual" as an additional species
#' @examples
#' add_residuals(species_cmps, mod_results[[1]], mod_results[[2]])
#' @export
add_residuals = function(species_cmps, model_dat, resid_dat){
  species_cmps = species_cmps[compound %in% model_dat[,compound]] #Let go of metabolites not measured
  all_comps = species_cmps[,unique(compound)]
  resid_dat[,Species:="Residual"]
  if("Resid" %in% names(resid_dat)) setnames(resid_dat, "Resid", "newValue")
  for(x in all_comps){
    if(!is.na(model_dat[compound==x,Slope])){
      species_cmps[compound==x, newValue:=CMP*model_dat[compound==x, Slope]]
    } else {
      species_cmps[compound==x, newValue:=NA]
    }
    species_cmps = rbind(species_cmps, resid_dat[compound==x], fill = T)
  }
  species_cmps = species_cmps[!is.na(newValue)]
  return(species_cmps)
}

#' Calculates contributions to variance from a contribution table
#'
#' @import data.table
#' @param species_contribution_table Table of species contributions (with residuals)
#' @param valueVar Column name to use for values
#' @return Data table of variance shares by each species in the original table for each metabolite
#' @examples
#' calculate_var_shares(species_cmps)
#' @export
calculate_var_shares = function(species_contribution_table, valueVar = "newValue"){ #generic, table of values for each speices and sample and compound
  spec_list = species_contribution_table[,unique(Species)]
  spec_table_wide = dcast(species_contribution_table, Sample+compound~Species, value.var = valueVar, fill = 0, fun.aggregate=sum)
  var_shares = rbindlist(lapply(spec_list, function(y){
    all1 = rbindlist(lapply(spec_list, function(x){
      foo = spec_table_wide[,cov(get(x), get(y), use="complete.obs"), by=compound]
      foo[,Species:=x]
      return(foo)
    }))
    all1[,Species2:=y]
  }))
  var_shares = var_shares[,sum(V1),by=list(compound, Species)]
  tot_vals = species_contribution_table[,sum(get(valueVar)), by = list(compound, Sample)]
  true_met_var = tot_vals[,list(var(V1), mean(V1)), by = compound]
  setnames(true_met_var, c("V1", "V2"), c("Var", "Mean"))
  var_shares = merge(var_shares, true_met_var, by="compound")
  var_shares[,VarShare:=V1/Var]
  return(var_shares)
}

#' Plot species contributions for a single metabolite
#'
#' @import ggplot2
#' @param varShares Dataset of contributions
#' @param metabolite Compound to plot
#' @param include_zeros Whether to plot taxa that do not have any contribution
#' @return plot of contributions
#' @examples
#' plot_contributions(varShares)
#' @export
plot_contributions = function(varShares, metabolite, metIDcol = "metID", include_zeros = F){
  plot_dat = varShares[get(metIDcol)==metabolite]
  if(!include_zeros) plot_dat = plot_dat[VarShare != 0]
  ggplot(plot_dat,  aes(y=VarShare, x = Species, fill = Species)) + geom_bar(stat = "identity") + scale_fill_viridis(option = "plasma", discrete = T) + geom_abline(intercept = 0, slope = 0, linetype = 2) +
    theme(strip.background = element_blank(), axis.ticks.y = element_blank(), axis.text.x = element_text(size=8), axis.text.y = element_blank(), legend.title = element_blank(), strip.text = element_blank(), axis.title.y = element_blank(), panel.spacing = unit(0.15, "inches"), plot.margin = margin(0.2, 0.4, 0.3, 0.1, "inches")) +
    ylab("Contribution to variance") + xlab("Taxon") +  coord_flip() + ggtitle(metabolite)#

}

#' Plot summary of metabolite-species contributions
#'
#' @import ggplot2
#' @import RColorBrewer
#' @param varShares Dataset of contributions
#' @param include_zeros Whether to plot taxa that do not have any contribution
#' @return plot of contributions
#' @examples
#' plot_summary_contributions(varShares)
#' @export
plot_summary_contributions = function(varShares, include_zeros = T){
  varShares[,metID:=met_names(as.character(compound))]
  met_order = varShares[Species=="Residual"][order(VarShare, decreasing = F), metID]
  varShares[,metID:=factor(metID, levels = met_order)]
  spec_order = varShares[Species != "Residual",length(VarShare[abs(VarShare) > 0.05]), by=Species][order(V1, decreasing = T), Species]
  spec_order = c(spec_order, "Residual")
  varShares[,Species:=factor(Species, levels = spec_order)]

  if(!include_zeros){
    varShares = varShares[VarShare != 0]
  }
  plot1 = ggplot(varShares, aes(x=metID, y = as.character(Species), fill = VarShare)) + geom_tile() + theme_minimal() +
    theme(axis.text.x = element_text(angle=90, hjust=1, vjust =0.5), axis.line = element_blank()) +
    scale_fill_gradient(low = brewer.pal(9,"Blues")[1], high = brewer.pal(9, "Blues")[9]) + xlab("Metabolite") + ylab("Taxon")
  return(plot1)
}

#' Read in files for a MIMOSA 2 analysis
#'
#' @import data.table
#' @param file_list List of shiny file inputs to load
#' @param configTable Table of configuration parameters
#' @param app Whether configurations are coming from the Shiny app or otherwise
#' @return list of processed abundance data.tables
#' @examples
#' read_mimosa2_files(input_file_list, config_table)
#' @export
read_mimosa2_files = function(file_list, configTable, app = T){
  if(configTable[V1=="database", V2!=get_text("database_choices")[4]]){
    if(app) species = fread(file_list[["file1"]]$datapath) else species = fread(file_list[["file1"]])
    species = spec_table_fix(species)
    #Filter species using default abundance values
    if(configTable[V1=="specNzeroFrac", is.numeric(V2)]){
      print(configTable[V1=="specNzeroFrac"])
      species = filter_species_abunds(species, filter_type = "fracNonzero", configTable[V1=="specNzeroFrac", V2])
    } else { #Use default values
      species = filter_species_abunds(species, filter_type = "fracNonzero")
    }
    if(configTable[V1=="specMinMean", is.numeric(V2)]){
      print(configTable[V1=="specMinMean"])
      species = filter_species_abunds(species, filter_type = "mean", configTable[V1=="specMinMean", V2])
    } else { #Use default values
      species = filter_species_abunds(species, filter_type = "mean")
    }
  } else {
    #Read metagenome file and hold it as species if we are using it for that
    if(app) species = fread(file_list[["metagenome"]]$datapath) else species = fread(file_list[["metagenome"]])
  }
  #Read metabolites
  if(app) mets = fread(file_list[["file2"]]$datapath) else mets = fread(file_list[["file2"]])
  met_nonzero_filt = ifelse(configTable[V1=="metNzeroFilter", is.numeric(V2)], configTable[V1=="metNzeroFilter", V2], 5)
  mets = met_table_fix(mets, met_nonzero_filt)
  shared_samps = intersect(names(species), names(mets))
  if(length(shared_samps) < 2) stop("Sample IDs don't match between species and metabolites")
  spec_colname = ifelse(configTable[V1=="database", V2==get_text("database_choices")[4]], "KO", "OTU")
  species = species[,c(spec_colname, shared_samps), with=F]
  mets = mets[,c("compound", shared_samps), with=F]
  dat_list = list(species, mets)
  names(dat_list) = c("species", "mets")
  for(extraFile in c("metagenome","netAddFile")){
    if(!(extraFile=="metagenome" & configTable[V1=="database", V2==get_text("database_choices")[4]])){ #SKip metagenome if already read in
      if(extraFile %in% names(file_list)){
        if(!is.null(file_list[[extraFile]])){
          if(file_list[[extraFile]] != F){
            if(app) dat = fread(file_list[[extraFile]]$datapath) else dat = fread(file_list[[extraFile]])
            dat_list[[length(dat_list)+1]] = dat
            names(dat_list[[length(dat_list)]]) = extraFile
          }
        }
      }
    } else {
      #save species as metagenome also if that's what's happening
      dat_list[[extraFile]] = species
    }
  }
  return(dat_list)
}


#' Build species-specific metabolic model for MIMOSA2 analysis
#'
#' @import data.table
#' @param species Table of species/taxon abundances
#' @param config_table Data.table of input files and settings for MIMOSA
#' @param netAdd Table of netowrk information to add to the model
#' @return Data.table of network model of genes and reactions for each species/taxon
#' @examples
#' build_metabolic_model(config_table)
#' @export
build_metabolic_model = function(species, config_table, netAdd = NULL){
  ### Get species to use for network if starting from seq vars
  if(config_table[V1=="database", V2==get_text("database_choices")[1]]){ ### Sequence variatns
    seq_list = species[,OTU]
    if(any(grepl("[0-9]+", seq_list)|grepl("[B|D-F|H-S|U-Z|b|d-f|h-s|u-z]+", seq_list))) stop("Feature IDs have non-nucleotide characters, but the sequence variant input option was selected. If the rows of your table are OTU IDs, select the option for their database source on the input page.")
    if(config_table[V1=="genomeChoices", V2==get_text("source_choices")[1]]) { ## Greengenes
      seq_results = map_seqvar(seq_list, repSeqDir = paste0(config_table[V1=="data_prefix", V2], "rep_seqs/"), repSeqFile = "gg_13_8_99_db.udb", add_agora_names = F, seqID = 0.99) #Run vsearch to get gg OTUs
      species[,seqID:=paste0("seq", 1:nrow(species))]
      samps = names(species)[!names(species) %in% c("OTU", "seqID")]
      new_species = merge(species, seq_results, by = "seqID", all.x=T)
      new_species = new_species[,lapply(.SD, sum), by=dbID, .SDcols = samps]
      setnames(new_species, "dbID", "OTU")
      new_species[is.na(OTU), OTU:=0] #Unassigned
      species = new_species
      mod_list = species[,OTU]
    } else if(config_table[V1=="genomeChoices", V2==get_text("source_choices")[2]]){ ## AGORA
      seq_results = map_seqvar(seq_list, repSeqDir = paste0(config_table[V1=="data_prefix", V2], "blastDB/"), repSeqFile = "agora_NCBI_16S.udb", method = "vsearch", file_prefix = "seqtemp", seqID = config_table[V1=="simThreshold", as.numeric(V2)], add_agora_names = T)
      species[,seqID:=paste0("seq", 1:nrow(species))]
      samps = names(species)[!names(species) %in% c("OTU", "seqID")]
      new_species = merge(species, seq_results, by = "seqID", all.x=T)
      new_species = new_species[,lapply(.SD, sum), by=AGORA_ID, .SDcols = samps]
      new_species[is.na(AGORA_ID), AGORA_ID:="Other"]
      setnames(new_species, "AGORA_ID", "OTU")
      mod_list = seq_results[,unique(AGORA_ID)]
      species = new_species

      # if(database != get_text("database_choices")[1]){
      #   seq_results = get_agora_from_otus(species_dat[,OTU], database = database)
      # } else {
      # }
      ### Convert species abundances to AGORA species IDs

    }  else stop("Model source option not found")
  } else if(config_table[V1=="database", V2==get_text("database_choices")[2]]){ ## GG OTUs
      #Nothing to do
    if(config_table[V1=="genomeChoices", V2==get_text("source_choices")[1]]){ #KEGG
      mod_list = species[,OTU]
    } else if(config_table[V1=="genomeChoices", V2==get_text("source_choices")[2]]){ ## AGORA
      species = otus_to_agora(species, gg_file_path = paste0(config_table[V1=="data_prefix", V2], "rep_seqs/gg_13_8_99_toAGORA_97_map.txt"))
      mod_list = species[!is.na(OTU),OTU]
    } else stop('Model option not implemented')
  } else if(config_table[V1=="database", V2==get_text("database_choices")[3]]){ # SILVA
    if(config_table[V1=="genomeChoices", V2==get_text("source_choices")[2]]){ # AGORA
      species = otus_to_agora(species, "SILVA", silva_file_path = paste0(config_table[V1=="data_prefix", V2], "rep_seqs/silva_132_99_toAGORA_97_map.txt"))
      mod_list = species[!is.na(OTU), OTU]
    } else stop("This combination of taxa format and reaction source is not implemented. Please choose a different option.")
  }
  ### Now build network from mod_list
  if(config_table[V1=="genomeChoices", V2==get_text("source_choices")[1]]){ #KEGG
    if(config_table[V1=="database", !V2 %in% get_text("database_choices")[c(1, 2, 4)]]){
      stop("Only Greengenes currently implemented")
    }
    if(config_table[V1=="database", V2==get_text("database_choices")[4]]){
      #Get network from metagenome KOs
      species = species[rowSums(species[,names(species) != "KO", with=F]) != 0]
      network_template = fread(paste0(config_table[V1=="kegg_prefix", V2], "/network_template.txt")) ##generate_network_template_kegg(kegg_paths[1], all_kegg = kegg_paths[2:3], write_out = F)
      network = generate_genomic_network(species[,unique(KO)], keggSource = "KeggTemplate", degree_filter = 0, rxn_table = network_template, return_mats = F)
    } else { ##database_choices 1 or 3
      network = get_kegg_network(mod_list, net_path = paste0(config_table[V1=="data_prefix", V2], "picrustGenomeData/indivModels/"))
    }
  } else if(config_table[V1=="genomeChoices", V2==get_text("source_choices")[2]]){ #AGORA
    network = build_species_networks_w_agora(mod_list, agora_path = paste0(config_table[V1=="data_prefix", V2], "AGORA/"))
  } else stop('Invalid model format specified')
  if(config_table[V1=="database", V2!=get_text("database_choices")[4]]){
    species = species[OTU %in% network[,OTU]]
  }
  if(config_table[V1=="netAdd", length(V2) != 0]){
    if(config_table[V1=="netAdd", V2!=F]){
      netAdd = fread(config_table[V1=="netAdd", V2])
      network = add_to_network(network, netAdd, data_path = config_table[V1=="data_prefix", V2], kegg_path = config_table[V1=="kegg_prefix", V2])
    }
  }
  # if(config_table[V1=="gapfill", V2 != F]){
  #   #Do stuff
  # }
  if(config_table[V1=="database", V2!=get_text("database_choices")[4]]) network = network[OTU %in% species[,OTU]]
  network = filter_currency_metabolites(network, degree_filter = 30)
  return(list(network, species))
}

#' Updated version of getting all single-species CMP scores for every compound and taxon
#'
#' @import data.table
#' @param species_table OTU abundance table (wide format)
#' @param network Species-specific network table, product of build_network functions
#' @param normalize Whether to normalize rows when making the network EMM
#' @return data.table of cmp scores for each taxon and compound
#' @examples
#' get_species_cmp_scores(species_data, network)
#' @export
get_species_cmp_scores = function(species_table, network, normalize = T, relAbund = T){
  network[is.na(stoichReac), stoichReac:=0] #solve NA problem
  network[is.na(stoichProd), stoichProd:=0]
  network[,stoichReac:=stoichReac*normalized_copy_number] #Add in copy num/16S normalization factor
  network[,stoichProd:=stoichProd*normalized_copy_number]
  spec_list = species_table[,unique(OTU)]
  species_table[,OTU:=as.character(OTU)]
  network[,OTU:=as.character(OTU)]
  species_table = melt(species_table, id.var = "OTU", variable.name = "Sample")
  species_table[,Sample:=as.character(Sample)]
  #Convert species to relative abundance if requested
  if(relAbund){
    species_table[,value:=as.double(value)]
    species_table[,value:=value/sum(value)*100, by=Sample]
    species_table[is.nan(value), value:=0] #Just in case of all-0 samples (although this is bad for other reasons)
  }
  if(length(intersect(spec_list, network[,unique(OTU)]))==0) stop("All taxa missing network information, is this the correct network model?")
  if(!all(spec_list %in% network[,unique(OTU)])) warning("Some taxa missing network information")
  network_reacs = network[,list(OTU, KO, Reac, stoichReac)]
  network_prods = network[,list(OTU, KO, Prod, stoichProd)]
  network_reacs[,stoichReac:=-1*stoichReac]
  setnames(network_reacs, c("Reac", "stoichReac"), c("compound", "stoich"))
  setnames(network_prods, c("Prod", "stoichProd"), c("compound", "stoich"))
  if(normalize){
    network_reacs[,stoich:=as.double(stoich)]
    network_prods[,stoich:=as.double(stoich)]
    network_reacs[,stoich:=stoich/abs(sum(stoich)), by=list(OTU, compound)]
    network_prods[,stoich:=stoich/sum(stoich), by=list(OTU, compound)]
  }
  net2 = rbind(network_reacs, network_prods, fill = T)
  net2 = net2[!is.na(compound)] #remove NAs
  spec_cmps = merge(species_table, net2, by = "OTU", allow.cartesian = T)
  spec_cmps[,CMP:=value*stoich]
  spec_cmps = spec_cmps[,sum(CMP), by=list(OTU, Sample, compound)]
  setnames(spec_cmps, c("OTU", "V1"), c("Species", "CMP"))
  all_comps = spec_cmps[,unique(compound)]
  if(length(intersect(all_comps, kegg_mapping[,KEGG])) < 2){ #If compounds are not KEGG IDs
    #Convert AGORA IDs to KEGG IDs
    spec_cmps[,KEGG:=agora_kegg_mets(compound)]
    spec_cmps = spec_cmps[!is.na(KEGG)]
    spec_cmps = spec_cmps[,sum(CMP), by=list(Species, KEGG, Sample)] #Check that this makes sense
    #separate internal/external?
    setnames(spec_cmps, c("KEGG", "V1"), c("compound", "CMP"))
  }

  return(spec_cmps)
}



#' Updated version of getting all sample-level CMP scores from a KO abundance table
#'
#' @import data.table
#' @param ko_table KO abundance table (wide format)
#' @param network Species-specific network table, product of build_network functions
#' @param normalize Whether to normalize rows when making the network EMM
#' @return data.table of cmp scores for each taxon and compound
#' @examples
#' get_cmp_scores_kos(ko_data, network)
#' @export
get_cmp_scores_kos = function(ko_table, network, normalize = T, relAbund = T){
  network[is.na(stoichReac), stoichReac:=0] #solve NA problem
  network[is.na(stoichProd), stoichProd:=0]
  #network[,stoichReac:=stoichReac*normalized_copy_number] #Add in copy num/16S normalization factor
  #network[,stoichProd:=stoichProd*normalized_copy_number]
  ko_table_melt = melt(ko_table, id.var = "KO", variable.name = "Sample")
  ko_table_melt[,Sample:=as.character(Sample)]
  if(relAbund){
    ko_table_melt[,value:=as.double(value)]
    ko_table_melt[,value:=value/sum(value)*100, by=Sample]
  }
  network_reacs = network[,list(KO, Reac, stoichReac)]
  network_prods = network[,list(KO, Prod, stoichProd)]
  network_reacs[,stoichReac:=-1*stoichReac]
  setnames(network_reacs, c("Reac", "stoichReac"), c("compound", "stoich"))
  setnames(network_prods, c("Prod", "stoichProd"), c("compound", "stoich"))
  if(normalize){
    network_reacs[,stoich:=as.double(stoich)]
    network_prods[,stoich:=as.double(stoich)]
    network_reacs[,stoich:=stoich/abs(sum(stoich)), by=compound]
    network_prods[,stoich:=stoich/sum(stoich), by=compound]
  }
  net2 = rbind(network_reacs, network_prods, fill = T)
  net2 = net2[!is.na(compound)]
  spec_cmps = merge(ko_table_melt, net2, by = "KO", allow.cartesian = T)
  spec_cmps[,CMP:=value*stoich]
  spec_cmps = spec_cmps[,sum(CMP), by=list(Sample, compound)]
  setnames(spec_cmps, "V1", "CMP")
  all_comps = spec_cmps[,unique(compound)]
  if(length(intersect(all_comps, kegg_mapping[,KEGG])) < 2){ #If compounds are not KEGG IDs
    #Convert AGORA IDs to KEGG IDs
    spec_cmps[,KEGG:=agora_kegg_mets(compound)]
    spec_cmps = spec_cmps[!is.na(KEGG)]
    spec_cmps = spec_cmps[,sum(CMP), by=list( KEGG, Sample)] #Check that this makes sense
    #separate internal/external?
    setnames(spec_cmps, c("KEGG", "V1"), c("compound", "CMP"))
  }
  spec_cmps[,Species:="TotalMetagenome"]
  return(spec_cmps)
}

#' Add reactions to a network. Will set stoichiometry and copy number to 1 if missing. Format can be either "KO, Rxn, Prod" for reaction IDs with all transformations and correct stoichiometry, or just "KO" but with reaction IDs that are defined in the KEGG network.
#'
#' @import data.table
#' @param network Data.table of taxa, genes and reactions
#' @param addTable Data.table of taxa, genes and/or reactions to add, or generic genes and reactions to be applied to all taxa
#' @param target_format Format of taxa, genes, and/or reactions to add - must be "KEGG" or "Cobra". If NULL, will try to guess
#' @param source_format Format of taxa, genes, and/or reactions to add - must be "KEGG" or "Cobra". If NULL, will try to guess
#' @param kegg_path
#' @param data_path
#' @return Expanded network table
#' @examples
#' add_to_network(network, netAddTable)
#' @export
add_to_network = function(network, addTable, target_format = NULL, source_format = NULL, kegg_path = "data/KEGGfiles/", data_path = "data/"){
  if(is.null(source_format)) source_format = get_compound_format(network[,unique(Reac)])
  if(is.null(target_format) & "Reac" %in% names(addTable) & "Prod" %in% names(addTable)) target_format = get_compound_format(addTable[,c(unique(Reac), unique(Prod))]) else target_format = "KEGG" #Just assume KEGG if only genes provided
  print(source_format)
  print(target_format)
  table_type = names(addTable)
  print(table_type)
  if("Species" %in% table_type){
    setnames(addTable, "Species", "OTU")
    table_type[table_type == "Species"] = "OTU"
  }
  if(all(c("KO", "Reac", "Prod") %in% table_type)){
    if(target_format != source_format){
      if(target_format == "Cobra"){
        #Convert
        addTable[,Reac:=kegg_agora_mets(Reac)]
        addTable[,Prod:=kegg_agora_mets(Prod)]
      } else {
        addTable[,Reac:=agora_kegg_mets(Reac)]
        addTable[,Prod:=agora_kegg_mets(Prod)]
      }
    }
    if(!all(c("stoichReac", "stoichProd") %in% table_type)){
      addTable[,stoichReac:=1]
      addTable[,stoichProd:=1]
    }
    if(!"normalized_copy_number" %in% table_type){
      if("OTU" %in% table_type){
        copyNums = fread(ifelse(target_format=="Cobra", paste0(data_path, "blastDB/agora_NCBItax_processed_nodups.txt"), paste0("gunzip -c ", data_path, "picrustGenomeData/16S_13_5_precalculated.tab.gz")))
        if(target_format == "Cobra"){
          copyNums = unique(copyNums[,list(AGORA_ID, CopyNum)])
          specID = "AGORA_ID"
        } else {
          setnames(copyNums, c("OTU", "CopyNum"))
          specID = "OTU"
        }
        print(copyNums)
        addTable = merge(addTable, copyNums, by.x = "OTU", by.y = specID, all.x = T, all.y=F)
        addTable[is.na(CopyNum), CopyNum:=1]
        addTable[,normalized_copy_number:=1/CopyNum]
        addTable[,CopyNum:=NULL]
      } else {
        addTable[,normalized_copy_number:=1]
      }
    }
  } else if("KO" %in% table_type){
    if(source_format != "KEGG"){
      warning("For non-KEGG based gene/reaction modifications without specified compounds, only reaction removals will be applied")
    } else {
      full_kegg_table = fread(paste0(kegg_path, "network_template.txt"))
      addTable = merge(addTable, full_kegg_table, by = "KO", all.x = F, all.y = F, allow.cartesian = F)
      print(addTable)
    }
  } else {
    stop("Invalid format for reaction addition table")
  }
  if(!"OTU" %in% table_type){  #No species netAdd
    ## Do removals first
    if("remove" %in% table_type){
      print("removing here")
      removeTable = addTable[remove == T]
      addTable = addTable[remove == F|is.na(remove)]
      addTable[,remove:=NULL]
      if(!all(c("Prod", "Reac") %in% table_type)){
        network = network[!KO %in% removeTable[,KO]]
      } else {
        for(j in 1:nrow(removeTable)){
          network = network[!(KO==removeTable[j,KO] & Prod==removeTable[j,Prod] & Reac==removeTable[j,Reac])]
        }
      }
    }
    #Now add in rxns for all species
    if("OTU" %in% names(network)){ #IF it's not a metagenome network
      all_spec = network[,unique(OTU)]
      new_net = data.table()
      for(spec in all_spec){ #Propagate to all species
        new_net1 = copy(addTable)
        new_net1[,OTU:=spec]
        new_net = rbind(new_net, new_net1)
      }
      addTable = new_net
    }
  } else if("remove" %in% table_type){
    #Do species-specific removals
    removeTable = addTable[remove == T]
    addTable = addTable[remove == F|is.na(remove)]
    addTable[,remove:=NULL]
    if(!all(c("Prod", "Reac") %in% table_type)){ #If just KOs
      for(j in 1:nrow(removeTable)){
        network = network[KO != removeTable[j,KO] & OTU != removeTable[j,OTU]]
      }
    } else {
      for(j in 1:nrow(removeTable)){
        network = network[!(KO==removeTable[j,KO] & Prod==removeTable[j,Prod] & Reac==removeTable[j,Reac] & OTU==removeTable[j,OTU])]
      }
    }
  }
  if("OTU" %in% names(network)){
    network[,OTU:=as.character(OTU)]
    addTable[,OTU:=as.character(OTU)]
  }
  network = rbind(network, addTable, fill = T)
  ## Remove duplicates
  setkey(network, NULL)
  network = unique(network)
  return(network)
}

#' Check configuration table formatting
#'
#' @import data.table
#' @param config_table Table of settings for MIMOSA analysis
#' @param data_path Path to MIMOSA2shiny data files
#' @param app Whether this is being called by the MIMOSA web app
#' @return Cleaned-up configuration table
#' @examples
#' check_config_table(table1)
#' @export
check_config_table = function(config_table, data_path = "data/", app = F){
  if(app){
    req_params = c("database", "genomeChoices")
  } else {
    req_params = c("file1", "file2", "database", "genomeChoices")
    if(config_table[V1=="database", V2==get_text("database_choices")[4]]){
      req_params[req_params=="file1"] = "metagenome"
    }
  }
  if(any(!req_params %in% config_table[,V1])) stop("Required parameters missing from configuration file")
  all_params = c(req_params, "metagenome","contribType", "gapfill", "metType", "netAdd", "simThreshold", "kegg_prefix", "data_prefix") #Move to package sysdata?
  config_table[V2=="", V2:=FALSE]
  if(length(all_params[!all_params %in% config_table[,V1]]) > 0){
    config_table = rbind(config_table, data.table(V1 = all_params[!all_params %in% config_table[,V1]], V2 = FALSE))
  }
  return(config_table)
}

#' Run a MIMOSA 2 analysis
#'
#' @import data.table
#' @param config_table Data.table of input files and settings for MIMOSA analysis
#' @return Scaling model and variance contribution results
#' @examples
#' run_mimosa2(species_file, met_file, config_file)
#' @export
run_mimosa2 = function(config_table){
  #process arguments
  config_table = check_config_table(config_table)
  file_list = as.list(config_table[grepl("file", V1, ignore.case = T)|V1=="metagenome", V2])
  names(file_list) = config_table[grepl("file", V1, ignore.case = T)|V1=="metagenome", V1]
  data_inputs = read_mimosa2_files(file_list, config_table, app = F)
  species = data_inputs$species
  mets = data_inputs$mets

  network_results = build_metabolic_model(species, config_table)
  network = network_results[[1]]
  species = network_results[[2]]
  if(!is.null(data_inputs$metagenome) & config_table[V1=="database", V2!=get_text("database_choices")[4]]){
    #If we are doing a comparison of the species network and the metagenome network
    #Metagenome data
    #Implement doing stuff with this later
    metagenome_network = build_metabolic_model(data_inputs$metagenome, config_table)
    # species2 = metagenome_data[[1]]
    # metagenome_network = metagenome_data[[2]]
    #Metagenome data
  }

  if(config_table[V1=="metType", V2 ==get_text("met_type_choices")[2]]){ #Assume it is KEGG unless otherwise specified
      mets = map_to_kegg(mets)
  }
  if(config_table[V1=="database", V2==get_text("database_choices")[4]]){
    indiv_cmps = get_cmp_scores_kos(species, network) #Use KO abundances instead of species abundances to get cmps
  } else {
    indiv_cmps = get_species_cmp_scores(species, network)
  }
  mets_melt = melt(mets, id.var = "compound", variable.name = "Sample")
  cmp_mods = fit_cmp_mods(indiv_cmps, mets_melt)
  indiv_cmps = add_residuals(indiv_cmps, cmp_mods[[1]], cmp_mods[[2]])
  var_shares = calculate_var_shares(indiv_cmps)
  if(!is.null(data_inputs$metagenome) & config_table[V1=="database", V2!=get_text("database_choices")[4]]){
    indiv_cmps2 = get_cmp_scores_kos(species, metagenome_network)
    cmp_mods2 = fit_cmp_mods(indiv_cmps2, mets_melt)
    indiv_cmps2 = add_residuals(indiv_cmps2, cmp_mods2[[1]], cmp_mods2[[2]])
    var_shares_metagenome = calculate_var_shares(indiv_cmps2)
    return(list(varShares = var_shares, modelData = cmp_mods[[1]], varSharesMetagenome = var_shares_metagenome, ModelDataMetagenome = cmp_mods2))
  } else {
    return(list(varShares = var_shares, modelData = cmp_mods[[1]]))
  }
  #Order dataset for plotting
  #met_order = var_shares[Species=="Residual"][order(VarShare, increasing = T), metID]
  #var_shares[,metID:=factor(metID, levels = metID)]

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
#' map_seqvar(seqs)
#' @export
map_seqvar = function(seqs, repSeqDir = "data/blastDB/", repSeqFile = "agora_NCBI_16S.udb", method = "vsearch", vsearch_path = "vsearch",
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

#' Convert metabolite name table to KEGG metabolite table
#'
#' @import MetaboAnalystR
#' @import data.table
#' @param met_table Table of metabolite abundances
#' @return A new table of metabolite abundances using KEGG compound IDs
#' @examples
#' new_mets = map_to_kegg(mets)
#' @export
map_to_kegg = function(met_table){
  mSet = InitDataObjects("NA", "utils", FALSE)
  cmpds = met_table[,compound]
  mSet = Setup.MapData(mSet, cmpds)
  mSet = CrossReferencing(mSet, "name", kegg = T)
  mSet = CreateMappingResultTable(mSet)
  match_table = data.table(mSet$dataSet$map.table)
  num_nas = nrow(match_table[is.na(KEGG)|KEGG=="NA"])
  if(num_nas > 0) warning(paste0(num_nas, " metabolites were not matched to KEGG IDs and will be ignored"))
  met_table = merge(met_table, match_table[,list(Query, KEGG)], by.x = "compound", by.y = "Query", all.x = T)[!is.na(KEGG) & KEGG != "NA"]
  met_table = met_table[,lapply(.SD, sum), by=KEGG, .SDcols = names(met_table)[!names(met_table) %in% c("compound", "KEGG")]]
  setnames(met_table, "KEGG", "compound")
  return(met_table)
}
