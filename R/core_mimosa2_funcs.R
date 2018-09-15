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
#' @param gg_path Path to Picrust ref files
#' @param netAdd Table of netowrk information to add to the model
#' @return Data.table of network model of genes and reactions for each species/taxon
#' @examples
#' build_metabolic_model(config_table)
#' @export
build_metabolic_model = function(species, config_table, gg_path =  "data/picrustGenomeData/indivModels/", netAdd = NULL){
  ### Get species to use for network if starting from seq vars
  if(config_table[V1=="database", V2==get_text("database_choices")[1]]){ ### Sequence variatns
    seq_list = species[,OTU]
    if(any(grepl("[0-9]+", seq_list)|grepl("[B|D-F|H-S|U-Z|b|d-f|h-s|u-z]+", seq_list))) stop("Feature IDs have non-nucleotide characters, but the sequence variant input option was selected. If the rows of your table are OTU IDs, select the option for their database source on the input page.")
    if(config_table[V1=="genomeChoices", V2==get_text("source_choices")[1]]) { ## Greengenes
      seq_results = get_otus_from_seqvar(seq_list, repSeqDir = "~/Documents/MIMOSA2shiny/data/rep_seqs/", repSeqFile = "gg_13_8_99_db.udb", add_agora_names = F, seqID = 0.99) #Run vsearch to get gg OTUs
      species[,seqID:=paste0("seq", 1:nrow(species))]
      print(seq_results)
      samps = names(species)[!names(species) %in% c("OTU", "seqID")]
      new_species = merge(species, seq_results, by = "seqID", all.x=T)
      new_species = new_species[,lapply(.SD, sum), by=dbID, .SDcols = samps]
      setnames(new_species, "dbID", "OTU")
      new_species[is.na(OTU), OTU:=0] #Unassigned
      species = new_species
      mod_list = species[,OTU]
    } else if(config_table[V1=="genomeChoices", V2==get_text("source_choices")[2]]){ ## AGORA
      seq_results = get_otus_from_seqvar(seq_list, repSeqDir = "~/Documents/MIMOSA2shiny/data/blastDB/", repSeqFile = "agora_NCBI_16S.udb", method = "vsearch", file_prefix = "seqtemp", seqID = config_table[V1=="simThreshold", as.numeric(V2)], add_agora_names = T)
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
      species = gg_to_agora(species)
      mod_list = species[!is.na(OTU),OTU]
    } else stop('Model option not implemented')
  }
  ### Now build network from mod_list
  if(config_table[V1=="genomeChoices", V2==get_text("source_choices")[1]]){ #KEGG
    print(mod_list)
    if(config_table[V1=="database", !V2 %in% get_text("database_choices")[c(1, 2, 4)]]){
      stop("Only Greengenes currently implemented")
    }
    if(config_table[V1=="database", V2==get_text("database_choices")[4]]){
      #Get network from metagenome KOs
      species = species[rowSums(species[,names(species) != "KO", with=F]) != 0]
      network_template = fread(paste0(config_table[V1=="kegg_prefix", V2], "/network_template.txt")) ##generate_network_template_kegg(kegg_paths[1], all_kegg = kegg_paths[2:3], write_out = F)
      network = generate_genomic_network(species[,unique(KO)], keggSource = "KeggTemplate", degree_filter = 0, rxn_table = network_template, return_mats = F)
    } else { ##database_choices 1 or 3
      network = get_kegg_network(mod_list, net_path = gg_path)
    }
  } else if(config_table[V1=="genomeChoices", V2==get_text("source_choices")[2]]){ #AGORA
    network = build_species_networks_w_agora(mod_list)
  } else stop('Invalid model format specified')
  if(config_table[V1=="database", V2!=get_text("database_choices")[4]]) species = species[OTU %in% network[,OTU]]
  if(config_table[V1=="netAdd", !is.null(V2)]){
    if(config_table[V1=="netAdd", V2!=F]){
      network = add_to_network(network, netAdd, target_format = config_table[V1=="genomeChoices"])
    }
  }
  if(config_table[V1=="gapfill", V2 != F]){
    #Do stuff
  }
  network = filter_currency_metabolites(network, degree_filter = 30)
  return(list(network, species))
}

#' Add reactions to a network. Will set stoichiometry and copy number to 1 if missing. Format can be either "KO, Rxn, Prod" for reaction IDs with all transformations and correct stoichiometry, or just "KO" but with reaction IDs that are defined in the KEGG network.
#'
#' @import data.table
#' @param network Data.table of taxa, genes and reactions
#' @param addTable Data.table of taxa, genes and/or reactions to add, or generic genes and reactions to be applied to all taxa
#' @param target_format Format of taxa, genes, and/or reactions to add - must be "KEGG" or "Cobra"
#' @param source_format Format of taxa, genes, and/or reactions to add - must be "KEGG" or "Cobra". If NULL, will try to guess
#' @return Expanded network table
#' @examples
#' add_to_network(network, netAddTable)
#' @export
add_to_network = function(network, addTable, target_format, source_format = NULL){
  if(is.null(source_format)) source_format = get_compound_format(network[,unique(Reac)])
  table_type = names(addTable)
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
      addTable[,normalized_copy_number:=1]
    }
  } else if("KO" %in% table_type){
    if(source_format != "KEGG") stop("Currently not implemented, if not KEGG format provide all compound IDs") else {
      full_kegg_table = generate_network_template_kegg("~/Documents/MIMOSA2shiny/data/KEGGfiles/reaction_mapformula.lst", c("~/Documents/MIMOSA2shiny/data/KEGGfiles/reaction_ko.list", "~/Documents/MIMOSA2shiny/data/KEGGfiles/reaction"))
      addTable = merge(addTable, full_kegg_table, by = "KO", all.x = F, all.y = F, allow.cartesian = F)
    }
  } else {
    stop("Invalid format")
  }
  if("remove" %in% table_type){
    removeTable = addTable[remove == T]
    addTable = addTable[remove == F]
  }
  if(!"Species" %in% table_type){
    all_spec = network[,unique(Species)]
    new_net = data.table()
    for(spec in all_spec){ #Propagate to all species
      new_net1 = copy(addTable)
      new_net1[,Species:=spec]
      new_net = rbind(new_net, new_net1)
    }
    addTable = new_net
  }
  network = rbind(network, addTable, fill = T)
  if("remove" %in% table_type){
    if(!all(c("Prod", "Reac") %in% table_type)){ #If just KOs
      network = network[!KO %in% removeTable[,KO]]
    } else {
      for(j in nrow(removeTable)){
        network = network[!(KO==removeTable[j,KO] & Prod==removeTable[j,Prod] & Reac==removeTable[j,Reac])]
      }
    }
  }
  return(network)
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
  config_table[V2=="", V2:=FALSE]
  file_list = as.list(config_table[grepl("file", V1, ignore.case = T)|V1=="metagenome", V2])
  names(file_list) = config_table[grepl("file", V1, ignore.case = T)|V1=="metagenome", V1]
  data_inputs = read_mimosa2_files(file_list, config_table, app = F)
  species = data_inputs$species
  mets = data_inputs$mets

  network_results = build_metabolic_model(species, config_table, gg_path = "~/Documents/MIMOSA2shiny/data/picrustGenomeData/indivModels/")
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

  if(config_table[V1=="metType", V2 !=met_type_choices[1]]){
      #mets = map_to_kegg(mets)
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
  #Order dataset for plotting
  #met_order = var_shares[Species=="Residual"][order(VarShare, increasing = T), metID]
  #var_shares[,metID:=factor(metID, levels = metID)]

  return(list(varShares = var_shares, modelData = cmp_mods[[1]]))
}

