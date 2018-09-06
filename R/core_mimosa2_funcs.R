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
#' @param varShares Dataset of contributions
#' @param include_zeros Whether to plot taxa that do not have any contribution
#' @return plot of contributions
#' @examples
#' plot_summary_contributions(varShares)
#' @export
plot_summary_contributions = function(varShares, include_zeros = T){
  varShares[,metID:=met_names(as.character(compound))]
  if(!include_zeros){
    varShares = varShares[VarShare != 0]
  }
  plot1 = ggplot(varShares, aes(x=metID, y = as.character(Species), fill = VarShare)) + geom_tile() + theme_minimal() +
    theme(axis.text.x = element_text(angle=90, hjust=0, vjust =0.5), axis.line = element_blank()) +
    scale_fill_gradient(low = brewer.pal(9,"Blues")[1], high = brewer.pal(9, "Blues")[9]) + xlab("Metabolite") + ylab("Taxon")
  return(plot1)
}

#' Read in files for a MIMOSA 2 analysis
#'
#' @import data.table
#' @param file_list List of shiny file inputs to load
#' @param configTable Table of configuration parameters
#' @return list of processed abundance data.tables
#' @examples
#' read_mimosa2_files(input_file_list, config_table)
#' @export
read_mimosa2_files = function(file_list, configTable){
  if(configTable[V1=="metagenome_use", V2 == F]){
    species = fread(file_list[["file1"]]$datapath)
    species = spec_table_fix(species)
  } else {
    metagenome = fread(file_list[["metagenome"]]$datapath)
    species = get_species_from_metagenome(metagenome)
  }
  mets = fread(file_list[["file2"]]$datapath)
  met_nonzero_filt = ifelse(configTable[V1=="metNzeroFilter", is.numeric(V2)], configTable[V1=="metNzeroFilter", V2], 5)
  mets = met_table_fix(mets, met_nonzero_filt)
  #Filter species using default abundance values
  if(configTable[V1=="specNzeroFrac", is.numeric(V2)]){
    species = filter_species_abunds(species, filter_type = "fracNonzero", configTable[V1=="specNzeroFrac", V2])
  } else { #Use default values
    species = filter_species_abunds(species, filter_type = "fracNonzero")
  }
  if(configTable[V1=="specMinMean", is.numeric(V2)]){
    species = filter_species_abunds(species, filter_type = "mean", configTable[V1=="specMinMean", V2])
  } else { #Use default values
    species = filter_species_abunds(species, filter_type = "mean")
  }
  species = filter_species_abunds(species, filter_type = "fracNonzero")
  species = filter_species_abunds(species, filter_type = "mean")
  shared_samps = intersect(names(species), names(mets))
  if(length(shared_samps) < 2) stop("Sample IDs don't match between species and metabolites")
  species = species[,c("OTU", shared_samps), with=F]
  mets = mets[,c("compound", shared_samps), with=F]
  dat_list = list(species, mets)
  names(dat_list) = c("species", "mets")
  for(extraFile in c("metagenome", "geneAddFile", "netAddFile")){
    print(file_list[[extraFile]])
    if(length(file_list[[extraFile]]) > 0 ){
      if(file_list[[extraFile]] != F){
        dat = fread(file_list[[extraFile]]$datapath)
        dat_list[[length(dat_list)+1]] = dat
        names(dat_list[[length(dat_list)]]) = extraFile
      }
    }
  }
  return(dat_list)
}


#' Build species-specific metabolic model for MIMOSA2 analysis
#'
#' @import data.table
#' @param species Table of species/taxon abundances
#' @param config_table Data.table of input files and settings for MIMOSA
#' @return Data.table of network model of genes and reactions for each species/taxon
#' @examples
#' build_metabolic_model(config_table)
#' @export
build_metabolic_model = function(species, config_table, geneAdd = NULL, netAdd = NULL){
  if(config_table[V1=="genomeChoices", V2==get_text("source_choices")[1]]){
    if(config_table[V1=="database", V2==get_text("database_choices")[1]]){
      seq_list = species[,OTU]
      if(any(grepl("[0-9]+", seq_list)|grepl("[B|D-F|H-S|U-Z|b|d-f|h-s|u-z]+", seq_list))) stop("Feature IDs have non-nucleotide characters, but the sequence variant input option was selected. If the rows of your table are OTU IDs, select the option for their database source on the input page.")
      species_table = get_otus_from_seqvar(seq_list, repSeqDir = "~/Documents/MIMOSA2shiny/data/rep_seqs/", repSeqFile = "gg_13_5.fasta.gz", add_agora_names = F, seqID = 0.97) #Run vsearch to get gg OTUs
    } else if(config_table[V1=="database", V2!= get_text("database_choices")[2]]){
      stop("Only Greengenes currently implemented")
    }
    if(!is.null(geneAdd)){
      contribution_table = generate_contribution_table_using_picrust(species, picrust_norm_file = "data/picrustGenomeData/16S_13_5_precalculated.tab.gz", picrust_ko_table_directory ="data/picrustGenomeData/indivGenomes/", picrust_ko_table_suffix = "_genomic_content.tab")
      contribution_table = add_genes_to_contribution_table(contribution_table, geneAdd)
      network = build_generic_network(contribution_table, kegg_paths = c("data/KEGGfiles/reaction_mapformula.lst", "data/KEGGfiles/reaction_ko.list", "data/KEGGfiles/reaction"))
    } else { #Just load in preprocessed
      network = get_kegg_network(species, net_path = "data/picrustGenomeData/indivModels/")
    }
  } else if(config_table[V1=="genomeChoices", V2==get_text("source_choices")[2]]){
    if(any(grepl("[0-9]+", species[,OTU])|grepl("[B|D-F|H-S|U-Z|b|d-f|h-s|u-z]+", species[,OTU]))) stop("Feature IDs have non-nucleotide characters, but the sequence variant input option was selected. If the rows of your table are OTU IDs, select the option for their database source on the input page.")
    network_results = build_species_networks_w_agora(species, config_table[V1=="database", V2], config_table[V1=="closest", V2], config_table[V1=="simThreshold", V2])
    species = network_results[[1]]
    network = network_results[[2]]
    species = species[OTU %in% network[,OTU]]
  }
  if(!is.null(netAdd)){
    network = add_rxns_to_network(network, netAdd)
    #This will need to map between metabolite IDs possibly
  }
  if(config_table[V1=="gapfill", V2 != F]){
    #Do stuff
  }
  return(list(network, species))
}

#' Add reactions to a network
#'
#' @import data.table
#' @param network Data.table of taxa, genes and reactions
#' @param netAdd Data.table of taxa, genes and reactions to add, or generic genes and reactions to be applied to all taxa
#' @return Expanded network table
#' @examples
#' add_rxns_to_network(network, netAddTable)
#' @export
add_rxns_to_network = function(network, netAddTable){
  
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
  data_inputs = read_mimosa2_files(config_table)
  species = data_inputs$species
  mets = data_inputs$mets

  if(!is.null(data_inputs$metagenome)){
    #Metagenome data
    #Implement this later
    metagenome_data = get_species_from_metagenome(data_inputs$metagenome)
    species2 = metagenome_data[[1]]
    metagenome_network = metagenome_data[[2]]
      #Metagenome data
      #Implement this later
  }
  network = build_metabolic_model(species, config_table)
  if(config_table[V1=="metType", V2 !=met_type_choices[1]]){
      #mets = map_to_kegg(mets)
      #Implement this later
  }
  indiv_cmps = get_species_cmp_scores(species, network)
  mets_melt = melt(mets, id.var = "compound", variable.name = "Sample")
  cmp_mods = fit_cmp_mods(indiv_cmps, mets_melt)
  indiv_cmps = add_residuals(indiv_cmps, cmp_mods[[1]], cmp_mods[[2]])
  var_shares = calculate_var_shares(indiv_cmps)
  return(list(varShares = var_shares, modelData = cmp_mods[[1]]))
}

