### MIMOSA linear model w/residual, then get contributions

#' Fits model scaling total CMPs to metabolite concentrations
#'
#' @import data.table
#' @param species_cmps Table of species contribution abundances
#' @param fake_mets_melt Table of metabolite concentrations
#' @return List of 2 data.tables - one with model summary results, one with model residuals
#' @examples
#' fit_cmp_mods(species_cmps, met_data)
#' @export
fit_cmp_mods = function(species_cmps, fake_mets_melt){
  tot_cmps = species_cmps[,sum(value), by=list(compound, Sample)]
  tot_cmps = merge(tot_cmps, fake_mets_melt[,list(compound, Sample, value)], by = c("compound", "Sample"))
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
  if(!"Resid" %in% names(resid_dat)) setnames(resid_dat, "Resid", "newValue")
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
plot_contributions = function(varShares, metabolite, include_zeros = F){
  plot_dat = varShares[compound==metabolite]
  if(!include_zeros) plot_dat = plot_dat[VarShare != 0]
  ggplot(plot_dat,  aes(y=VarShare, x = Species, fill = Species)) + geom_bar(stat = "identity") + scale_fill_viridis(option = "plasma", discrete = T) + geom_abline(intercept = 0, slope = 0, linetype = 2) + theme(strip.background = element_blank(), axis.ticks.y = element_blank(), axis.text.x = element_text(size=8), axis.text.y = element_blank(), legend.title = element_blank(), strip.text = element_blank(), axis.title.y = element_blank(), panel.spacing = unit(0.15, "inches"), plot.margin = margin(0.2, 0.4, 0.3, 0.1, "inches")) + ylab("Contribution to variance") + coord_flip()#

}


#' Run a MIMOSA 2 analysis
#'
#' @import data.table
#' @param species_file Path to species abundances
#' @param met_file Path to metabolite abundances
#' @param config_file Path to config file
#' @return Scaling model and variance contribution results
#' @examples
#' run_mimosa2(species_file, met_file, config_file)
#' @export
run_mimosa2 = function(species_file, met_file, config_file){
  #process arguments
  species = fread(species_file)
  species = spec_table_fix(species)
  mets = fread(met_file)
  met_col_name = names(mets)[names(mets) %in% c("compound", "KEGG", "Compound", "metabolite", "Metabolite")]
  if(length(met_col_name) != 1) stop("Ambiguous metabolite ID column name, must be one of Compound/KEGG/Metabolite")
  setnames(mets, met_col_name, "compound")
  shared_samps = intersect(names(species), names(mets))
  if(length(shared_samps) < 2) stop("Sample IDs don't match between species and metabolites")
  species = species[,c("OTU", shared_samps), with=F]
  mets = mets[,c("compound", shared_samps), with=F]
  configs = fread(config_file, header = F, fill = T, sep = "\t")
  if("metagenome" %in% configs[,V1]){
      #Metagenome data
      #Implement this later
  }
  if(configs[V1=="genomeChoices", V2=="Assign KOs with PICRUSt"]){
      if(configs[V1=="database", V2=="Sequence variants (recommended for AGORA)"]){
        seq_list = species[,OTU]
        species_table = get_otus_from_seqvar(seq_list, repSeqDir = "~/Documents/MIMOSA2shiny/data/rep_seqs/", repSeqFile = "gg_13_5.fasta.gz", add_agora_names = F, seqID = configs[V1=="simThreshold", V2]) #Run vsearch to get gg OTUs
      } else if(configs[V1=="database", V2 != "Greengenes 13_5 or 13_8"]){
        stop("Only Greengenes currently implemented")
      }
      contribution_table = generate_contribution_table_using_picrust(species, picrust_norm_file = "data/picrustGenomeData/16S_13_5_precalculated.tab.gz", picrust_ko_table_directory ="data/picrustGenomeData/indivGenomes/", picrust_ko_table_suffix = "_genomic_content.tab")
      contribution_table = contribution_table[contribution != 0]
      if("geneAdd" %in% configs[,V1]){
        contribution_table = add_genes_to_contribution_table(contribution_table, geneAddFile)
      }
    }
    if(configs[V1=="modelTemplate", V2=="Generic KEGG metabolic model"]){
      kegg_prefix = configs[V1=='kegg_prefix', V2]
      network = build_generic_network(contribution_table, kegg_paths = c(paste0(kegg_prefix, "reaction_mapformula.lst"), paste0(kegg_prefix, "reaction_ko.list"), paste0(kegg_prefix, "reaction")))
    }else{
      network_results = build_species_networks_w_agora(species, configs[V1=="database", V2], closest = configs[V1=="closest", V2] != "", simThreshold = configs[V1=="simThreshold", V2], filterAbund = T)
      species = network_results[[1]]
      network = network_results[[2]]
      species = species[OTU %in% network[,OTU]]
    }
    if("netAdd" %in% configs[,V1]){
      network = add_rxns_to_network(network, configs[V1=="netAddFile", V2])
      #This will need to map between metabolite IDs possibly
    }
    if("gapFill" %in% configs[,V1]){
      #Do stuff
    }
    if(metType!="KEGG Compound IDs"){
      #mets = map_to_kegg(mets)
      #Implement this later
    }
    indiv_cmps = get_species_cmp_scores(species, network)
    if(configs[V1=="modelTemplate", V2=="AGORA metabolic models (recommended)"]){ #Switch to KEGG IDs at this point
      kegg_mapping = fread("../MIMOSA2shiny/data/KEGGfiles/AGORA_KEGG_met_mappings.txt") #This should probably be part of the package data
      indiv_cmps = merge(indiv_cmps, kegg_mapping, by.x = "compound", by.y = "met", all.x=T)
      indiv_cmps = indiv_cmps[,sum(CMP), by=list(Species, KEGG, Sample)] #Check that this makes sense
      setnames(indiv_cmps, c("KEGG", "V1"), c("compound", "CMP"))
    }
    tot_cmps = indiv_cmps[,sum(CMP), by=list(compound, Sample)]
    setnames(tot_cmps, "V1", "value")
    mets_melt = melt(mets, id.var = "compound", variable.name = "Sample")
    cmp_mods = fit_cmp_mods(tot_cmps, mets_melt)
    indiv_cmps = add_residuals(indiv_cmps, cmp_mods[[1]], cmp_mods[[2]])
    var_shares = calculate_var_shares(indiv_cmps)
    return(list(varShares = var_shares, modelData = cmp_mods[[1]]))
}

