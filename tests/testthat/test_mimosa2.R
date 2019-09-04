library(data.table)
library(RColorBrewer)
library(cowplot)
options(stringsAsFactors = F)
context("MIMOSA2 tests")


test_config_file1 = "test_config_seq_agora.txt"
test_config_file2 = "../../../MIMOSA2shiny/data/exampleData/configs_example_clean.txt"
# config1 = fread(test_config_file1, header = F, fill = T)
# config1[V1=="genomeChoices", V2:=get_text("source_choices")[2]]
# config1[V1=="database", V2:=get_text("database_choices")[1]]
# config1 = config1[!V1 %in% c("modelTemplate", "gapfill")]
# config1 = rbind(config1, data.table(V1 = "vsearch_path", V2 = "~/Documents/MIMOSA2shiny/bin/vsearch"))
# config1 = rbind(config1, data.table(V1 = "metType", V2 = get_text("met_type_choices")[1]))
# write.table(config1, file = test_config_file1, row.names = F, col.names = F, quote=F, sep = "\t")
# 
test_results_normal = function(config_table, file_prefix, make_plots = F){
  expect_silent(check_config_table(config_table))
  config_table = check_config_table(config_table)
  file_list = as.list(config_table[grepl("file", V1, ignore.case = T), V2])
  names(file_list) = config_table[grepl("file", V1, ignore.case = T), V1]
  data_inputs = read_mimosa2_files(file_list, config_table, app = F)
  expect_silent(read_mimosa2_files(file_list, config_table, app = F))
  expect_error(read_mimosa2_files(file_list, config_table, app = T))
  species = data_inputs$species
  mets = data_inputs$mets
  expect_gt(nrow(species), 0)
  expect_gt(nrow(mets), 0)
  if(config_table[V1=="database", !V2 %in% get_text("database_choices")[4:5]]) expect_true("OTU" %in% names(species))
  expect_true("compound" %in% names(mets))
  if(config_table[V1=="metType", V2!=get_text("met_type_choices")[2]]) expect_equal(get_compound_format(mets[,compound]), "KEGG")
  if(config_table[V1=="database", !V2 %in% get_text("database_choices")[4:5]]) expect_equal(nrow(species[is.na(OTU)]), 0) else {
    if(config_table[V1=="database", V2==get_text("database_choices")[4]]){
      expect_equal(nrow(species[is.na(KO)]), 0)
    } else {
      expect_equal(nrow(species[is.na(Gene)]), 0)
      expect_equal(nrow(species[is.na(Sample)]), 0)
      expect_equal(nrow(species[is.na(OTU)]), 0)
    }
  }
  expect_equal(nrow(mets[is.na(compound)]), 0)
  network_results = build_metabolic_model(species, config_table, netAdd = data_inputs$netAdd)
  network = network_results[[1]]
  species = network_results[[2]]
  if(config_table[V1=="database", V2 != get_text("database_choices")[4]]) expect_true("OTU" %in% names(species))
  expect_gt(nrow(species), 0)
  if(config_table[V1=="database", V2 != get_text("database_choices")[4]]) {
    expect_true(all(c("OTU", "KO", "Reac", "Prod", "stoichReac", "stoichProd", "normalized_copy_number") %in% names(network)))
  } else {
    print(head(network))
    expect_true(all(c("KO", "Reac", "Prod", "stoichReac", "stoichProd") %in% names(network)))
  }
  if(config_table[V1=="database", V2 == get_text("database_choices")[5]]) {
    expect_true(all(c("Gene", "OTU", "CountContributedByOTU") %in% names(species)))
  }
  if(config_table[V1=="database", V2 != get_text("database_choices")[4]]) expect_setequal(network[,unique(as.character(OTU))], species[,as.character(OTU)])
  expect_known_output(network, file = paste0(file_prefix, "_net.rda"))
  # if(!is.null(data_inputs$metagenome) & config_table[V1=="database", V2!=get_text("database_choices")[4]]){
  #   metagenome_network = build_metabolic_model(data_inputs$metagenome, config_table, netAdd = data_inputs$netAdd)
  #   expect_equal(metagenome_network[,unique(OTU)], "TotalMetagenome")
  #   expect_true(all(c("KO", "Reac", "Prod", "stoichReac", "stoichProd", "normalized_copy_number") %in% names(metagenome_network)))
  # }
  
  if(config_table[V1=="metType", V2 ==get_text("met_type_choices")[2]]){ #Assume it is KEGG unless otherwise specified
    expect_warning(map_to_kegg(mets))
    mets = map_to_kegg(mets)
    expect_equal(get_compound_format(mets[,compound]), "KEGG")
  }

  #Get CMP scores
  if("rxnEdit" %in% config_table[,V1]){
    rxn_param = T
    cat("Will refine reaction network\n")
  } else rxn_param = F
  if("rankBased" %in% config_table[,V1]){
    rank_based = T
    cat("Will use rank-based/robust regression\n")
    if("rank_type" %in% config_table[,V1]){
      rank_type = config_table[V1=="rank_type", V2]
    } else {
      rank_type = "rfit"
    }
    cat(paste0("Regression type is ", rank_type, "\n"))
  } else rank_based = F
  if(config_table[V1=="database", V2==get_text("database_choices")[4]]){
    no_spec_param = T
    humann2_param = F
  } else if(config_table[V1=="database", V2==get_text("database_choices")[5]]){
    no_spec_param = F
    humann2_param = T
  } else {
    no_spec_param = F
    humann2_param = F
  }
  if("manualAGORA" %in% config_table[,V1]){
    agora_param = T
    cat("Manual AGORA models\n")
  } else {
    agora_param = F
  }
  # if("revRxns" %in% config_table[,V1]){ #Whether to add reverse of reversible-annotated rxns - mainly for agora networks
  #   network = add_rev_rxns(network, sameID = T) # Give reverse the same rxn ID
  #   cat("Will add reverse of reversible reactions\n")
  # }
  # if("refine" %in% config_table[,V1]){
  #   network = refine_rev_rxns(network, )
  # }
  if("met_transform" %in% config_table[,V1]){
    met_transform = config_table[V1=="met_transform", V2]
    cat(paste0("Will transform metabolite values, transform is ", met_transform))
  } else met_transform = ""
  if("score_transform" %in% config_table[,V1]){
    score_transform = config_table[V1=="score_transform", V2]
    cat(paste0("Will transform CMP values, transform is ", score_transform))
  } else score_transform = ""
  
  if(config_table[V1=="compare_only", V2==T]){
    compare_only = T
  } else {
    compare_only = F
  }
  
  #indiv_cmps = get_cmp_scores_kos(species, network) #Use KO abundances instead of species abundances to get cmps
  mets_melt = melt(mets, id.var = "compound", variable.name = "Sample")
  if(met_transform != ""){
    mets_melt = transform_mets(mets_melt, met_transform)
  }
  if(rxn_param){
    cmp_mods =  fit_cmp_net_edit(network, species, mets_melt, manual_agora = agora_param, rank_based = rank_based)
    network = cmp_mods[[3]] #Revised network
    indiv_cmps = cmp_mods[[4]]
    #Will have to report nice summary of rxns removed, rxns direction switched, etc
  } else {
    indiv_cmps = get_species_cmp_scores(species, network, normalize = !rxn_param, leave_rxns = rxn_param, manual_agora = agora_param, kos_only = no_spec_param, humann2 = humann2_param)
    if(score_transform != ""){
      indiv_cmps = transform_cmps(indiv_cmps, score_transform)
    }
    indiv_cmps = indiv_cmps[compound %in% mets[,compound]]
    cmp_mods = fit_cmp_mods(indiv_cmps, mets_melt, rank_based = rank_based, rank_type = rank_type)
  }
  
  expect_gt(nrow(indiv_cmps), 0)
  print(head(indiv_cmps))
  expect_setequal(names(indiv_cmps), c("Species", "Sample","compound", "CMP")) #, "NumSynthGenes", "NumSynthSpecies","NumSynthSpecGenes", "NumDegGenes","NumDegSpecies","NumDegSpecGenes"))
  #expect_setequal(indiv_cmps[,unique(Sample)], names(mets)[names(mets) != "compound"])
  expect_true(indiv_cmps[,all(is.numeric(CMP))])
  expect_equal(nrow(indiv_cmps[is.na(compound)]), 0)
  expect_equal(nrow(indiv_cmps[is.na(Species)]), 0)
  expect_equal(nrow(indiv_cmps[is.na(Sample)]), 0)
  
  print(cmp_mods)
  #This one is no longer true if there are metabolites with too few nonzero CMPs
  expect_lte(nrow(cmp_mods[[1]]), length(intersect(mets[,compound], indiv_cmps[,unique(compound)])))
  expect_gt(nrow(cmp_mods[[1]]), 0)
  expect_equal(nrow(cmp_mods[[1]])*indiv_cmps[,length(unique(Sample))], nrow(cmp_mods[[2]]))
  expect_true(cmp_mods[[2]][,all(is.numeric(Resid))])
  expect_equal(nrow(cmp_mods[[2]][is.na(Resid)]), 0)
  
  if(!compare_only & !no_spec_param){ #Option to skip contributions
    if(!humann2_param){
      spec_dat = melt(species, id.var = "OTU", variable.name = "Sample")[,list(value/sum(value), OTU), by=Sample] #convert to relative abundance
      bad_spec = spec_dat[,list(length(V1[V1 != 0])/length(V1), max(V1)), by=OTU]
      bad_spec = bad_spec[V1 < 0.1 & V2 < 0.1, OTU] #Never higher than 10% and absent in at least 90% of samples
      print(bad_spec)
    } else bad_spec = NULL
    if("signifThreshold" %in% config_table[,V1]){
      signifThreshold = config_table[V1 == "signifThreshold", as.numeric(V2)]
    } else {
      signifThreshold = 0.2
    }
    time1 = Sys.time()
    var_shares = calculate_var_shares(indiv_cmps, met_table = mets_melt, model_results = cmp_mods, config_table = config_table, species_merge = bad_spec, signif_threshold = signifThreshold)
    cat(paste0("Contribution calculation time: ", Sys.time() - time1, "\n"))
  } else {
    var_shares = NULL
  }
  expect_known_output(var_shares, file = paste0(file_prefix, "_var_shares.rda"))
  summary_cmps = get_cmp_summary(species, network, normalize = !rxn_param, manual_agora = agora_param, kos_only = no_spec_param, 
                                 humann2 = humann2_param, met_subset = mets[,compound], contrib_sizes = var_shares)
  expect_gt(nrow(summary_cmps$CompLevelSummary), 0)
  cmp_mods[[1]] = merge(cmp_mods[[1]], summary_cmps$CompLevelSummary, by = "compound", all.x = T)
  #Add species/rxn info
  
  if(length(summary_cmps) > 1) var_shares = merge(var_shares, summary_cmps$SpeciesLevelSummary, by = c("compound", "Species"), all.x = T)
  cmp_mods[[1]][,compound:=as.character(compound)]
  indiv_cmps[,compound:=as.character(compound)]
  if(!is.null(var_shares)){
    var_shares[,compound:=as.character(compound)]
    var_shares[,Species:=as.character(Species)]
    var_shares[,MetaboliteName:=met_names(as.character(compound))]
    var_shares[is.na(MetaboliteName), MetaboliteName:=compound]
    var_shares = var_shares[,list(compound, MetaboliteName, Rsq, VarDisp, ModelPVal, Slope, Intercept, Species, VarShare, PosVarShare, NumSynthGenes, SynthGenes, NumDegGenes, DegGenes)]
  }
  if(make_plots){
    CMP_plots = plot_all_cmp_mets(cmp_table = indiv_cmps, met_table = mets_melt, mod_results = cmp_mods[[1]])
    
    if(!compare_only){
      comp_list = var_shares[!is.na(VarShare), unique(as.character(compound))]
      comp_list = comp_list[!comp_list %in% var_shares[Species == "Residual" & VarShare == 1, as.character(compound)]]
      all_contrib_taxa = var_shares[compound %in% comp_list & !is.na(VarShare) & Species != "Residual", as.character(unique(Species))]
      getPalette = colorRampPalette(brewer.pal(12, "Paired"))
      if(var_shares[compound %in% comp_list & Species != "Residual", length(unique(Species[VarShare != 0])), by=compound][,any(V1 > 10)]){
        contrib_color_palette = c("gray", getPalette(length(all_contrib_taxa))) #"black",  #Work with plotting function filters
        names(contrib_color_palette) = c( "Other", all_contrib_taxa) #"Residual",
      } else {
        contrib_color_palette = getPalette(length(all_contrib_taxa)) #"black", 
        names(contrib_color_palette) = all_contrib_taxa #"Residual",
      }
      met_contrib_plots = lapply(comp_list, function(x){
        print(x)
        if(is.na(met_names(x))){
          met_id = x
        } else { met_id = met_names(x)}
        plot_contributions(var_shares, met_id, metIDcol = "MetaboliteName", color_palette = contrib_color_palette, include_residual = F, merge_threshold = 0.01)
      })
      #Contribution Legend
      leg_dat = data.table(V1 = factor(names(contrib_color_palette), levels = c(all_contrib_taxa, "Other"))) #, "Residual"
      setnames(leg_dat, "V1", "Contributing Taxa")
      print(leg_dat)
      legend_plot = ggplot(leg_dat, aes(fill = `Contributing Taxa`, x=`Contributing Taxa`)) + geom_bar() + scale_fill_manual(values = contrib_color_palette, name = "Contributing Taxa")# + theme(legend.text = element_text(size = 10))
      contrib_legend = tryCatch(get_legend(legend_plot), error = function(){ return(NULL)}) 
    } else {
      met_contrib_plots = NULL
      contrib_legend = NULL
    }
  }
  if(config_table[V1=="genomeChoices", V2 != get_text("source_choices")[1]]){
    network[,KEGGReac:=agora_kegg_mets(Reac)]
    network[,KEGGProd:=agora_kegg_mets(Prod)]
  } 
  analysis_summary = get_analysis_summary(input_species = data_inputs[[1]], species = species, mets = mets, network = network, indiv_cmps = indiv_cmps, cmp_mods = cmp_mods, var_shares = var_shares, config_table = config_table)
  
  #Write a test that the contributinos add up to the R-squared
  if(!is.null(var_shares)){
    contribs_all = var_shares[Species != "Residual",sum(VarShare), by=list(compound, Rsq)]
    expect_true(contribs_all[,all(abs(V1-Rsq) < 10e-10)])
    contribs_all2 = var_shares[,sum(VarShare), by=compound] #Including residual
    expect_true(contribs_all2[,all(abs(V1-1) < 10e-10)])
  }
  #Get shared legend
  # if(!config_table[V1 == "compare_only", identical(V2, TRUE)]){
  #   for(i in 1:length(met_contrib_plots)){
  #     print(comp_list[i])
  #     if(!is.null(met_contrib_plots[[i]])){
  #       if(!exists("contrib_legend")){ #Get legend from first non-nulll compound
  #         contrib_legend = get_legend(met_contrib_plots[[i]])
  #       }
  #     }
  #   }
  # }
  #plot_summary_contributions(plotData, include_zeros = T, remove_resid_rescale = F) - add this next
  expect_output(run_mimosa2(config_table))
}


test_that("Seq var -> AGORA species", {
  config1 = fread(test_config_file1, header = F, fill = T)
  test_results_normal(config1, file_prefix = "test_seq_agora")
  config1 = rbind(config1, data.table(V1 = "rankBased", V2 = T))
  test_results_normal(config1, "test_seq_agora_rank")
})


test_that("Seq var -> Greengenes OTUs, species-rxn KEGG mods", {
  config1 = fread(test_config_file1, header = F, fill = T)
  config1[V1=="genomeChoices", V2:=get_text("source_choices")[1]]
  config1[V1=="netAdd", V2:="test_netAdd_species_rxns_KEGG.txt"]
  test_results_normal(config1, file_prefix = "test_seq_gg")
  config1 = rbind(config1, data.table(V1 = "rankBased", V2 = T))
  test_results_normal(config1, "test_seq_agora_rank")
})

test_that("GG OTUs -> network", {
  config1 = fread(test_config_file1, header = F, fill = T)
  config1[V1=="database", V2:=get_text("database_choices")[2]]
  config1[V1=="genomeChoices", V2:=get_text("source_choices")[1]]
  config1[V1=="file1", V2:="test_gg.txt"]
  config1[V1=="netAdd", V2:="test_netAdd_species_rxns_KEGG.txt"]
  test_results_normal(config1, file_prefix = "test_otus_gg")
  config1 = rbind(config1, data.table(V1 = "rankBased", V2 = T))
  test_results_normal(config1, "test_seq_agora_rank")
  
})


test_that("Greengenes OTUs -> AGORA species", {
  config1 = fread(test_config_file1, header = F, fill = T)
  config1[V1=="genomeChoices", V2:=get_text("source_choices")[2]]
  config1[V1=="database", V2:=get_text("database_choices")[2]]
  config1[V1=="file1", V2:="test_gg.txt"]
  config1[V1=="netAdd", V2:="test_netAdd_species_rxns_AGORA.txt"]
  test_results_normal(config1, file_prefix = "gg_agora_addAgora")
  config1 = rbind(config1, data.table(V1 = "rankBased", V2 = T))
  test_results_normal(config1, "test_seq_agora_rank")
  
})

test_that("Greengenes OTUs -> AGORA species, KEGG add", {
  config1 = fread(test_config_file1, header = F, fill = T)
  config1[V1=="netAdd", V2:="test_netAdd_species_rxns_KEGG.txt"]
  test_results_normal(config1, file_prefix = "gg_agora_addKEGG")
  config1 = rbind(config1, data.table(V1 = "rankBased", V2 = T))
  test_results_normal(config1, "test_seq_agora_rank")
  
})


#Decide what to do about this
# test_that("Greengenes OTUs -> AGORA species, mixed add", {
#   config1 = fread(test_config_file1, header = F, fill = T)
#   config1[V1=="netAdd", V2:="test_netAdd_species_rxns_KEGG2.txt"]
#   expect_error(test_results_normal(config1, file_prefix = "gg_agora_addMixed"))
# })


test_that("Non-KEGG metabolites", {
  config1 = fread(test_config_file1, header = F, fill = T)
  config1[V1=="file2", V2:="test_mets_names.txt"]
  config1[V1=="metType", V2:=get_text("met_type_choices")[2]]
  test_results_normal(config1, file_prefix = "mets_names")
  config1 = rbind(config1, data.table(V1 = "rankBased", V2 = T))
  test_results_normal(config1, "test_seq_agora_rank")
  
})


test_that("Species-gene modifications work, KEGG", {
  config1 = fread(test_config_file1, header = F, fill = T)
  config1[V1=="file1", V2:="test_gg.txt"]
  config1[V1=="genomeChoices", V2:=get_text("source_choices")[1]]
  config1[V1=="database", V2:=get_text("database_choices")[2]]
  config1[V1=="netAdd", V2:="test_netAdd_species_genes_KEGG.txt"]
  test_results_normal(config1, file_prefix = "test_gg_addSpecGenes")
  config1 = rbind(config1, data.table(V1 = "rankBased", V2 = T))
  test_results_normal(config1, "test_seq_agora_rank")
  
})

test_that("Gene modifications work, KEGG", {
  config1 = fread(test_config_file1, header = F, fill = T)
  config1[V1=="file1", V2:="test_gg.txt"]
  config1[V1=="genomeChoices", V2:=get_text("source_choices")[1]]
  config1[V1=="database", V2:=get_text("database_choices")[2]]
  config1[V1=="netAdd", V2:="test_netAdd_genes_KEGG.txt"]
  test_results_normal(config1, file_prefix = "test_gg_addGenes")
  config1 = rbind(config1, data.table(V1 = "rankBased", V2 = T))
  test_results_normal(config1, "test_seq_agora_rank")
  
})

test_that("Rxn modifications work, KEGG", {
  config1 = fread(test_config_file1, header = F, fill = T)
  config1[V1=="file1", V2:="test_gg.txt"]
  config1[V1=="genomeChoices", V2:=get_text("source_choices")[1]]
  config1[V1=="database", V2:=get_text("database_choices")[2]]
  config1[V1=="netAdd", V2:="test_netAdd_rxns_KEGG.txt"]
  test_results_normal(config1, file_prefix = "test_gg_addRxns")
})

test_that("Rxn modifications work, AGORA", {
  config1 = fread(test_config_file1, header = F, fill = T)
  config1[V1=="file1", V2:="test_seqs.txt"]
  config1[V1=="genomeChoices", V2:=get_text("source_choices")[2]]
  config1[V1=="database", V2:=get_text("database_choices")[1]]
  config1[V1=="netAdd", V2:="test_netAdd_rxns_AGORA.txt"]
  test_results_normal(config1, file_prefix = "test_agora_addRxns")
  config1 = rbind(config1, data.table(V1 = "rankBased", V2 = T))
  test_results_normal(config1, "test_seq_agora_rank")
  
})

test_that("Gene modifications work, AGORA", {
  config1 = fread(test_config_file1, header = F, fill = T)
  config1[V1=="file1", V2:="test_seqs.txt"]
  config1[V1=="genomeChoices", V2:=get_text("source_choices")[2]]
  config1[V1=="database", V2:=get_text("database_choices")[1]]
  config1[V1=="netAdd", V2:="test_netAdd_genes_AGORA.txt"]
  test_results_normal(config1, file_prefix = "test_agora_addGenes")
  config1 = rbind(config1, data.table(V1 = "rankBased", V2 = T))
  test_results_normal(config1, "test_seq_agora_rank")
  
})



test_that("species-gene modifications work, AGORA", {
  config1 = fread(test_config_file1, header = F, fill = T)
  config1[V1=="file1", V2:="test_seqs.txt"]
  config1[V1=="genomeChoices", V2:=get_text("source_choices")[2]]
  config1[V1=="database", V2:=get_text("database_choices")[1]]
  config1[V1=="netAdd", V2:="test_netAdd_species_genes_AGORA.txt"]
  test_results_normal(config1, file_prefix = "test_agora_addSpecGenes")
  config1 = rbind(config1, data.table(V1 = "rankBased", V2 = T))
  test_results_normal(config1, "test_seq_agora_rank")
  
})


test_that("Metagenome option works", {
  config1 = fread(test_config_file1, header = F, fill = T)
  config1[V1=="database", V2:=get_text("database_choices")[4]]
  config1[V1=="genomeChoices", V2:=get_text("source_choices")[1]]
  config1[V1 == "file1", V2:="test_metagenome.txt"]
  #config1 = rbind(config1, data.table(V1=c("metagenome", "metagenome_format"), V2 = c("test_metagenome.txt", get_text("metagenome_options")[1])))
  config1[V1=="file2", V2:="test_metagenome_mets.txt"]
  config1[V1=="netAdd", V2:="test_netAdd_rxns_KEGG.txt"]
  test_results_normal(config1, file_prefix = "test_metagenome")
  config1 = rbind(config1, data.table(V1 = "rankBased", V2 = T))
  test_results_normal(config1, "test_seq_agora_rank")
  
})

test_that("Metagenome stratified option works", {
  config1 = fread(test_config_file1, header = F, fill = T)
  config1[V1=="database", V2:=get_text("database_choices")[5]]
  config1[V1=="genomeChoices", V2:=get_text("source_choices")[1]]
  config1[V1 == "file1", V2:="test_contributions.txt"]
  config1[V1=="file2", V2:="test_mets.txt"]
  config1[V1=="netAdd", V2:="test_netAdd_rxns_KEGG.txt"]
  test_results_normal(config1, file_prefix = "test_metagenome_stratified")
  config1 = rbind(config1, data.table(V1 = "rankBased", V2 = T))
  test_results_normal(config1, "test_seq_agora_rank")
  
})

test_that("Rank contribution timing", {
  config1 = fread(test_config_file1, header = F, fill = T)
  config1[V1=="database", V2:=get_text("database_choices")[2]]
  config1[V1=="genomeChoices", V2:=get_text("source_choices")[1]]
  config1[V1=="file1", V2:="test_gg.txt"]
  config1[V1=="netAdd", V2:="test_netAdd_species_rxns_KEGG.txt"]
  config1 = rbind(config1, data.table(V1 = "rankBased", V2 = T))
  foo = run_mimosa2(config1)
})
