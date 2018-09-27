library(data.table)
options(stringsAsFactors = F)
context("MIMOSA2 tests")


test_config_file1 = "test_config_seq_agora.txt"
config1 = fread(test_config_file1, header = F, fill = T)


test_results_normal = function(config_table, file_prefix){
  expect_silent(check_config_table(config_table))
  config_table = check_config_table(config_table)
  file_list = as.list(config_table[grepl("file", V1, ignore.case = T)|V1=="metagenome", V2])
  names(file_list) = config_table[grepl("file", V1, ignore.case = T)|V1=="metagenome", V1]
  data_inputs = read_mimosa2_files(file_list, config_table, app = F)
  expect_silent(read_mimosa2_files(file_list, config_table, app = F))
  expect_error(read_mimosa2_files(file_list, config_table, app = T))
  species = data_inputs$species
  mets = data_inputs$mets
  expect_gt(nrow(species), 0)
  expect_gt(nrow(mets), 0)
  if(config_table[V1=="database", V2 != get_text("database_choices")[4]]) expect_true("OTU" %in% names(species))
  expect_true("compound" %in% names(mets))
  if(config_table[V1=="metType", V2!=get_text("met_type_choices")[2]]) expect_equal(get_compound_format(mets[,compound]), "KEGG")
  if(config_table[V1=="database", V2 != get_text("database_choices")[4]]) expect_equal(nrow(species[is.na(OTU)]), 0) else {
    expect_equal(nrow(species[is.na(KO)]), 0)
  }
  expect_equal(nrow(mets[is.na(compound)]), 0)
  network_results = build_metabolic_model(species, config_table)
  network = network_results[[1]]
  species = network_results[[2]]
  if(config_table[V1=="database", V2 != get_text("database_choices")[4]]) expect_true("OTU" %in% names(species))
  expect_gt(nrow(species), 0)
  if(config_table[V1=="database", V2 != get_text("database_choices")[4]]) {
    expect_true(all(c("OTU", "KO", "Reac", "Prod", "stoichReac", "stoichProd", "normalized_copy_number") %in% names(network)))
  } else {
    expect_true(all(c("KO", "Reac", "Prod", "stoichReac", "stoichProd", "normalized_copy_number") %in% names(network)))
  }
  if(config_table[V1=="database", V2 != get_text("database_choices")[4]]) expect_setequal(network[,unique(as.character(OTU))], species[,as.character(OTU)])
  expect_known_output(network, file = paste0(file_prefix, "_net.rda"))
  if(!is.null(data_inputs$metagenome) & config_table[V1=="database", V2!=get_text("database_choices")[4]]){
    metagenome_network = build_metabolic_model(data_inputs$metagenome, config_table)
    expect_equal(metagenome_network[,unique(OTU)], "TotalMetagenome")
    expect_true(all(c("KO", "Reac", "Prod", "stoichReac", "stoichProd", "normalized_copy_number") %in% names(metagenome_network)))
  }
  if(config_table[V1=="metType", V2 ==get_text("met_type_choices")[2]]){ #Assume it is KEGG unless otherwise specified
    expect_warning(map_to_kegg(mets))
    mets = map_to_kegg(mets)
    expect_equal(get_compound_format(mets[,compound]), "KEGG")
  }
  if(config_table[V1=="database", V2==get_text("database_choices")[4]]){
    indiv_cmps = get_cmp_scores_kos(species, network) #Use KO abundances instead of species abundances to get cmps
  } else {
    indiv_cmps = get_species_cmp_scores(species, network)
  }
  expect_gt(nrow(indiv_cmps), 0)
  expect_setequal(names(indiv_cmps), c("Species", "compound", "Sample", "CMP"))
  expect_setequal(indiv_cmps[,unique(Sample)], names(mets)[names(mets) != "compound"])
  expect_true(indiv_cmps[,all(is.numeric(CMP))])
  expect_equal(nrow(indiv_cmps[is.na(compound)]), 0)
  expect_equal(nrow(indiv_cmps[is.na(Species)]), 0)
  expect_equal(nrow(indiv_cmps[is.na(Sample)]), 0)
  mets_melt = melt(mets, id.var = "compound", variable.name = "Sample")
  cmp_mods = fit_cmp_mods(indiv_cmps, mets_melt)
  expect_equal(nrow(cmp_mods[[1]]), length(intersect(mets[,compound], indiv_cmps[,unique(compound)])))
  expect_equal(nrow(cmp_mods[[1]])*indiv_cmps[,length(unique(Sample))], nrow(cmp_mods[[2]]))
  expect_true(cmp_mods[[2]][,all(is.numeric(Resid))])
  expect_equal(nrow(cmp_mods[[2]][is.na(Resid)]), 0)
  indiv_cmps = add_residuals(indiv_cmps, cmp_mods[[1]], cmp_mods[[2]])
  expect_equal(nrow(indiv_cmps[Species=="Residual"]), nrow(cmp_mods[[2]]))
  expect_equal(nrow(indiv_cmps[is.na(newValue)]), 0)
  expect_true(indiv_cmps[,all(is.numeric(newValue))])
  var_shares = calculate_var_shares(indiv_cmps)
  expect_known_output(var_shares, file = paste0(file_prefix, "_var_shares.rda"))
  expect_output(run_mimosa2(config_table))
}


test_that("Seq var -> AGORA species", {
  config1 = fread(test_config_file1, header = F, fill = T)
  test_results_normal(config1, file_prefix = "test_seq_agora")
})


test_that("Seq var -> Greengenes OTUs, species-rxn KEGG mods", {
  config1 = fread(test_config_file1, header = F, fill = T)
  config1[V1=="genomeChoices", V2:=get_text("source_choices")[1]]
  config1[V1=="netAdd", V2:="test_netAdd_species_rxns_KEGG.txt"]
  test_results_normal(config1, file_prefix = "test_seq_gg")
})

test_that("GG OTUs -> network", {
  config1 = fread(test_config_file1, header = F, fill = T)
  config1[V1=="database", V2:=get_text("database_choices")[2]]
  config1[V1=="genomeChoices", V2:=get_text("source_choices")[1]]
  config1[V1=="file1", V2:="test_gg.txt"]
  config1[V1=="netAdd", V2:="test_netAdd_species_rxns_KEGG.txt"]
  test_results_normal(config1, file_prefix = "test_otus_gg")
})


test_that("Greengenes OTUs -> AGORA species", {
  config1 = fread(test_config_file1, header = F, fill = T)
  config1[V1=="genomeChoices", V2:=get_text("source_choices")[2]]
  config1[V1=="database", V2:=get_text("database_choices")[2]]
  config1[V1=="file1", V2:="test_gg.txt"]
  config1[V1=="netAdd", V2:="test_netAdd_species_rxns_AGORA.txt"]
  test_results_normal(config1, file_prefix = "gg_agora_addAgora")
})

test_that("Greengenes OTUs -> AGORA species, KEGG add", {
  config1 = fread(test_config_file1, header = F, fill = T)
  config1[V1=="netAdd", V2:="test_netAdd_species_rxns_KEGG.txt"]
  test_results_normal(config1, file_prefix = "gg_agora_addKEGG")
})


#Decide what to do about this
test_that("Greengenes OTUs -> AGORA species, mixed add", {
  config1 = fread(test_config_file1, header = F, fill = T)
  config1[V1=="netAdd", V2:="test_netAdd_species_rxns_KEGG2.txt"]
  expect_error(test_results_normal(config1, file_prefix = "gg_agora_addMixed"))
})


test_that("Non-KEGG metabolites", {
  config1 = fread(test_config_file1, header = F, fill = T)
  config1[V1=="file2", V2:="test_mets_names.txt"]
  config1 = rbind(config1, data.table(V1="metType", V2=get_text("met_type_choices")[2]))
  test_results_normal(config1, file_prefix = "mets_names")
})


test_that("Species-gene modifications work, KEGG", {
  config1 = fread(test_config_file1, header = F, fill = T)
  config1[V1=="file1", V2:="test_gg.txt"]
  config1[V1=="genomeChoices", V2:=get_text("source_choices")[1]]
  config1[V1=="database", V2:=get_text("database_choices")[2]]
  config1[V1=="netAdd", V2:="test_netAdd_species_genes_KEGG.txt"]
  test_results_normal(config1, file_prefix = "test_gg_addSpecGenes")
})

test_that("Gene modifications work, KEGG", {
  config1 = fread(test_config_file1, header = F, fill = T)
  config1[V1=="file1", V2:="test_gg.txt"]
  config1[V1=="genomeChoices", V2:=get_text("source_choices")[1]]
  config1[V1=="database", V2:=get_text("database_choices")[2]]
  config1[V1=="netAdd", V2:="test_netAdd_genes_KEGG.txt"]
  test_results_normal(config1, file_prefix = "test_gg_addGenes")
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
})

test_that("Gene modifications work, AGORA", {
  config1 = fread(test_config_file1, header = F, fill = T)
  config1[V1=="file1", V2:="test_seqs.txt"]
  config1[V1=="genomeChoices", V2:=get_text("source_choices")[2]]
  config1[V1=="database", V2:=get_text("database_choices")[1]]
  config1[V1=="netAdd", V2:="test_netAdd_genes_AGORA.txt"]
  test_results_normal(config1, file_prefix = "test_agora_addGenes")
})



test_that("species-gene modifications work, AGORA", {
  config1 = fread(test_config_file1, header = F, fill = T)
  config1[V1=="file1", V2:="test_seqs.txt"]
  config1[V1=="genomeChoices", V2:=get_text("source_choices")[2]]
  config1[V1=="database", V2:=get_text("database_choices")[1]]
  config1[V1=="netAdd", V2:="test_netAdd_species_genes_AGORA.txt"]
  test_results_normal(config1, file_prefix = "test_agora_addSpecGenes")
})


test_that("Metagenome option works", {
  config1 = fread(test_config_file1, header = F, fill = T)
  config1[V1=="database", V2:=get_text("database_choices")[4]]
  config1[V1=="genomeChoices", V2:=get_text("source_choices")[1]]
  config1 = config1[V1 != "file1"]
  config1 = rbind(config1, data.table(V1="metagenome", V2 = "test_metagenome.txt"))
  config1[V1=="file2", V2:="test_metagenome_mets.txt"]
  config1[V1=="netAdd", V2:="test_netAdd_rxns_KEGG.txt"]
  test_results_normal(config1, file_prefix = "test_metagenome")
})


