#potential steps:
#get formatted KO abundances for one species
#get single species CMP scores
#do contributions analysis for all species
##Revision - aggregate to genus level 11/24/2015
##Revision 4/20/2016: single species MUSiCC
#NOTE: need to use OTUs normalized by 16S copy #
#Reorganization: making each of these a function in core_functions

#' Evaluate species contributors for a single metabolite with qPCR/species abundance data
#'
#' @import data.table
#' @param j metabolite # (usually from lapply/sapply)
#' @param cmps_sub_good CMP scores for metabolites with abundance data
#' @param all_rxns list of relevant reactions for each metabolite
#' @param subjects vector of subjects
#' @param norm_kos data.table of gene abundances
#' @param ko_net network created by generate_genomic_network
#' @param spec_abunds original species abundances
#' @param ref_kos gene abundances for each species
#' @param cor_with whether to look at the correlation of CMP scores of each species by itself with the metabolite, or of the whole community with that species removed
#' @return List of length two. Item 1 is a table of all relevant species and correlation between the scores alone and the true scores, while item 2 is a vector of taxa that were classified as potential key taxa.
#' @examples
#' sapply(1:length(metabolites), cmp_species_contributions, cmp_scores, all_rxns,
#' all_subjects, ko_abunds, ko_net, spec_abunds, ref_kos)
#' @export
cmp_species_contributions = function(j, cmps_sub_good, all_rxns, subjects, norm_kos, ko_net, spec_abunds, ref_kos, cor_with=T){
  if(!is.null(all_rxns[[j]])){
    if(!cor_with){
      compound = cmps_sub_good[j,compound]
      kos_involved = unique(all_rxns[[j]][Reversible==0,KO])
      #plan - look at correlation between old cmp score and cmp score w/out each species
      species_involved = ref_kos[KO %in% kos_involved, unique(Species)]

      species_cmp_cors = sapply(species_involved, function(x){
        spec_without = spec_abunds[,names(spec_abunds)!=x, with=F]
        vals_without = kos_from_species(spec_without, ref_kos)
        vals_without = vals_without[KO %in% names(ko_net[[1]])]
        net_mod = ko_net[[1]][,which(names(ko_net[[1]]) %in% vals_without[,KO])]
        cmp_without = data.matrix(net_mod[compound,vals_without[,KO]])%*%data.matrix(vals_without[,subjects,with=F])
        return(try(cor(as.vector(unlist(cmps_sub_good[compound,subjects,with=F])),as.vector(cmp_without), method="pearson")))
      })
      spec_cors = data.table(Species=species_involved, Cor=species_cmp_cors)
      #save KOs that have a major effect
      species_good = spec_cors[is.na(Cor)|(Cor < 0.5),Species]
      return(list(spec_cors, species_good))
      #  return(ko_cors)
    }else{
      compound = cmps_sub_good[j,compound]
      kos_involved = unique(all_rxns[[j]][Reversible==0,KO])
      species_involved = ref_kos[KO %in% kos_involved, unique(Species)]
      species_cmp_cors = sapply(species_involved, function(x){
        spec_alone = spec_abunds[Species == x]
        vals_alone = kos_from_species(spec_alone, ref_kos)
        vals_alone = vals_alone[KO %in% names(ko_net[[1]])]
        net_mod = ko_net[[1]][,which(names(ko_net[[1]]) %in% vals_alone[,KO])]
        cmp_alone = data.matrix(net_mod[compound,vals_alone[,KO]])%*%data.matrix(vals_alone[,subjects,with=F])
        return(try(cor(as.vector(unlist(cmps_sub_good[compound,subjects,with=F])),as.vector(cmp_alone), method="pearson")))
      })
      spec_cors = data.table(Species=species_involved, Cor=species_cmp_cors, compound = compound)
      spec_cors[,Pass:=ifelse(Cor > 0.5, 1, 0)]
      return(spec_cors)
    }
  } else stop("all_rxns must be a list with reaction tables for each metabolite")
}

#' Evaluate species contributors for a single metabolite with OTU and PICRUSt data, typically called by get_spec_contribs
#'
#' @import data.table
#' @param j metabolite ID # (usually from lapply/sapply)
#' @param cmps_sub_good CMP scores for metabolites with abundance data
#' @param all_rxns list of relevant reactions for each metabolite
#' @param subjects vector of subjects
#' @param norm_kos data.table of gene abundances
#' @param ko_net network created by generate_genomic_network
#' @param all_taxa vector of OTUs
#' @param single_spec_cmps single-species CMP scores calculated from get_spec_contribs function
#' @param cor_with whether to look at the correlation of CMP scores of each species by itself with the metabolite, or of the whole community with that species removed
#' @param comparison "mets" or "cmps"; whether to compare with community CMPs or metabolites themselves
#' @param met_data Optional, metabolite concentration table if comparison = "mets"
#' @return list of 2-item lists for every metabolite - 1st item is data.table of OTUs and correlations, second item is vector of OTUs with correlations > 0.5
#' @examples
#' lapply(1:length(metabolites), cmp_species_contributions_picrust, cmp_scores, all_rxns,
#' all_subjects, ko_abunds, ko_net, spec_abunds, ref_kos)
#' @export
cmp_species_contributions_picrust = function(j, cmps_sub_good, all_rxns, subjects, norm_kos, ko_net, all_taxa, single_spec_cmps, cor_with=T, comparison = "cmps", met_data = ""){
  if(!is.null(all_rxns[[j]])){
    compound = cmps_sub_good[j,compound]
    kos_involved = unique(all_rxns[[j]][Reversible==0,KO])
      if(cor_with){
        if(comparison == "cmps"){
          species_cmp_cors = sapply(1:length(all_taxa), function(x){
            if(nrow(cmps_sub_good[compound]) > 1){stop("Duplicate metabolite IDs, please fix")}
            if(!identical(single_spec_cmps[[x]], NA)){
              if(compound %in% single_spec_cmps[[x]][,compound]){
                return(try(cor(as.vector(unlist(cmps_sub_good[compound,subjects,with=F])),as.vector(unlist(single_spec_cmps[[x]][compound,subjects,with=F])), method="pearson", use = "complete.obs")))
              } else {
                return(NA)
              }
            } else { return(NA)}
          })
        } else {
          if(identical(met_data, "")) stop("Need to supply metabolite data!")
          if("KEGG" %in% names(met_data)) setnames(met_data, "KEGG", "compound")
          setkey(met_data, "compound")
          species_cmp_cors = sapply(1:length(all_taxa), function(x){
            if(!is.na(single_spec_cmps[[x]])){
              return(try(cor(as.vector(unlist(met_data[compound,subjects,with=F])),as.vector(unlist(single_spec_cmps[[x]][compound,subjects,with=F])), method="pearson")))
            } else { return(NA)}
          })
        }
        spec_cors = data.table(Species=all_taxa, Cor=species_cmp_cors, compound = compound)
        spec_cors[,Pass:=ifelse(Cor > 0.5, 1, 0)]
        return(spec_cors)
      } else {
        stop("cor_without not currently implemented")
      }
  } else stop("No reactions found for this compound")
}

#' Get variance share contributions to CMP scores for all compounds
#'
#' @import data.table
#' @param cmps_sub_good CMP scores for metabolites with abundance data
#' @param all_rxns list of relevant reactions for each metabolite
#' @param subjects vector of subjects
#' @param norm_kos data.table of gene abundances
#' @param ko_net network created by generate_genomic_network
#' @param all_taxa vector of OTUs
#' @param single_spec_cmps single-species CMP scores calculated from get_spec_contribs function
#' @param contrib_threshold Threshold for binary contribution
#' @param scaled_pos_contrib Whether to use scaled positive contributions
#' @param cov_shares Whether to calculate contribution to covariance with mets instead of variance (requires supplying metabolite data)
#' @param met_data Metabolite concentration data (if calculating contribution to covariance)
#' @return data.table of variance contributions to CMP scores
#' @examples
#' var_shares_cmps(cmps_total, all_rxns, subjects, norm_kos, ko_net, all_taxa, single_spec_cmps, contrib_threshold)
#' @export
var_shares_cmps = function(cmps_sub_good, all_rxns, subjects, norm_kos, ko_net, all_taxa, single_spec_cmps, contrib_threshold = 0.45, cov_shares = F, met_data = NULL){
  if(length(single_spec_cmps) != length(all_taxa)){
    stop("Single-species CMPs not consistent with taxa list!")
  }
  full_spec_cmps = rbindlist(lapply(1:length(single_spec_cmps), function(x){
    new_dat = single_spec_cmps[[x]][,c(subjects, "compound"),with=F]
    new_dat[,Species:=all_taxa[x]]
    return(new_dat)
  }))
  comps = full_spec_cmps[,unique(compound)]
  full_spec_cmps_wide = rbindlist(lapply(comps, function(x){
    dat = data.table(t(as.matrix(full_spec_cmps[compound==x,subjects,with=F])), Sample = subjects, compound = x)
    setnames(dat, c(full_spec_cmps[compound==x, Species], "Sample", "compound"))
    return(dat)
  }), fill = T)
  full_spec_cmps_wide = full_spec_cmps_wide[,c("compound", "Sample", names(full_spec_cmps_wide)[!names(full_spec_cmps_wide) %in% c("compound", "Sample")]), with=F]
  for(j in 1:length(all_taxa)){
    full_spec_cmps_wide[is.na(get(as.character(all_taxa[j]))), as.character(all_taxa[j]):=0]
  }

  #Calculate var shares of CMP scores
  if(cov_shares == F){
    var_shares = rbindlist(lapply(all_taxa, function(y){
      all1 = rbindlist(lapply(all_taxa, function(x){
        #      if(full_spec_cmps_wide[,!all(is.na(get(as.character(x))))] & full_spec_cmps_wide[,!all(is.na(get(as.character(y))))]){
        foo = full_spec_cmps_wide[,cov(get(as.character(x)), get(as.character(y)), use="complete.obs"), by = compound]
        #      }
        foo[,Species:=x]
        return(foo)
      }))
      all1[,Species2:=y]
    }))
    var_shares = var_shares[,sum(V1, na.rm = T),by=list(compound, Species)]
    #these are different b/c of the network generation for each species - need to make sure we don't do that. probably mimosa will move to totaling up species-level scores anyway so will be a moot point.
    #evidently not the only reason - KOs are in relative abundance space for community-level cmps. hmm.
    true_met_var = melt(cmps_sub_good, id.var = "compound")[,list(var(value), mean(value)), by = compound]
    setnames(true_met_var, c("V1", "V2"), c("TrueVar", "Mean"))
    var_shares = merge(var_shares, true_met_var, by="compound")
    var_shares[,VarShare:=V1/TrueVar]
    #if(!var_shares[,sum(VarShare), by = compound][,all(V1==1)]) stop("Inconsistency between total CMPs and sum of species contributors!")
    var_shares[,Pass:=ifelse(abs(VarShare) > contrib_threshold, 1,0)]
  } else {
    if(is.null(met_data)){
      stop("Must supply metabolite data")
    }
    met_data = melt(met_data, id.var = "compound", variable.name = "Sample")
    full_spec_cmps_wide = merge(full_spec_cmps_wide, met_data, by=c("Sample", "compound"))
    cov_shares = rbindlist(lapply(all_taxa, function(y){
      foo = full_spec_cmps_wide[,cov(get(as.character(y)), value, use = "complete.obs"), by=compound]
      foo[,Species:=y]
      return(foo)
    }))
    cmps_melt = melt(cmps_sub_good, id.var = "compound", variable.name ="Sample", value.name = "CMP")
    cmps_melt = merge(cmps_melt, met_data, by = c("compound", "Sample"))
    true_met_cov =  cmps_melt[,cov(CMP, value), by = compound]
    setnames(true_met_cov, "V1", "TrueCov")
    cov_shares = merge(cov_shares, true_met_cov, by = "compound")
    cov_shares[,VarShare:=V1/TrueCov]
    cov_shares[,Pass:=ifelse(abs(VarShare) > contrib_threshold, 1, 0)]
    var_shares = cov_shares
  }

  return(var_shares)
}


#' Get contributors for every metabolite and save
#'
#' @param contribs data.table of OTU, genes, samples and abundances
#' @return contribs with abundances normalized to copy number
#' @export
single_spec_musicc = function(contribs){
  ##Get OTU abundances from single-copy contribs for consistency
  #otus = contribs[GeneCountPerGenome == 1, list(OTU, Sample, CountContributedByOTU)]
  #setkey(otus, NULL)
  #otus = unique(otus)
  #otus[,OTURelAbund:=CountContributedByOTU/sum(CountContributedByOTU), by = Sample]
  #otus[,CountContributedByOTU:=NULL]
  #all_samps = unique(contribs[,Sample])
  contribs[,OTUAbund:=CountContributedByOTU/GeneCountPerGenome]
  otu_abunds = contribs[,list(Sample, OTU, OTUAbund)]
  #Get single line per OTU-sample
  otu_abunds = otu_abunds[,mean(OTUAbund), by=list(Sample, OTU)]
  otu_abunds[,OTURelAbund:=V1/sum(V1), by=Sample] #Ok, now this is back to original if that's what's provided
  contribs = merge(contribs, otu_abunds[,list(Sample, OTU, OTURelAbund)], by=c("OTU", "Sample"), all.x=T)
  all_otus = sort(unique(contribs[,OTU]))
  #cat(all_otus)
  all_koAbunds_byOTU = vector("list", length(all_otus))
  all_kos = contribs[,unique(Gene)]
  contribs[,singleMusicc:=OTURelAbund*GeneCountPerGenome] #Relative abundance of OTU times copy number
  return(contribs)
}

#' Get functional relative abundances from a single species' abundances over time
#'
#' @param contribs metagenome_contributions output of PICRUSt
#' @param otu_id specific OTU ID or "all"
#' @param out_prefix file path to write to
#' @return gene abundance matrix based on that species only
#' @export
make_unnormalized_single_spec = function(contribs, otu_id = "all", out_prefix){
  all_samps = unique(contribs[,Sample])
  contribs[,RelAbundSample:=CountContributedByOTU/sum(CountContributedByOTU),by=Sample] #relative abundance of this gene from this taxon in all the genes
  if(otu_id!="all"){
    contribs = contribs[OTU==otu_id]
    contribs = data.table::dcast.data.table(contribs, Gene~Sample, value.var = "RelAbundSample")
    #contribs = dcast.data.table(contribs, Gene~Sample, value.var = "ContributionPercentOfSample")
    for(j in 1:length(all_samps)){
      if(!(all_samps[j] %in% names(contribs))){
        contribs[,eval(all_samps[j]):=rep(0)]
      }
    }
    setnames(contribs, "Gene","KO")
    contribs = contribs[,c("KO",sort(all_samps)),with=F]
    write.table(contribs, file = paste0(out_prefix,otu_id,"_KOabund.txt"), quote=F, row.names=F, sep = "\t")
  }
  return(contribs)
}

#' Sum OTU data to the genus level
#'
#' Needs Greengenes 97_otu_taxonomy.txt to be in the same directory
#'
#' @param contribs metagenome_contributions output of PICRUSt
#' @param valueVar "relAbundSample" or "singleMusicc", abundance metric to use for single taxon gene abundances
#' @param taxonomy table of taxonomic classifications for OTUs, with either a "Genus" column or a Greengenes-formatted taxonomy column
#' @return contribs added to the genus level
#' @export
sum_to_genus = function(contribs, valueVar, taxonomy){ #Value var is relAbundSample or singleMusicc
  if(!"Genus" %in% names(taxonomy) & ncol(taxonomy)==2){ #Must be unprocessed file
    setnames(taxonomy, c("OTU", "Taxonomy"))
    taxonomy[,Genus:=gsub("; s__*","", Taxonomy)]
  }
  contribs = merge(contribs, taxonomy, by="OTU", all.x=T, all.y=F)
  contribs = contribs[,sum(get(valueVar)), by=list(Gene, Sample, Genus)]
  data.table::setnames(contribs, "Genus", "OTU") #now just treat genera like otus
  data.table::setnames(contribs, "V1", valueVar)
  return(contribs)
}

#' Get a list of relevant gene abundances for each species and each metabolite, optionally save
#'
#' @param contribs data.table of OTU, genes, samples and abundances
#' @param valueVar "relAbundSample" or "singleMusicc", abundance metric to use for single taxon gene abundances
#' @param out_prefix directory to write to
#' @return list of relevant gene abundances for each species and each metabolite
#' @export
contribs_by_species_list = function(contribs, valueVar, out_prefix, write_out = T){
  #Separate contribs by species and save
  all_otus = sort(unique(contribs[,OTU]))
  all_samps = unique(contribs[,as.character(Sample)])
  #cat(all_otus)
  all_koAbunds_byOTU = vector("list", length(all_otus))
  for(k in 1:length(all_otus)){
    contribs_sub = contribs[OTU==all_otus[k]]
    contribs_sub = data.table::dcast.data.table(contribs_sub, Gene~Sample, value.var = valueVar, fill = 0, fun.aggregate = sum)
    #contribs_sub = dcast.data.table(contribs_sub, Gene~Sample, value.var = "ContributionPercentOfSample")
    #Add in any samples that are all 0 for this taxon
    for(j in 1:length(all_samps)){
      if(!(all_samps[j] %in% names(contribs_sub))){
        contribs_sub[,(as.character(all_samps[j])):=rep(0)]
      }
    }
    data.table::setnames(contribs_sub, "Gene","KO")
    all_koAbunds_byOTU[[k]] = contribs_sub[,c("KO",sort(all_samps)),with=F]
  }
  if(write_out) save(all_koAbunds_byOTU, file = paste0(out_prefix,"_allKOAbundsByOTU.rda"))
  return(all_koAbunds_byOTU)
}

#' Get single-species CMP scores for every OTU and metabolite
#'
#' @param all_otus vector of OTUs or taxa
#' @param all_koAbunds_byOTU list of KO abundances for each OTU
#' @param valueVar "relAbundSample" or "singleMusicc", abundance metric to use for single taxon gene abundances
#' @param out_prefix Path of directory for loading and saving output
#' @param rxn_table Community network template
#' @param degree_filter Filter compounds connected to more than this number of KOs in the metabolic network model (default 30)
#' @param write_out Whether to save output
#' @param normalize Whether to normalize rows when making the network EMM
#' @param ko_net Optional, community network model. if provided the function will not generate one for this specific species.
#' @return list of matrices of cmp scores for each taxon
#' @export
get_all_singleSpec_cmps = function(all_otus, all_koAbunds_byOTU, valueVar, out_prefix, rxn_table, degree_filter = 30, write_out = T, normalize = T, ko_net = ""){
  cat(paste0("Getting all single-species CMP scores, using ",valueVar, " gene abundances\n"))
  if(length(all_otus) != length(all_koAbunds_byOTU)) stop("Problem! OTU lists don't match")
  cmps_alone = vector("list", length(all_otus))
  for(k in 1:length(all_otus)){
    if(nrow(rxn_table[KO %in% all_koAbunds_byOTU[[k]][,KO]]) > 0){
      if(identical(ko_net, "")){
        sub_ko_net = generate_genomic_network(all_koAbunds_byOTU[[k]][,KO], keggSource = "KeggTemplate", degree_filter = degree_filter, rxn_table = rxn_table, normalize = normalize)
        if(is.null(sub_ko_net[[1]])){
          cmps_alone[[k]] = NA
        } else {
          cmps_alone[[k]] = get_cmp_scores(sub_ko_net[[1]], all_koAbunds_byOTU[[k]])
        }
      } else {
        cmps_alone[[k]] = get_cmp_scores(ko_net[[1]], all_koAbunds_byOTU[[k]])
      }
    } else {
      cmps_alone[[k]] = NA
    }
  }
  if(write_out) save(cmps_alone, file = paste0(out_prefix, "_allCMPsAloneByOTU.rda"))
  return(cmps_alone)
}




#' Perform a series of steps related to identifying single-species contributors.
#'
#' @param contrib_file File where metagenome_contributions PICRUSt output is located
#' @param data_dir Directory where OTU data is located
#' @param results_file File where run_all_metabolites output is located
#' @param out_prefix File path where output will be written too
#' @param otu_id "All" to evaluate all detected OTUs, otherwise, a vector of OTU IDs to test for contribution
#' @param otu_file OTU abundances file
#' @param valueVar type of functional abundances to use for single species calculations - either "singleMusicc" or "RelAbundSample"
#' @param sum_to_genus Whether to sum over OTUs to the genus level (requires Greengenes taxonomy info)
#' @param comparison "mets" or "cmps"; whether to compare with community CMPs or metabolites themselves
#' @param met_data Optional, metabolite concentration table if comparison = "mets"
#' @param var_share If selected, will calculate variance share contributions instead of correlations
#' @param cov_share If selected, will calculate covariance share contributions instead of correlations (cannot select both var and cov). Must provide metabolite concentration data.
#' @param alreadyNormalized Don't re-normalize contribution table at all
#' @return table of species, metabolites, whether that species is a contributor for that metabolite and the correlation strength
#' @examples
#' read_files(gene_file, met_file)
#' @export
get_spec_contribs = function(contrib_file, data_dir, results_file, out_prefix, otu_id = "all", otu_file, valueVar = "singleMusicc", sum_to_genus, write_out = T, comparison = "cmps", met_data = "", taxonomy_file = "", var_share = F, cov_share = F, alreadyNormalized = F){
  #devtools::load_all()
  if(!valueVar %in% c("RelAbundSample", "singleMusicc")) stop("Invalid abundance metric, must be RelAbundSample or singleMusicc")
  if(!data.table::is.data.table(contrib_file)){ #Allow supplying the contrib table directly
    contribs = data.table::fread(contrib_file, stringsAsFactors = F)
  } else {
    contribs = contrib_file
  }
  contribs = contribs[CountContributedByOTU != 0]
  all_otus = sort(unique(contribs[,OTU]))
    if(valueVar == "RelAbundSample"){ #Using relative abundance out of all genes
      valueVar = "RelAbundSample"
      if(!alreadyNormalized) contribs = make_unnormalized_single_spec(contribs, otu_id, out_prefix)
      if(sum_to_genus){
        if(taxonomy_file == "") stop("Must provide taxonomy file") else {
          taxonomy = fread(taxonomy_file, fill = T) #In case of blank taxa
        }
        contribs = sum_to_genus(contribs, valueVar = valueVar, taxonomy)
        all_otus = sort(contribs[,unique(OTU)]) #Get new set of OTUs
        out_prefix = paste0(out_prefix, "_genus_") #For automatic file saving
      }
      all_koAbunds_byOTU = contribs_by_species_list(contribs, valueVar = valueVar, out_prefix, write_out)
    } else if(valueVar == "singleMusicc"){
      if(!alreadyNormalized) contribs = single_spec_musicc(contribs)
      if(sum_to_genus){
        if(taxonomy_file == "") stop("Must provide taxonomy file") else {
          taxonomy = fread(taxonomy_file)
        }
        contribs = sum_to_genus(contribs, valueVar = valueVar, taxonomy)
        all_otus = sort(contribs[,unique(OTU)]) #Get new set of OTUs
        out_prefix = paste0(out_prefix, "_genus_") #For automatic file saving
      }
      all_koAbunds_byOTU = contribs_by_species_list(contribs, valueVar = valueVar, out_prefix, write_out)
    }
    load(results_file)
    subjects_ko = names(norm_kos)[names(norm_kos)!="KO"]
    subjects = intersect(subjects_ko, unique(as.character(contribs[,Sample])))
    if(length(subjects) != length(subjects_ko)){ cat("Non-matching sample IDs between genes/metabolites and contributions, only using consensus\n")
      contribs = contribs[Sample %in% subjects]
    }
    if(length(subjects)==0){ stop("Sample IDs not consistent between contributions and genes/metabolites")}
    if(var_share == F & cov_share == F){
      #get CMP scores
      if(otu_id != "all"){
        kos_alone = all_koAbunds_byOTU[which(all_otus == otu_id)]
        kos_alone = kos_alone[,c(names(kos_alone)[names(kos_alone)!="KO"],"KO"),with=F]
        load(results_file)
        cmps_alone = get_cmp_scores(ko_net[[1]], kos_alone)
        if(write_out) write.table(cmps_alone, file = paste0(out_prefix, otu_id, "_cmps.txt"), quote=F, row.names=F, sep = "\t")
      }else{
        cmps_alone = get_all_singleSpec_cmps(all_otus, all_koAbunds_byOTU, valueVar, out_prefix = out_prefix, rxn_table = ko_net[[3]], write_out = write_out) #ko_net from results file
      }
      ##Analyze contributions
      #  if(!identical(sort(subjects), sort(unique(as.character(contribs[,Sample]))))){ stop("Samples not consistent between contributions and genes/metabolites")}
      cmps_sub_good = get_cmp_scores(ko_net[[1]], norm_kos) # Full community cmp scores
      cmp_sub_good = cmps_sub_good[node_data[,compound]]
      all_rxns = lapply(node_data[,compound], function(x){ return(ko_net[[3]][Reac==x|Prod==x])})
      all_rxns = lapply(all_rxns,get_non_rev_rxns, by_species = F)
    #  if(!all(is.na(match(c("cmps_alone", "all_otus"),ls())))){ #If all of these do not already exist
    #contribs = data.table::fread(contrib_file)
      spec_contrib=rbindlist(lapply(1:length(node_data[,compound]),cmp_species_contributions_picrust, cmps_sub_good = cmp_sub_good, all_rxns = all_rxns, subjects = subjects, norm_kos = norm_kos, ko_net = ko_net, all_taxa = all_otus, single_spec_cmps = cmps_alone, cor_with=T, comparison, met_data))
      spec_contrib = spec_contrib[!is.na(Cor)]
    } else{
      cmps_alone = get_all_singleSpec_cmps(all_otus, all_koAbunds_byOTU, valueVar, out_prefix = out_prefix, rxn_table = ko_net[[3]], write_out = write_out, ko_net = ko_net)
      cmps_sub_good = get_cmp_scores(ko_net[[1]], norm_kos) # Full community cmp scores
      cmp_sub_good = cmps_sub_good[node_data[,compound]]
      all_rxns = lapply(node_data[,compound], function(x){ return(ko_net[[3]][Reac==x|Prod==x])})
      all_rxns = lapply(all_rxns,get_non_rev_rxns, by_species = F)
      
      if(var_share == T & cov_share == F){
        spec_contrib = var_shares_cmps(cmps_sub_good = cmp_sub_good, all_rxns = all_rxns, subjects = subjects, norm_kos = norm_kos, 
                                       ko_net = ko_net, all_taxa = all_otus, single_spec_cmps = cmps_alone)
      } else if(cov_share ==T & var_share == F){
        spec_contrib = var_shares_cmps(cmps_sub_good = cmp_sub_good, all_rxns = all_rxns, subjects = subjects, norm_kos = norm_kos, 
                                       ko_net = ko_net, all_taxa = all_otus, single_spec_cmps = cmps_alone, 
                                       cov_shares = T, met_data = met_data)
      }
    }
  ## Save output
  if(write_out) write.table(spec_contrib, file = paste0(out_prefix, "_specContrib.txt"), quote=F, row.names = F, sep = "\t")
  #save(spec_contrib, file = paste0(out_dir, "/picrust_spec_contrib_",valueVar,".rda"))
  return(spec_contrib)
}

