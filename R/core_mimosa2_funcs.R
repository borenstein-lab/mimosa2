### MIMOSA linear model w/residual, then get contributions

#' Fits model scaling total CMPs to metabolite concentrations
#'
#' @import data.table
#' @param species_cmps Table of species contribution abundances
#' @param mets_melt Table of metabolite concentrations
#' @param rank_based Whether to use Rfit instead of standard linear regression
#' @param rank_type Type of robust regression to use (Rfit or mblm)
#' @param cmp_nonzero_min Filter metabolites that have fewer than this number of nonzero CMP scores
#' @return List of 2 data.tables - one with model summary results, one with model residuals
#' @examples
#' fit_cmp_mods(species_cmps, met_data)
#' @export
fit_cmp_mods = function(species_cmps, mets_melt, rank_based = F, rank_type = "Rfit", cmp_nonzero_min = 4){
    tot_cmps = species_cmps[,sum(CMP), by=list(compound, Sample)]
    #Fill in any missing 0s
    tot_cmps = melt(dcast(tot_cmps, compound~Sample, value.var = "V1", fill = 0), id.var = "compound", variable.name = "Sample")
    tot_cmps = merge(tot_cmps, mets_melt[,list(compound, Sample, value)], by = c("compound", "Sample"))
    setnames(tot_cmps, c("value.x", "value.y"), c("V1", "value"))
    if(nrow(tot_cmps) < 2) stop("Insufficient data")
    bad_mets = tot_cmps[,length(V1[V1 != 0]), by=compound][V1 < cmp_nonzero_min, compound]
    tot_cmps = tot_cmps[!compound %in% bad_mets]
    all_comps = tot_cmps[,unique(compound)]
    model_dat = data.table(compound = all_comps, Intercept = 0, Slope = 0, Rsq = 0, PVal = 0) #Make all other columns numeric
    resid_dat = data.table(expand.grid(compound = all_comps, Sample = tot_cmps[,unique(Sample)]))
    for(x in 1:length(all_comps)){
      if(!rank_based){
        scaling_mod = tryCatch(tot_cmps[compound==all_comps[x], lm(value~V1)], error=function(e){ NA})
      } else {
        if(rank_type == "mblm"){
          scaling_mod = tryCatch(tot_cmps[compound==all_comps[x], mblm::mblm(value~V1)], error=function(e){ NA})
        } else if (rank_type == "svm"){
          scaling_mod = tryCatch(tot_cmps[compound==all_comps[x], e1071::svm(value~V1, kernel = "linear")], error=function(e){ NA})
        } else {
          scaling_mod = tryCatch(tot_cmps[compound==all_comps[x], Rfit::rfit(value~V1)], error=function(e){ NA})
        }
      }
      if(!identical(scaling_mod, NA)){
        scaling_coefs = coef(scaling_mod)
        scaling_resids = resid(scaling_mod)
        model_dat[x,Intercept:=scaling_coefs[1]]
        model_dat[x,Slope:=scaling_coefs[2]]
        if(!rank_based){
          model_dat[x,Rsq:=summary(scaling_mod)$r.squared]
          model_dat[x,PVal:=anova(scaling_mod)[["Pr(>F)"]][1]]
        } else {
          if(rank_type == "mblm"){
            model_dat[x,Rsq:=0]
            model_dat[x,PVal:=tryCatch(mblm::summary.mblm(scaling_mod)$coefficients[2,4], error=function(e){ NA})]
          } else { 
            mod_rsq = tryCatch(Rfit::summary.rfit(scaling_mod, overall.test = "drop")$R2, error=function(e){ NA})
            mod_pval = tryCatch(Rfit::drop.test(scaling_mod)$p.value, error=function(e){ NA})
            model_dat[x,RsqAdj:=mod_rsq]
            model_dat[x,PVal:=mod_pval]
            model_dat[x,Rsq:=1-(scaling_mod$D1/scaling_mod$D0)]
            model_dat[x,Tauhat:=scaling_mod$tauhat]
          }
        }
        if(length(scaling_resids) != nrow(resid_dat[compound==all_comps[x]])){
          print(scaling_resids)
          print(scaling_mod)
          print(model_dat)
          print(resid_dat[compound==all_comps[x]])
          print(tot_cmps[compound==all_comps[x]])
          stop("Missing residuals")
        } 
        resid_dat[compound==all_comps[x], Resid:=scaling_resids]
      }
    }
    return(list(model_dat, resid_dat))
}

#' Calculate general summary statistics on the results of a MIMOSA analysis
#'
#' @param input_species
#' @param species 
#' @param mets 
#' @param network
#' @param indiv_cmps 
#' @param cmp_mods 
#' @param var_shares 
#' @param config_table 
#'
#' @return Table of summary statistics (1st column is description, 2nd column is value)
#' @export
#'
#' @examples
#' get_analysis_summary(input_species, species, mets, network, cmps, mod_results, contributions, config_table)
get_analysis_summary = function(input_species, species, mets, network, indiv_cmps, cmp_mods, var_shares, config_table, pval_threshold = 0.1){
  #Sample size
  if(identical(config_table[V1=="file1_type", V2], get_text("database_choices")[5])){
    sample_size = length(intersect(species[,unique(Sample)], names(mets)))
  } else {
    sample_size = length(intersect(names(species), names(mets)))
  }
  #1 get number of species mapped
  if(identical(config_table[V1=="file1_type", V2], get_text("database_choices")[5])){
    species_mapped = species[,length(unique(OTU))]
    species_orig = input_species[,length(unique(OTU))]
  } else {
    species_mapped = nrow(species)
    species_orig = nrow(input_species)
  }
  # get number of metabolites in network
  mets_orig = nrow(mets)
  if("KEGGReac" %in% names(network)){ # If bigg IDs
    num_mets_network = length(mets[compound %in% network[grepl("[e]", Reac, fixed = T),KEGGReac]|compound %in% network[grepl("[e]", Prod, fixed = T), KEGGProd], compound])
  } else {
    num_mets_network = length(mets[compound %in% network[,Reac]|compound %in% network[,Prod], compound])
  }
  #Get number of metabolites with pred
  mets_pred = indiv_cmps[,length(unique(compound))]
  #num metabolite models fit
  if(nrow(cmp_mods[[1]][!is.na(Rsq) & Rsq != 0]) > 0){
    mets_mod = cmp_mods[[1]][!is.na(Rsq) & Rsq != 0, length(unique(compound))]
    #num signif metabolites
    mets_signif = cmp_mods[[1]][PVal < pval_threshold, length(compound)]
    mets_signif_pos = cmp_mods[[1]][PVal < pval_threshold & Slope > 0, length(compound)]
  } else { 
    mets_mod = 0 
    mets_signif = 0
    mets_signif_pos = 0
  }
  #Num metabolites contributors analyzed
  if(!is.null(var_shares)){
    mets_contrib = var_shares[,length(unique(compound))]
    #Num unique contributor taxa
    taxa_contrib = var_shares[VarShare != 0 & Species != "Residual",length(unique(Species))]
  } else {
    mets_contrib = 0
    taxa_contrib = 0
  }
  summary_table = data.table(Stat = c("Sample size with complete data", "Original number of taxa (or KOs for unstratified metagenome)", "Number of mapped taxa (or KOs for unstratified metagenome)", 
                                      "Original number of metabolites", "Number of metabolites in network model", "Number of metabolites with CMP scores", 
                                      "Number of metabolites with successful model fits", paste0("Number of significant (p <", pval_threshold, ") metabolites"),
                                      paste0("Number of significant (p <", pval_threshold, ") metabolites with positive model slope"),
                                      "Number of metabolites with analyzed taxa contributors", "Number of contributing taxa"), 
                             Value = c(sample_size, species_orig, species_mapped, mets_orig, num_mets_network, mets_pred, mets_mod, mets_signif, mets_signif_pos,  
                                       mets_contrib, taxa_contrib))
  return(summary_table)
}

#' Adjusted R-squared for a rank-based model
#' 
#' @import Rfit
#' @param null_disp Null dispersion value
#' @param model_disp Model dispersion value
#' @param tauhat Estimate of tau scale parameter
#' @param df2 Degrees of freedom
#' @param df1 Model df (assumes 1)
#' @return Adjusted r-squared
#' @examples
#' j.disp.rsq.adj(null_disp, model_disp, tau, df2)
#' @export
j.disp.rsq.adj = function(null_disp, model_disp, tauhat, df2, df1 = 1){
  ## Adjusted R-squared from rfit model: 
  Fstat = ((null_disp-model_disp))/(tauhat/2)
  R2 = (df1/df2 * Fstat)/(1 + df1/df2 * Fstat)
  return(R2)
}

#' Jaeckel's dispersion from the mean (null)
#' 
#' @import Rfit
#' @param y Vector of response values
#' @param scrs Optional score values (otherwise will use Wilcoxon)
#' @return Value of Jaeckel's dispersion from the mean for this set of values
#' @examples
#' j.disp.mean(y_vals)
#' @export
j.disp.mean = function(y, scrs = NULL){
  if(is.null(scrs)) scrs = getWScores(length(y))
  e = y - mean(y, na.rm = T)
  return(drop(crossprod(e[order(e)], scrs)))
}

#' Get Wilcoxon scores for Jaeckel's dispersion calculations
#' 
#' @import Rfit
#' @param nSamps number of samples
#' @return Vector of scores (length of nSamps)
#' @examples
#' getWScores(10)
#' @export
getWScores = function(nSamps){
  return(getScores(Rfit::wscores, seq_len(nSamps)/(nSamps + 1)))
}

#' Jaeckel's dispersion from model predictions
#' 
#' @import Rfit
#' @param yhat Model predictions
#' @param y True response values
#' @param scrs Optional score values (otherwise will use Wilcoxon)
#' @return Value of Jaeckel's dispersion between predictions and the model
#' @examples
#' j.disp.fit(preds, y_vals)
#' @export
j.disp.fit = function(yhat, y, scrs = NULL){
  if(is.null(scrs)) scrs = getWScores(length(y))
  e = y - yhat
  return(drop(crossprod(e[order(e)], scrs))) ###Note this is sum(e[order(e)]*scrs)
}


#' Fits single model scaling total CMPs to concentrations of a metabolite
#'
#' @import data.table
#' @param met1 metabolite ID to fit
#' @param cmp_rxn_dat Table of species-reaction contribution abundances
#' @param mets_melt Table of metabolite concentrations
#' @param rank_based Whether to use rfit instead of lm
#' @return lm model object
#' @examples
#' fit_single_scaling_mod(met1, species_cmps, met_data)
#' @export
fit_single_scaling_mod = function(met1, cmp_rxn_dat, mets_melt, rank_based = F, rank_type = "mblm"){ #For a single metabolite
  tot_cmps = cmp_rxn_dat[compound==met1, sum(CMP), by=list(Sample, compound)]
  tot_cmps = merge(tot_cmps, mets_melt, by = c("Sample", "compound"))
  if(rank_based == F){
    scaling_mod = tryCatch(tot_cmps[compound==met1, lm(value~V1)], error=function(e){ NA})
  } else {
    if(rank_type == "mblm"){
      scaling_mod = tryCatch(tot_cmps[compound==met1, mblm::mblm(value ~ V1)], error=function(e){ NA})
    } else {
      scaling_mod = tryCatch(tot_cmps[compound==met1, Rfit::rfit(value~V1)], error=function(e){ NA})
    }
  }
  return(scaling_mod)
}

#' Greedy algorithm to test network improvements /refine rxn network
#'
#' @import data.table
#' @param network Edge list reaction table (including reversibility info)
#' @param species Table of species abundances
#' @param mets_melt Table of metabolite concentrations
#' @param manual_agora Whether to leave in agora/bigg format
#' @param rsq_factor Improvement in R squared needed to keep going
#' @param min_rsq Minimum r squared needed to accept an improvement
#' @param min_rxns Minimum number of reactions remaining to continue iterating
#' @param max_rxns Maximum number of reactions to try
#' @param rank_based Whether to use rank-based regression
#' @param species_specific Whether to adjust rxns in a species-specific or whole-community manner
#' @return A list of model fitting results as well as revised versions of the community metabolic network model and the updated CMP scores
#' @examples
#' get_best_rxn_subset(met1, species, met_data)
#' @export
fit_cmp_net_edit = function(network, species, mets_melt, manual_agora = F, rsq_factor = 1.15, min_rsq = 0.1, min_rxns = 3, max_rxns_test = 40, rank_based = F, species_specific = T){
  #Get rxn cmp scores
  species_cmps = get_species_cmp_scores(species, network, normalize = F, manual_agora = manual_agora, leave_rxns = T)[compound %in% mets_melt[,unique(compound)]]
  if(species_specific){
    species_cmps[,SpecRxn:=paste0(Species, "_", KO)]
    network[,SpecRxn:=paste0(OTU, "_", KO)]
    rxn_var = "SpecRxn"
  } else {
    rxn_var = "KO"
    species_cmps = species_cmps[,sum(CMP), by=list(Sample, KO, compound)]
    setnames(species_cmps, "V1", "CMP")
  }
  comp_order = mets_melt[,sd(value)/mean(value), by=compound][order(V1, decreasing = T), compound]
  model_dat = data.table(compound = comp_order)
  resid_dat = data.table(expand.grid(compound = comp_order, Sample = species_cmps[,unique(Sample)]))
  all_rxns_removed = c()
  #Go through in order of largest to smallest coefficient of variation

  for(j in 1:length(comp_order)){
    met1 = comp_order[j]
    print(met1)
    cmp_dat = species_cmps[compound==met1]
    met_net = network[Reac==met1|Prod==met1]

    if(nrow(cmp_dat)==0){ #If not actually in network (environmental metabolite or whatever)
      next
    }
    uniq_rxns = cmp_dat[,unique(get(rxn_var))]
    print(length(uniq_rxns))
    if(length(uniq_rxns) > max_rxns_test){ #Only test most variable
      rxn_vary = cmp_dat[,var(CMP), by=rxn_var]
      uniq_rxns = rxn_vary[order(V1, decreasing = T)][1:max_rxns_test, get(rxn_var)]
    }
    orig_scaling_mod = fit_single_scaling_mod(met1, cmp_dat, mets_melt, rank_based = rank_based)
    new_scaling_mod = orig_scaling_mod
    new_rsq = 1
    new_slope = coef(new_scaling_mod)[2]
    rsqs = rep(1, length(uniq_rxns))
    rxns_to_keep = c()
    rxns_to_remove = c()
    old_rsq = ifelse(is.na(summary(orig_scaling_mod)$r.squared), 0, summary(orig_scaling_mod)$r.squared)
    old_slope = coef(orig_scaling_mod)[2]
    if(is.na(old_slope)| length(old_slope)==0){ old_slope = 0}
    while(any(!uniq_rxns %in% rxns_to_keep) & ((length(uniq_rxns) > min_rxns & new_rsq > rsq_factor*old_rsq & new_rsq > min_rsq )|(old_slope < 0 & new_slope < 0))){
      #If new one isn't much better htan old one we're done
      orig_scaling_mod = new_scaling_mod
      old_rsq = ifelse(is.na(summary(orig_scaling_mod)$r.squared), 0, summary(orig_scaling_mod)$r.squared)
      old_slope = coef(orig_scaling_mod)[2]
      if(is.na(old_slope)| length(old_slope)==0){ old_slope = 0}

      #update full model
      #print(uniq_rxns)

      all_scaling_mods = list()
      for(k in 1:length(uniq_rxns)){
        cmp_dat1 = cmp_dat[get(rxn_var) != uniq_rxns[k]]
        all_scaling_mods[[k]] = tryCatch(fit_single_scaling_mod(met1, cmp_dat1, mets_melt, rank_based = rank_based), error=function(e){ NA})
      }
      slopes = sapply(all_scaling_mods, function(x){
        if(!identical(x, NA)){
          return(x$coefficients[2])
        } else return(0)
      }) #Rule out changes that make this negative - oh - this isn't working
      rsqs = sapply(all_scaling_mods, function(x){ 
        if(!identical(x, NA)){
        return(summary(x)$r.squared)} else {
          return(NA)
        }
        }) #Really?
      if(old_slope > 0 & length(rsqs[!is.na(rsqs)]) > 0){
        rxns_to_keep = c(rxns_to_keep, uniq_rxns[slopes <= 0]) #Keep in rxns taht would break it if removed
        rxns_to_remove = c(rxns_to_remove, sample(uniq_rxns[which(rsqs==max(rsqs, na.rm = T))], 1)) #randomly pick 1 if multiple
        print(rxns_to_remove)
        cmp_dat = cmp_dat[!SpecRxn %in% rxns_to_remove]
        new_scaling_mod = fit_single_scaling_mod(met1, cmp_dat, mets_melt, rank_based = rank_based)
        #print(summary(new_scaling_mod))
        new_rsq = summary(new_scaling_mod)$r.squared
        #Update stuff for next loop
        met_net = met_net[!SpecRxn %in% rxns_to_remove]
        uniq_rxns = uniq_rxns[!uniq_rxns %in% c(rxns_to_remove, rxns_to_keep)]
      } else {
        if(any(slopes[!is.na(slopes)] > 0)){
          rxns_to_remove = c(rxns_to_remove, sample(uniq_rxns[which(slopes > 0 & rsqs==max(rsqs[slopes > 0], na.rm = T))], 1))
          cmp_dat = cmp_dat[!SpecRxn %in% rxns_to_remove]
          new_scaling_mod = fit_single_scaling_mod(met1, cmp_dat, mets_melt, rank_based = rank_based)
          #print(summary(new_scaling_mod))
          new_rsq = summary(new_scaling_mod)$r.squared
          #Update stuff for next loop
          met_net = met_net[!SpecRxn %in% rxns_to_remove]
          uniq_rxns = uniq_rxns[!uniq_rxns %in% c(rxns_to_remove, rxns_to_keep)]
        } else { #give up?
          break
        }
      }

    }
    #Do we recalculate CMPs now? I guess so?
    if(length(rxns_to_remove) > 0){
      if(any(comp_order[(j+1):length(comp_order)] %in% c(network[SpecRxn %in% rxns_to_remove,unique(Reac)],network[SpecRxn %in% rxns_to_remove,unique(Prod)]))){
        #Recalculate only if necessary
        species_cmps = get_species_cmp_scores(species, network[!SpecRxn %in% rxns_to_remove], normalize = F, manual_agora = T, leave_rxns = T)[compound %in% mets_melt[,unique(compound)]]
        species_cmps[,SpecRxn:=paste0(Species, "_", KO)]
      }
      #Remove from full network
      network = network[!SpecRxn %in% rxns_to_remove]
      all_rxns_removed = rbind(all_rxns_removed, data.table(compound = met1, Rxn = rxns_to_remove))
    }
    print(all_rxns_removed)
    #Model data
    scaling_coefs = coef(orig_scaling_mod)
    scaling_resids = resid(orig_scaling_mod)
    model_dat[compound==met1 ,Intercept:=scaling_coefs[1]]
    model_dat[compound==met1, Slope:=scaling_coefs[2]]
    model_dat[compound==met1, Rsq:=summary(orig_scaling_mod)$r.squared]
    if(length(scaling_resids) != nrow(resid_dat[compound==met1])) stop("Missing residuals")
    resid_dat[compound==met1, Resid:=scaling_resids]
  }
  new_cmps = get_species_cmp_scores(species, network, normalize = F, manual_agora = manual_agora)

  return(list(model_dat, resid_dat, network, new_cmps))
}

#
#' Returns a melted data table of species and their scaled cmps (based on model coefficient), along with residual contributions
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
  print(model_dat)

  #all_comps = model_dat[!is.na(Rsq), compound]
  resid_dat[,Species:="Residual"]
  if("Resid" %in% names(resid_dat)) setnames(resid_dat, "Resid", "newValue")
  species_cmps[,newValue:=-1e7] #Set column type as numeric, but want to remove NAs later
  for(x in all_comps){
    if(!is.na(model_dat[compound==x,Slope])){
      species_cmps[compound==x, newValue:=as.numeric(CMP*model_dat[compound==x, Slope])]
    } else {
      species_cmps[compound==x, newValue:=NA]
    }
    species_cmps = rbind(species_cmps, resid_dat[compound==x], fill = T)
  }
  species_cmps = species_cmps[newValue != -1e7 & !is.na(newValue)]
  return(species_cmps)
}

#' Wrapper function to calculate contributions to variance, specific method depends on configuration
#'
#' @import data.table
#' @param species_contribution_table Table of species contributions (with residuals)
#' @param met_table Long-form table of metabolite concentrations (Only used for permutation-based calculation)
#' @param model_results Model fitting results from fit_cmp_mods
#' @param valueVar Column name to use for values
#' @param config_table Table of analysis settings
#' @param signif_threshold Only analyze metabolites with a model p-value less than this value (for rank-based only)
#' @return Data table of variance shares by each species in the original table for each metabolite
#' @examples
#' calculate_var_shares(species_cmps)
#' @export
calculate_var_shares = function(species_contribution_table, met_table, model_results, config_table, valueVar = "newValue", species_merge = NULL, 
                                signif_threshold = 0.2, add_summary = T){
  #Option to merge low abundance species
  if(length(species_merge) > 1){ #Merge low abundance species into "Rare/Low-abundance"
    cat(paste0("Merging ", length(species_merge), " taxa for contributor analysis\n"))
    species_contribution_table[,NewSpecies:=ifelse(Species %in% species_merge, "Rare/Low-abundance", Species)]
    species_contribution_table = species_contribution_table[,sum(CMP), by=list(compound, NewSpecies, Sample)]
    setnames(species_contribution_table, c("V1", "NewSpecies"), c("CMP", "Species"))
  }
  
  if(!"rankBased" %in% config_table[,V1]){
    #Add residuals here
    species_contribution_table = add_residuals(species_contribution_table, model_dat = model_results[[1]], 
                                               resid_dat = model_results[[2]])
    var_share_results = calculate_var_shares_linear(species_contribution_table, model_cov = T)
  } else {
    ## add option to merge low-abundance species
    var_share_results = rank_based_rsq_contribs(species_contribution_table, met_table, config_table, cmp_mods = model_results, 
                                               add_residual = T, signif_threshold = signif_threshold) #Default 5*M orderings
  }
  if(!is.null(var_share_results)){
    var_share_results[,Species:=as.character(Species)]
    var_share_results[Species != "Residual", PosVarShare:=ifelse(VarShare > 0, VarShare/sum(VarShare[VarShare > 0], na.rm=T),0), by=compound]
  } 
  #Consistnet column names
  if(!is.null(var_share_results)){
    var_share_results = merge(var_share_results, model_results[[1]], by="compound", all.x = T)
    if("Var" %in% names(var_share_results)) setnames(var_share_results, "Var", "VarDisp")
    if("NullDisp" %in% names(var_share_results)) setnames(var_share_results, "NullDisp", "VarDisp")
    #Compound-level values first, then species-level contributions
    var_share_results = var_share_results[,list(compound, Rsq, VarDisp, PVal, Slope, Intercept, Species, VarShare, PosVarShare)]
    setnames(var_share_results, "PVal", "ModelPVal")
  }
  return(var_share_results)
}


#' Calculates contributions to variance from a contribution table
#'
#' @import data.table
#' @param species_contribution_table Table of species contributions (with residuals)
#' @param valueVar Column name to use for values
#' @param model_cov Just calculate covariances with metabolite/response data
#' @param metVar If model_cov = T, variable name for metabolite concentrations
#' @return Data table of variance shares by each species in the original table for each metabolite
#' @examples
#' calculate_var_shares_linear(species_cmps)
#' @export
calculate_var_shares_linear = function(species_contribution_table, valueVar = "newValue", model_cov = F){ #generic, table of values for each speices and sample and compound
  if(!model_cov){
    spec_list = species_contribution_table[,unique(as.character(Species))]
    spec_table_wide = dcast(species_contribution_table, Sample+compound~Species, 
                            value.var = valueVar, fill = 0, fun.aggregate=sum)
    var_shares = rbindlist(lapply(spec_list, function(y){
      all1 = rbindlist(lapply(spec_list, function(x){
        foo = spec_table_wide[,cov(get(x), get(y), use="complete.obs"), by=compound]
        foo[,Species:=x]
        return(foo)
      }))
      all1[,Species2:=y]
    }))
    var_shares = var_shares[,sum(V1),by=list(compound, Species)]
  } else {
    species_contribution_table[,totVar:=sum(get(valueVar)), by=list(compound, Sample)] #Total concentrations (must include residuals)
    var_shares = species_contribution_table[,cov(get(valueVar), totVar), by=list(compound, Species)]
  }
  tot_vals = species_contribution_table[,sum(get(valueVar)), by = list(compound, Sample)]
  true_met_var = tot_vals[,list(var(V1), mean(V1)), by = compound]
  setnames(true_met_var, c("V1", "V2"), c("Var", "Mean"))
  var_shares = merge(var_shares, true_met_var, by="compound")
  var_shares[,VarShare:=V1/Var]
  var_shares[,Species:=as.character(Species)]
  return(var_shares)
}

#' Calculates contributions to variance using permutation/subset Shapley analysis for regression fits
#'
#' @import data.table
#' @param species_cmps Table of species contributions (without residuals)
#' @param mets_melt Table of metabolite concentrations
#' @param config_table Configuration table including model-fitting settings
#' @param nperm Number of permutations per decomposition
#' @param signif_threshold Will only analyze models with a p-value below this threshold
#' @param return_perm Whether to return full table of permutation results
#' @return Data table of variance shares by each species in the original table for each metabolite
#' @examples
#' run_shapley_contrib_analysis(species_cmps)
#' @export
run_shapley_contrib_analysis = function(species_cmps, mets_melt, config_table, nperm = 15000, signif_threshold = 0.01, return_perm = F){
  if("rankBased" %in% config_table[,V1]){
    rank_based = T
    if("rank_type" %in% config_table[,V1]){
      rank_type = config_table[V1=="rank_type", V2]
    } else rank_type = "rfit"
  } else rank_based = F
  
  cmp_mods = fit_cmp_mods(species_cmps, mets_melt, rank_based = rank_based, rank_type = rank_type)
  mod_fit_true = cmp_mods[[1]][!is.na(Rsq) & PVal < signif_threshold,list(compound, Rsq)]
  species_cmps = species_cmps[compound %in% mod_fit_true[,compound]]
  #Fill in 0s
  species_cmps = melt(dcast(species_cmps, compound+Sample~Species, value.var = "CMP", fill = 0), id.vars = c("compound", "Sample"), variable.name = "Species", value.name = "CMP")
    
  spec_list = species_cmps[,sort(unique(as.character(Species)))]
  nspec = length(spec_list)
  R1 = sapply(1:nperm, function(x){
    sample.int(nspec)
  }) # Matrix of permutations #fread(perm_file, header = F)[,get(paste0("V",perm_id))]
  
  allContribs = data.table()
  for(perm_id in 1:nperm){
    cumulMetVars = copy(mod_fit_true)
    setnames(cumulMetVars, "Rsq", "TrueRsq")
    spec_order = R1[,perm_id]
    for(j in 1:nspec){
      if(j < nspec){
        perm_dat = species_cmps[!Species %in% spec_list[spec_order[1:j]]]
        #Skip compounds where this species changes nothing - automatically 0
        zero_comps = species_cmps[Species == spec_list[spec_order[j]], length(CMP[CMP != 0]), by=compound][V1==0, compound]
        perm_dat = perm_dat[!compound %in% zero_comps]
        #Fill in 0s for all species even if they don't do anything for that compound
        #fit model under permutation
        mod_fit1 = suppressWarnings(suppressMessages(fit_cmp_mods(perm_dat, mets_melt, rank_based = rank_based, rank_type = rank_type)))
        cumulMetVar = mod_fit1[[1]][,list(compound, Rsq)]
        cumulMetVars = merge(cumulMetVars, cumulMetVar, by = "compound", all.x = T)
        if(j==1){
          cumulMetVars[compound %in% zero_comps, NewRsq:=TrueRsq] #If removing this species did nothing, just keep Rsq from previous species' removal
        } else {
          cumulMetVars[compound %in% zero_comps, NewRsq:=get(spec_list[spec_order[j-1]])]
        }
        cumulMetVars[is.na(Rsq), Rsq:=0]
        setnames(cumulMetVars, "Rsq", spec_list[spec_order[j]])
      } else {
        cumulMetVars[,(spec_list[spec_order[j]]):=0] 
      }
      if(j > 1){
        cumulMetVars[,paste0("Marg_", spec_list[spec_order[j]]):=get(spec_list[spec_order[j-1]]) - get(spec_list[spec_order[j]])]
      } else {
        cumulMetVars[,paste0("Marg_", spec_list[spec_order[j]]):=TrueRsq - get(spec_list[spec_order[j]])]
      }
    }
    cumulMetVars[,OrderID:=perm_id]
    cumulMetVars = cumulMetVars[,c("compound","TrueRsq", sort(names(cumulMetVars)[3:ncol(cumulMetVars)])), with=F]
    allContribs = rbind(allContribs, cumulMetVars, fill = T)
  }
  allContribs_mean = allContribs[,lapply(.SD, mean), by=compound, .SDcols = paste0("Marg_", spec_list)]
  setnames(allContribs_mean, gsub("Marg_", "", names(allContribs_mean)))
  allContribs_mean = melt(allContribs_mean, variable.name = "Species", id.var = "compound")
  allContribs_mean = merge(allContribs_mean, mod_fit_true, by = "compound", all = T)
  if(return_perm){
    return(list(PermMat = allContribs, Contribs = allContribs_mean))
  } else {
    return(allContribs_mean)
  }
}

#' Calculates contributions to variance using permutation/subset Shapley analysis for regression fits
#'
#' @import data.table
#' @param species_cmps Table of species contributions (without residuals)
#' @param mets_melt Table of metabolite concentrations
#' @param config_table Configuration table including model-fitting settings
#' @param nperm Number of permutations per decomposition
#' @param signif_threshold Will only analyze models with a p-value below this threshold
#' @param return_perm Whether to return full table of permutation results
#' @return Data table of variance shares by each species in the original table for each metabolite
#' @examples
#' run_samp_shapley_analysis(species_cmps)
#' @export
run_samp_shapley_analysis = function(species_cmps, mets_melt, config_table, nperm = 15000, signif_threshold = 0.01, return_perm = F){
  if("rankBased" %in% config_table[,V1]){
    rank_based = T
    if("rank_type" %in% config_table[,V1]){
      rank_type = config_table[V1=="rank_type", V2]
    } else rank_type = "rfit"
  } else rank_based = F
  
  cmp_mods = fit_cmp_mods(species_cmps, mets_melt, rank_based = rank_based, rank_type = rank_type)
  preds_true = species_cmps[compound %in% cmp_mods[[1]][!is.na(Rsq) & PVal < signif_threshold, compound]] ### Wait this doesn't make sense!!
  species_cmps = species_cmps[compound %in% mod_fit_true[,compound]]
  #Fill in 0s
  species_cmps = melt(dcast(species_cmps, compound+Sample~Species, value.var = "CMP", fill = 0), id.vars = c("compound", "Sample"), variable.name = "Species", value.name = "CMP")
  
  spec_list = species_cmps[,sort(unique(as.character(Species)))]
  nspec = length(spec_list)
  R1 = sapply(1:nperm, function(x){
    sample.int(nspec)
  }) # Matrix of permutations #fread(perm_file, header = F)[,get(paste0("V",perm_id))]
  
  allContribs = data.table()
  for(perm_id in 1:nperm){
    cumulMetVars = copy(mod_fit_true)
    setnames(cumulMetVars, "Rsq", "TrueRsq")
    spec_order = R1[,perm_id]
    for(j in 1:nspec){
      if(j < nspec){
        perm_dat = species_cmps[!Species %in% spec_list[spec_order[1:j]]]
        #Skip compounds where this species changes nothing - automatically 0
        zero_comps = species_cmps[Species == spec_list[spec_order[j]], length(CMP[CMP != 0]), by=compound][V1==0, compound]
        perm_dat = perm_dat[!compound %in% zero_comps]
        #Fill in 0s for all species even if they don't do anything for that compound
        #fit model under permutation
        mod_fit1 = suppressWarnings(suppressMessages(fit_cmp_mods(perm_dat, mets_melt, rank_based = rank_based, rank_type = rank_type)))
        cumulMetVar = mod_fit1[[1]][,list(compound, Rsq)]
        cumulMetVars = merge(cumulMetVars, cumulMetVar, by = "compound", all.x = T)
        cumulMetVars[compound %in% zero_comps, Rsq:=ifelse(j==1, TrueRsq, get(spec_list[spec_order[j-1]]))] #If removing this species did nothing, just keep Rsq from previous species' removal
        cumulMetVars[is.na(Rsq), Rsq:=0]
        setnames(cumulMetVars, "Rsq", spec_list[spec_order[j]])
      } else {
        cumulMetVars[,(spec_list[spec_order[j]]):=0] 
      }
      if(j > 1){
        cumulMetVars[,paste0("Marg_", spec_list[spec_order[j]]):=get(spec_list[spec_order[j-1]]) - get(spec_list[spec_order[j]])]
      } else {
        cumulMetVars[,paste0("Marg_", spec_list[spec_order[j]]):=TrueRsq - get(spec_list[spec_order[j]])]
      }
    }
    cumulMetVars[,OrderID:=perm_id]
    cumulMetVars = cumulMetVars[,c("compound","TrueRsq", sort(names(cumulMetVars)[3:ncol(cumulMetVars)])), with=F]
    allContribs = rbind(allContribs, cumulMetVars, fill = T)
  }
  allContribs_mean = allContribs[,lapply(.SD, mean), by=compound, .SDcols = paste0("Marg_", spec_list)]
  setnames(allContribs_mean, gsub("Marg_", "", names(allContribs_mean)))
  allContribs_mean = melt(allContribs_mean, variable.name = "Species", id.var = "compound")
  allContribs_mean = merge(allContribs_mean, mod_fit_true, by = "compound", all = T)
  if(return_perm){
    return(list(PermMat = allContribs, Contribs = allContribs_mean))
  } else {
    return(allContribs_mean)
  }
}


#' Calculates contributions to variance for Rfit regression models
#'
#' @import data.table
#' @param species_cmps Table of species contributions (without residuals)
#' @param mets_melt Table of metabolite concentrations
#' @param config_table Configuration table including model-fitting settings
#' @param cmp_mods Output of fit_cmp_mods
#' @param adj_rsq Whether to decompose standard disp rsq or adjusted
#' @param nperm Number of permutations per decomposition
#' @param signif_threshold Will only analyze models with a p-value below this threshold
#' @param return_perm Whether to return full table of permutation results
#' @param add_residual Whether to add the residual as an additional contributor at the end
#' @param nperm_max Maximum # of random orderings to do, regardless of the # of taxa (default 1500)
#' @return Data table of variance shares by each species in the original table for each metabolite
#' @examples
#' rank_based_rsq_contribs(species_cmps, mets_melt, config_table)
#' @export
rank_based_rsq_contribs = function(species_cmps, mets_melt, config_table, cmp_mods = NULL, adj_rsq = F, nperm_start = 0, nperm_taxa = 5, signif_threshold = 0.1, return_perm = F, 
                                  add_residual = F, nperm_max = 300){ #Shapley contribs without re-evaluating the model every time
  #Fit model if not provided
  if(is.null(cmp_mods)){
    cmp_mods = fit_cmp_mods(species_cmps, mets_melt, rank_based = T, rank_type = "rfit")
  } 
  #Two options for statistic to decompose
  if(!adj_rsq){
    mod_fit_true = cmp_mods[[1]][!is.na(Rsq) & PVal < signif_threshold,list(compound, Rsq)]
    setnames(mod_fit_true, "Rsq", "TrueRsq")
  } else {
    mod_fit_true = cmp_mods[[1]][!is.na(Rsq) & PVal < signif_threshold,list(compound, RsqAdj, Tauhat)]
    setnames(mod_fit_true, "RsqAdj", "TrueRsq")
    mod_fit_true[,DF2:=sapply(compound, function(x){
      mets_melt[compound==x & !is.na(value), length(unique(Sample))]-2
    })]
  }
  #Calculate total null dispersion for each compound
  mod_fit_true[,NullDisp:=sapply(compound, function(x){
    j.disp.mean(mets_melt[compound==x, value])
  })]
  
  cat(paste0("Analyzing contributors to ", mod_fit_true[,length(unique(compound))], " metabolites with a model p-value less than ", signif_threshold, "\n"))
  species_cmps = species_cmps[compound %in% mod_fit_true[,compound]]
  if(mod_fit_true[,length(unique(compound))]==0){
    warning("No significant metabolites, contributors will not be analyzed")
    return(NULL)
  }
  #Fill in 0s
  #Only if doing all compounds at once
  species_cmps = melt(dcast(species_cmps, compound+Species~Sample, value.var = "CMP", fill = 0), id.vars = c("compound", "Species"), variable.name = "Sample", value.name = "CMP")
  #Get predictions from full model
  species_cmps[,Species:=as.character(Species)]
  species_cmps[,Sample:=as.character(Sample)]
  species_cmps = merge(species_cmps, cmp_mods[[1]][,list(compound, Intercept, Slope)], by = "compound", all.x = T)
  species_cmps[,PredVal:=CMP*Slope] #+Intercept]
  species_cmps = merge(species_cmps, mets_melt, by = c("compound", "Sample"), all.x = T)
  
  spec_list = species_cmps[,sort(unique(Species))]
  nspec = length(spec_list)
  if(nspec == 1){
    allContribs = mod_fit_true
    allContribs[,Species:=spec_list]
    allContribs[,VarShare:=TrueRsq]
    if(add_residual){
      mod_fit_true[,ResidualContrib:=1-TrueRsq]
      resid_dat = data.table(mod_fit_true[,list(compound, TrueRsq, 1-TrueRsq)], Species = "Residual")
      setnames(resid_dat, "V3", "VarShare")
      allContribs = rbind(allContribs, resid_dat, fill = T)
    }
    return(allContribs)
  } else {
    comp_list = mod_fit_true[,compound]
    #Create all columns to begin with
    # cumulMetVarTemplate = copy(mod_fit_true)
    # cumulMetVarTemplateMarg = copy(mod_fit_true)
    # for(m in 1:length(spec_list)){
    #   cumulMetVarTemplate[,(spec_list[m]):=0]
    #   cumulMetVarTemplateMarg[,(spec_list[m]):=0]
    # }
    mod_fit_true[,NSpec:=sapply(compound, function(x){
      species_cmps[compound==x & CMP != 0, length(unique(Species))]
    })]
    mod_fit_true[,NCombo:=factorial(NSpec)]
    mod_fit_true[,NPerm:=ifelse(NCombo < nperm_taxa*NSpec, NCombo, nperm_taxa*NSpec)]
    mod_fit_true[NPerm > nperm_max, NPerm:=nperm_max]
    mod_fit_true[,AllPerm:=ifelse(NPerm==NCombo, T, F)]

    allContribs_mean_all = data.table()
    for(k in 1:length(comp_list)){
      cmps1 = species_cmps[compound == comp_list[k]]
      new_spec_list = spec_list[spec_list %in% cmps1[CMP != 0, Species]]
      cmps1 = cmps1[Species %in% new_spec_list]
      nspec = length(new_spec_list)
      cmps1a = as.matrix(dcast(cmps1, Sample~Species, value.var = "PredVal")[,2:(nspec+1)])
      cat(nspec, "taxa for compound", comp_list[k], "\n")
      null_disp_comp = mod_fit_true[compound == comp_list[k], NullDisp]
      true_rsq = mod_fit_true[compound == comp_list[k], TrueRsq]
      met_vals = species_cmps[compound==comp_list[k] & Species==new_spec_list[1]]$value #same order this way
      if(nspec==1){
        allContribs = data.table(compound = comp_list[k])
        allContribs[, (new_spec_list):= true_rsq]
      } else {
        nperm  = mod_fit_true[compound == comp_list[k], NPerm]
        all_perm = mod_fit_true[compound == comp_list[k], AllPerm]
        # if(factorial(nspec) < nperm_taxa*nspec & factorial(nspec) < nperm_max){
        #   cat("Calculating contributions across all possible subsets\n")
        #   nperm = factorial(nspec)
        #   all_perm = T
        # } else if(nperm_start == 0){ #Can specify either a hard number of random orderings or as a function of the number of taxa
        #   nperm = nperm_taxa*nspec
        #   all_perm = F
        # }
        # if(nperm > nperm_max){
        #   cat("Recommended number of permutations exceeds maximum\n")
        #   nperm = nperm_max
        #   all_perm = F
        # }
        cat("Running ", nperm, " permutations for metabolite ", comp_list[k], "\n")
        if(all_perm){
          R1 = t(gtools::permutations(nspec, nspec))
        } else {
          R1 = t(as.matrix(permute::shuffleSet(nspec, nperm)))
          # R1 = sapply(1:nperm, function(x){
          #   sample.int(nspec)
          # }) 
        }
        allContribs = matrix(nrow = nperm, ncol = length(new_spec_list))
        
        #For each random ordering, calculate marginal effect of each species
        for(perm_id in 1:nperm){
          #cmps1 = species_cmps[compound == comp_list[k]]
          cumulMetVars = rep(0, nspec)
          cumulMetVarsMarg = rep(0, nspec)
          # cumulMetVars = copy(cumulMetVarTemplate[compound == comp_list[k]])
          # cumulMetVarsMarg = copy(cumulMetVarTemplateMarg[compound == comp_list[k]])
          spec_order = R1[,perm_id]
          for(j in 1:nspec){
            if(j < nspec){
              #we are removing species spec_order[j] at each iteration
              #print(cmps1)
              #cmps1 = cmps1[Species != new_spec_list[spec_order[j]]]
              
              #Skip compounds where this species changes nothing - automatically 0
              #zero_comps = cmps1[Species == spec_list[spec_order[j]], length(CMP[CMP != 0]), by=compound][V1==0, compound]
              #perm_dat = perm_dat[!compound %in% zero_comps]
              #Fill in 0s for all species even if they don't do anything for that compound
              #fit model under permutation
              #Calculate new disp
              
              #tot_pred = cmps1[,sum(PredVal), by=Sample]$V1
              #print(tot_pred)
              if(j < nspec-1){
                tot_pred = rowSums(cmps1a[,-spec_order[1:j]])
              } else { #one column, nothing to sum
                tot_pred = cmps1a[,-spec_order[1:j]]
              }
              new_disp = j.disp.fit(tot_pred, met_vals) #, by=compound]
              #cumulMetVars = merge(cumulMetVars, new_disp, by = "compound", all.x = T)
              if(!adj_rsq){
                cumulMetVars[spec_order[j]] = (1-new_disp/null_disp_comp)
                #set(cumulMetVars, 1, 3+spec_order[j], (1-new_disp/null_disp_comp))
                #cumulMetVars[,(new_spec_list[spec_order[j]]):=(1-new_disp/null_disp_comp)]
              } else {
                cumulMetVars[,(new_spec_list[spec_order[j]]):=j.disp.rsq.adj(NullDisp, model_disp = V1, tauhat = Tauhat, df2 = DF2)]
              }
              #cumulMetVars[is.na(NewRsq), NewRsq:=0]
              #setnames(cumulMetVars, "NewRsq", spec_list[spec_order[j]])
              #cumulMetVars[,V1:=NULL]
              #zero_comp_rows = cumulMetVars[,which(compound %in% zero_comps)]
              #nonzero_rows = cumulMetVars[,which(!compound %in% zero_comps)]
              
              #set(cumulMetVars, zero_comp_rows, paste0("Marg_", spec_list[spec_order[j]]), 0)
              #cumulMetVars[compound %in% zero_comps, paste0("Marg_", spec_list[spec_order[j]]):=0] #If removing this species did nothing, just keep Rsq from previous species' removal
              if(j > 1){
                #set(cumulMetVars, zero_comp_rows, spec_list[spec_order[j]], spec_list[spec_order[j-1]])
                #cumulMetVars[compound %in% zero_comps, (spec_list[spec_order[j]]):=get(spec_list[spec_order[j-1]])] #Same as previous
                #set(cumulMetVars, nonzero_rows, paste0("Marg_", spec_list[spec_order[j]]), spec_list[spec_order[j-1]])
                #cumulMetVars[!compound %in% zero_comps,paste0("Marg_", spec_list[spec_order[j]]):=get(spec_list[spec_order[j-1]]) - get(spec_list[spec_order[j]])]
                #cumulMetVars[,paste0("Marg_", new_spec_list[spec_order[j]]):=get(new_spec_list[spec_order[j-1]]) - get(new_spec_list[spec_order[j]])]
                #set(cumulMetVarsMarg, 1, 3+spec_order[j], cumulMetVars[1, 3+spec_order[j]-1]-cumulMetVars[1, 3+spec_order[j]])
                cumulMetVarsMarg[spec_order[j]] = cumulMetVars[spec_order[j-1]] - cumulMetVars[spec_order[j]]
              } else {
                #cumulMetVars[,paste0("Marg_", new_spec_list[spec_order[j]]):=TrueRsq - get(new_spec_list[spec_order[j]])]
                #set(cumulMetVarsMarg, 1, 3+j, true_rsq-cumulMetVars[1, 3+j])
                cumulMetVarsMarg[spec_order[j]] = true_rsq - cumulMetVars[spec_order[j]]
                # cumulMetVars[compound %in% zero_comps, (spec_list[spec_order[j]]):=TrueRsq] #Same as orig
                # cumulMetVars[!compound %in% zero_comps,paste0("Marg_", spec_list[spec_order[j]]):=TrueRsq - get(spec_list[spec_order[j]])]
              }
            } else { #Final iteration
              #cumulMetVars[,(spec_list[spec_order[j]]):=0] 
              #cumulMetVars[,paste0("Marg_", new_spec_list[spec_order[j]]):=get(new_spec_list[spec_order[j-1]])]
              #set(cumulMetVarsMarg, 1, 3+j, cumulMetVars[1, 3+j-1])
              cumulMetVarsMarg[spec_order[j]] = cumulMetVars[spec_order[j-1]]
            }
          }
          #cumulMetVars[,OrderID:=perm_id]
          #cumulMetVars = cumulMetVars[,c("compound","TrueRsq", sort(names(cumulMetVars)[3:ncol(cumulMetVars)])), with=F]
          #cumulMetVarsMarg[,TrueRsq:=true_rsq]
          allContribs[perm_id,] = cumulMetVarsMarg
        }
        allContribs = data.table(compound = comp_list[k], allContribs)
        setnames(allContribs, names(allContribs)[2:ncol(allContribs)], new_spec_list)
      }
      allContribs_mean = allContribs[,lapply(.SD, mean), by=compound, .SDcols = new_spec_list]#paste0("Marg_", spec_list)]
      #setnames(allContribs_mean, gsub("Marg_", "", names(allContribs_mean)))
      allContribs_mean = melt(allContribs_mean, variable.name = "Species", id.var = "compound")
      allContribs_mean_all = rbind(allContribs_mean_all, allContribs_mean, fill = T) # we will have to set NAs to 0
    }
    allContribs_mean_all = melt(dcast(allContribs_mean_all, compound~Species, fill = 0), id.var = "compound", variable.name = "Species")
    allContribs_mean_all = merge(allContribs_mean_all, mod_fit_true, by = "compound", all = T)
    if(adj_rsq) allContribs_mean_all[,DF2:=NULL]
    
    ## Format like variance contributions
    #print(allContribs_mean_all)
    #allContribs_mean_all[,VarShare:=value]
    setnames(allContribs_mean_all, "value", "VarShare")
    if(add_residual){
      mod_fit_true[,ResidualContrib:=1-TrueRsq]
      resid_dat = data.table(mod_fit_true[,list(compound, TrueRsq, 1-TrueRsq)], Species = "Residual")
      setnames(resid_dat, "V3", "VarShare")
      allContribs_mean_all = rbind(allContribs_mean_all, resid_dat, fill = T)
    }
    if(return_perm){
      return(list(PermMat = allContribs, Contribs = allContribs_mean))
    } else {
      return(allContribs_mean_all)
    }
  } 
}



#Commonly used color scheme
getPalette = colorRampPalette(brewer.pal(12, "Paired"))

#' Plot species contributions for a single metabolite
#'
#' @import ggplot2
#' @import data.table
#' @param varShares Dataset of contributions
#' @param metabolite Compound to plot
#' @param include_zeros Whether to plot taxa that do not have any contribution
#' @param include_residual Whether to include the residual variation
#' @param merge_threshold Threshold at which low-contributing species are merged into "other"
#' @param spec_threshold Threshold number of minimum contributing species for merging into "other" to be applied
#' @param color_palette Color palette to use for taxa in plot
#' @param order_spec Whether to sort species by contribution size
#' @param taxa_max If this compound has more than this many contributing taxa, only this number will be plotted
#' @return plot of contributions
#' @examples
#' plot_contributions(varShares)
#' @export
plot_contributions = function(varShares, metabolite, metIDcol = "metID", include_zeros = F, include_residual = F, merge_threshold = 0.02, spec_threshold = 10, color_palette = NULL, order_spec = T, taxa_max = 15){
  plot_dat = varShares[get(metIDcol)==metabolite]
  if(nrow(plot_dat)==0){
    warning(paste0("No taxonomic contribution data for ", metabolite, "\n"))
    return(NULL)
  }  #Skip if nothing
  if(!include_zeros) plot_dat = plot_dat[VarShare != 0]
  if(!include_residual) plot_dat = plot_dat[Species != "Residual"]
  if(merge_threshold > 0 & plot_dat[,length(unique(Species))] > spec_threshold){ #Only merge if more species than a threshold
    merge_species = plot_dat[abs(VarShare) < merge_threshold, unique(Species)]
    if(length(merge_species) > 1){
      plot_dat[Species %in% merge_species, Species:="Other"]
      plot_dat = plot_dat[,sum(VarShare), by=Species]
      setnames(plot_dat, "V1", "VarShare")
    }
  }
  if(plot_dat[,length(unique(Species))] > taxa_max){ #If there are too many taxa, just leave some out
    top_spec = plot_dat[Species != "Other"][order(abs(VarShare), decreasing = T)][1:15, Species]
    plot_dat = plot_dat[Species %in% top_spec]
  }
  if(is.null(color_palette)) {
    color_pals = c( "black", "gray", sample(getPalette(plot_dat[,length(unique(Species))-2])))
  } else {
    color_pals = color_palette
  }
  if(order_spec){
    spec_order = plot_dat[!Species %in% c("Residual", "Other")][order(abs(VarShare), decreasing = T), as.character(unique(Species))]
    if(nrow(plot_dat[as.character(Species)=="Other"]) > 0){
      spec_order = c("Other", spec_order)
    }
    if(include_residual){
      spec_order = c("Residual", "Other", spec_order)
    }
    print(spec_order)
    plot_dat[,Species:=factor(Species, levels = spec_order)]
  }
  print(plot_dat)
  ggplot(plot_dat,  aes(y=VarShare, x = Species, fill = Species)) + geom_bar(stat = "identity") + scale_fill_manual(values = color_pals) + geom_abline(intercept = 0, slope = 0, linetype = 2) +
    theme(strip.background = element_blank(), axis.ticks.y = element_blank(), axis.text.x = element_text(size=7), axis.text.y = element_blank(),
          legend.title = element_blank(), strip.text = element_blank(), axis.title.y = element_blank(), axis.title.x = element_text(size = 10), 
          legend.text = element_text(size=10), plot.title = element_text(face = "plain", size = 9),
          panel.spacing = unit(0.15, "inches"), plot.margin = margin(0.2, 0.4, 0.3, 0.1, "inches"), legend.position = "bottom", legend.key.size = unit(0.25, "inches")) +
    ylab("Contribution to variance") + xlab("Taxon") +  coord_flip() + ggtitle(metabolite)#

}

#' Plot summary of metabolite-species contributions
#'
#' @import ggplot2
#' @import RColorBrewer
#' @import cowplot
#' @import data.table
#' @param varShares Dataset of contributions
#' @param include_zeros Whether to plot taxa that do not have any contribution
#' @param remove_resid_rescale Whether to remove the residual and rescale contributions
#' @param met_id_col Column name with metabolite IDs
#' @return plot of contributions
#' @examples
#' plot_summary_contributions(varShares)
#' @export
plot_summary_contributions = function(varShares, include_zeros = T, remove_resid_rescale = F, met_id_col = "metID", include_rsq = T){
  if(!include_zeros){
    varShares = varShares[VarShare != 0]
    bad_mets = varShares[Species != "Residual",all(is.na(VarShare)), by=compound][V1==T, compound]
    varShares = varShares[!compound %in% bad_mets]
  }
  ## Should have a function to calculate font sizes depending on the # of features
  if(met_id_col == "metID" & !"metID" %in% names(varShares)){
    varShares[,metID:=met_names(as.character(compound))]
  } else {
    varShares[,metID:=get(met_id_col)]
    varShares[is.na(metID), metID:=compound]
  }
  #Deal with IDs not in the naming database (should probably update this as well at some point)
  #Order metabolites
  varShares[,PosNeg:=factor(ifelse(Slope > 0, "Positive", "Negative"), levels = c("Positive", "Negative"))]
  resid_dat = unique(varShares[,list(metID, Rsq, PosNeg)])
  met_order = resid_dat[order(Rsq, decreasing = F), metID]
  varShares = varShares[metID %in% met_order]
  varShares[,metID:=factor(metID, levels = met_order)]
  
  
  if(include_rsq){
    resid_plot = ggplot(resid_dat, aes(x=metID, y = Rsq)) + geom_bar(stat = "identity") + scale_y_continuous(expand = c(0,0))+ theme_minimal() +
      theme(axis.text.x = element_text(angle=90, hjust=0, vjust =0.5), axis.line = element_blank(), axis.title.x = element_blank()) + ylab("Model R-squared") +
      facet_grid(~PosNeg, scales = "free_x", space = "free_x")
  }
  varShares = varShares[Species != "Residual"]
  spec_order = varShares[,length(VarShare[abs(VarShare) > 0.05]), by=Species][order(V1, decreasing = F), Species]
  varShares[,Species2:=factor(as.character(Species), levels = spec_order)]

  if(remove_resid_rescale){
    varShares[,VarShareNoResid:=VarShare/sum(VarShare), by=metID]
    plot_var = "VarShareNoResid"
    color_lab = "Scaled contribution to variance"
  } else {
    plot_var = "VarShare"
    color_lab = "Contribution to variance"
  }
  plot1 = ggplot(varShares, aes(x=metID, y = Species2)) + geom_tile(aes_string(fill = plot_var)) + theme_minimal() +
    theme(axis.text.x = element_text(angle=90, hjust=1, vjust =0.5), axis.line = element_blank(), legend.position = "bottom", legend.text = element_text(size = 7)) +
    scale_fill_gradient2(low = brewer.pal(9,"Reds")[9], mid = "white",  high = brewer.pal(9, "Blues")[9], midpoint = 0, name = color_lab) +
      ylab("Taxon")+xlab("Metabolite") + facet_grid(~PosNeg, scales = "free_x", space = "free_x")
  if(!include_rsq){
    return(plot1)
  } else {
    plot_all = plot_grid(resid_plot, plot1 + theme(axis.text.x = element_blank()), nrow = 2, align = "v", axis = "lr", rel_heights = c(1, 2.5))
    return(plot_all)
  }
}

#' Plot CMP values vs metabolite concentrations for a single metabolite
#'
#' @import ggplot2
#' @import RColorBrewer
#' @import cowplot
#' @import data.table
#' @param cmp_table Dataset of CMP scores (long format, can be species-specific or added already) for a single metabolite
#' @param met_table Dataset of metabolite concentrations (long format) for a single metabolite
#' @param mod_results Optional additional info on model comparison results
#' @param sample_col Column name for sample column (shared between the tables)
#' @return plot of contributions
#' @examples
#' cmp_met_plot(cmp_scores, mets_melt, mod_results = mods_out[[1]])
#' @export
cmp_met_plot = function(cmp_table, met_table, mod_results = NULL, sample_col = "Sample", met_title = NULL){
  if("Species" %in% names(cmp_table)){
    cmp_table = cmp_table[,sum(CMP), by=c(sample_col, "compound")]
    setnames(cmp_table, "V1", "CMP")
  }
  m_compare_mets = merge(cmp_table, met_table, by=c(sample_col, "compound"))
  setnames(m_compare_mets, "value", "Met")
  if("V1" %in% names(m_compare_mets)) setnames(m_compare_mets, "V1", "CMP")
  if(m_compare_mets[,length(unique(compound))] > 1) stop("Can only plot 1 metabolite at a time")

  #Separate out by synth only, deg only, both
  synth_deg = m_compare_mets[,list(sum(CMP > 0), sum(CMP < 0)), by=compound]
  setnames(synth_deg, c("V1", "V2"), c("SynthSamples", "DegSamples"))
  m_compare_mets = merge(m_compare_mets, synth_deg, by="compound")
  m_compare_mets[,CMPType:=ifelse(SynthSamples > 0 & DegSamples==0, "Synth", "Both")]
  m_compare_mets[,CMPType:=ifelse(DegSamples > 0 & SynthSamples==0, "Deg", CMPType)]
  
  if(!is.null(mod_results)){
    m_compare_mets = merge(m_compare_mets, mod_results, by = "compound", all = T)
    m_compare_mets[,annotResult:=round(Rsq, 3)] 
    m_compare_mets[,M2Signif:=ifelse(PVal < 0.01, T, F)]
  }
  if(is.null(met_title)|is.na(met_title)) met_title = m_compare_mets[,unique(compound)]
  
  met_plot = ggplot(m_compare_mets, aes(x=CMP, y = Met)) + theme_cowplot() + 
    theme(axis.text = element_text(size = 4), axis.title = element_text(size = 9), plot.title = element_text(face = "plain", size = 9))  + 
    ggtitle(met_title) + ylab("Metabolite")
  if(is.null(mod_results)){
    met_plot = met_plot + geom_point(alpha = 0.6, size = 1.2)
  } else {
    met_plot = met_plot + geom_point(alpha = 0.6, size = 1.2) + 
      geom_abline(slope = mod_results[,Slope], intercept = mod_results[,Intercept], linetype = 2, alpha = 0.6, color = "black") +
      annotate(geom = "text", label = m_compare_mets[1,annotResult], x = Inf, y = Inf, hjust = 1, vjust = 1, size = 2.5, fontface = "bold") 
  } #, color = ifelse(m_compare_mets[1,M2Signif], "darkred", "darkblue"
  return(met_plot)
}

#' Generate a list of plots of CMP values vs metabolite concentrations for all metabolites
#'
#' @import ggplot2
#' @import RColorBrewer
#' @import cowplot
#' @import data.table
#' @param cmp_table Dataset of CMP scores (long format, can be species-specific or added already) for a single metabolite
#' @param met_table Dataset of metabolite concentrations (long format) for a single metabolite
#' @param mod_results Optional additional info on model comparison results
#' @param sample_col Column name for sample column (shared between the tables)
#' @return plot of contributions
#' @examples
#' plot_all_cmp_mets(varShares)
#' @export
plot_all_cmp_mets = function(cmp_table, met_table, mod_results){
  tot_cmps = cmp_table[Species != "Residual", sum(CMP), by=list(compound, Sample)]
  setnames(tot_cmps, "V1", "CMP")
  comp_list = mod_results[!identical(Rsq, NA) &  Rsq != 0][order(PVal), compound]
  all_cmp_plots = vector("list", length(comp_list))
  for(i in 1:length(comp_list)){
    all_cmp_plots[[i]] = tryCatch(cmp_met_plot(tot_cmps[compound==comp_list[i]], met_table[compound==comp_list[i]], 
                         mod_results = mod_results[compound==comp_list[i]], met_title = met_names(comp_list[i])), 
                         error=function(e){ NA})
  }
  names(all_cmp_plots) = comp_list
  return(all_cmp_plots)
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
  if(configTable[V1=="file1_type", !V2 %in% get_text("database_choices")[4:5]]){ #Not metagenome-only
    if(app){
      if(!is.null(file_list[["file1"]]$datapath)){
        species = fread(file_list[["file1"]]$datapath, header = T)
      } else {
        stop("Missing species file, wrong option specified?")
      }
    } else species = fread(file_list[["file1"]], header = T)
    species = spec_table_fix(species)
    #Filter species using default abundance values
    if(configTable[V1=="specNzeroFrac", is.numeric(V2)]){
      print(configTable[V1=="specNzeroFrac"])
      species = filter_species_abunds(species, filter_type = "fracNonzero", minSampFrac = configTable[V1=="specNzeroFrac", V2])
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
    #Read metagenome file and hold it as species 
    if(app){
      species = fread(file_list[["file1"]]$datapath, header = T)
    } else {
      species = fread(file_list[["file1"]], header = T)
    }
    if(configTable[V1=="file1_type", V2==get_text("database_choices")[5]]){ #stratified
        if(!all(c("OTU", "Gene","Sample", "CountContributedByOTU") %in% names(species))){ #Assume it is picrust format or fix if not
           species = humann2_format_contributions(species, file_read = T)
        }
        #Fill in samples missing for some gene-species combos
        species = melt(dcast(species, OTU+Gene~Sample, value.var = "CountContributedByOTU", fill = 0), id.var = c("OTU", "Gene"), variable.name = "Sample", value.name = "CountContributedByOTU")
      } else { #Unstratified
        if(!"KO" %in% names(species)){
          if("function" %in% names(species)){
            setnames(species, "function", "KO")
          } else {
            stop("Invalid column name in KO data, must be either KO or function")
          }
        }
      }
    }
  #Save option to use for everything else
  humann2_metagenome = ifelse(configTable[V1=="file1_type", V2==get_text("database_choices")[5]], T, F)
  #Read metabolites
  if(app) mets = fread(file_list[["file2"]]$datapath, header = T) else mets = fread(file_list[["file2"]], header = T)
  met_nonzero_filt = ifelse(configTable[V1=="metNzeroFilter", is.numeric(V2)], configTable[V1=="metNzeroFilter", V2], 5)
  mets = met_table_fix(mets, met_nonzero_filt)
  if(humann2_metagenome == F) shared_samps = intersect(names(species), names(mets)) else {
    shared_samps = intersect(names(mets), species[,unique(Sample)])
  }
  if(length(shared_samps) < 2) stop("Sample IDs don't match between species and metabolites")
  spec_colname = ifelse(configTable[V1=="file1_type", V2==get_text("database_choices")[4]], "KO", "OTU")
  if(humann2_metagenome == F ) species = species[,c(spec_colname, shared_samps), with=F] else {
    species = species[Sample %in% shared_samps]
  }
  mets = mets[,c("compound", shared_samps), with=F]
  dat_list = list(species, mets)
  names(dat_list) = c("species", "mets")
  for(extraFile in c("netAddFile")){
    #if(!(extraFile=="metagenome" & configTable[V1=="file1_type", V2==get_text("database_choices")[4]])){ #SKip metagenome if already read in
      if(extraFile %in% names(file_list)){
        if(!is.null(file_list[[extraFile]])){
          if(file_list[[extraFile]] != F){
            if(app) dat = fread(file_list[[extraFile]]$datapath) else dat = fread(file_list[[extraFile]])
            dat_list[[length(dat_list)+1]] = dat
            names(dat_list[[length(dat_list)]]) = extraFile
          }
        }
     # }
    } 
    # else {
    #   #save species as metagenome also if that's what's happening
    #   dat_list[[extraFile]] = species
    # }
  }
  # if(configTable[V1=="file1_type", V2 %in% get_text("database_choices")[4:5]]){ # TO do either check col name or ignore this option
  #   if(configTable[V1=="file1_type", V2==get_text("database_choices")[5]]){
  #     if(!all(c("OTU", "Gene","Sample", "CountContributedByOTU") %in% names(species))){ #Assume it is picrust format or fix if not
  #       dat_list$metagenome = humann2_format_contributions(dat_list$metagenome, file_read = T)
  #     }
  #   }
  # }
  # if("metagenome" %in% names(dat_list) & configTable[V1=="file1_type", V2 != get_text("database_choices")[4]]){ #extra metagenome
  #   metagenome_samps = ifelse(configTable[V1=="metagenome_format", V1==get_text("metagenome_options")[2]], dat_list$metagenome[,unique(Sample)], names(dat_list$metagenome))
  #   if(sum(shared_samps %in% metagenome_samps) < 3){
  #     warning("Metagenome sample IDs not compatible with species and metabolites, will be ignored")
  #     dat_list$metagenome = NULL
  #   }
  # }
  return(dat_list)
}


#' Build species-specific metabolic model for MIMOSA2 analysis
#'
#' @import data.table
#' @param species Table of species/taxon abundances
#' @param config_table Data.table of input files and settings for MIMOSA
#' @param netAdd Table of netowrk information to add to the model
#' @param manual_agora Option to provide AGORA species directly (for simulation data)
#' @param degree_filt Filter currency metabolites linked to more reactions than this value
#' @return Data.table of network model of genes and reactions for each species/taxon
#' @examples
#' build_metabolic_model(config_table)
#' @export
build_metabolic_model = function(species, config_table, netAdd = NULL, manual_agora = F, degree_filt = 0){
  ### Get species to use for network if starting from seq vars
  if(!manual_agora){
    if(config_table[V1=="file1_type", V2==get_text("database_choices")[1]]){ ### Sequence variatns
      seq_list = species[,OTU]
      if(any(grepl("[0-9]+", seq_list)|grepl("[B|D-F|H-S|U-Z|b|d-f|h-s|u-z]+", seq_list))) stop("Feature IDs have non-nucleotide characters, but the sequence variant input option was selected. If the rows of your table are OTU IDs, select the option for their database source on the input page.")
      if(config_table[V1=="ref_choices", V2==get_text("source_choices")[1]]) { ## Greengenes
        seq_results = map_seqvar(seq_list, repSeqDir = paste0(config_table[V1=="data_prefix", V2], "/picrustGG/"), repSeqFile = "gg_13_8_99_db.udb", add_agora_names = F, seqID = 0.99, vsearch_path = ifelse("vsearch_path" %in% config_table[,V1], config_table[V1=="vsearch_path", V2], "vsearch")) #Run vsearch to get gg OTUs
        species[,seqID:=paste0("seq", 1:nrow(species))]
        samps = names(species)[!names(species) %in% c("OTU", "seqID")]
        new_species = merge(species, seq_results, by = "seqID", all.x=T)
        new_species = new_species[,lapply(.SD, sum), by=dbID, .SDcols = samps]
        setnames(new_species, "dbID", "OTU")
        new_species[is.na(OTU), OTU:=0] #Unassigned
        species = new_species
        mod_list = species[,OTU]
      } else if(config_table[V1=="ref_choices", V2==get_text("source_choices")[2]]){ ## AGORA
        seq_results = map_seqvar(seq_list, repSeqDir = paste0(config_table[V1=="data_prefix", V2], "/AGORA/"), repSeqFile = "agora_NCBI_16S.udb", method = "vsearch", file_prefix = "seqtemp", seqID = config_table[V1=="simThreshold", as.numeric(V2)], add_agora_names = T, vsearch_path = ifelse("vsearch_path" %in% config_table[,V1], config_table[V1=="vsearch_path", V2], "vsearch"))
        species[,seqID:=paste0("seq", 1:nrow(species))]
        samps = names(species)[!names(species) %in% c("OTU", "seqID")]
        new_species = merge(species, seq_results, by = "seqID", all.x=T)
        new_species = new_species[,lapply(.SD, sum), by=dbID, .SDcols = samps]
        new_species[is.na(dbID), dbID:="Other"]
        setnames(new_species, "dbID", "OTU")
        mod_list = seq_results[,unique(dbID)]
        species = new_species

        # if(database != get_text("database_choices")[1]){
        #   seq_results = get_agora_from_otus(species_dat[,OTU], database = database)
        # } else {
        # }
        ### Convert species abundances to AGORA species IDs

      }  else if(config_table[V1=="ref_choices", V2==get_text("source_choices")[3]]){ ## embl_gems
        seq_results = map_seqvar(seq_list, repSeqDir = paste0(config_table[V1=="data_prefix", V2], "/embl_gems/"), repSeqFile = "all_16S_seqs.udb",
                                 method = "vsearch", file_prefix = "seqtemp", seqID = config_table[V1=="simThreshold", as.numeric(V2)],
                                 add_agora_names = F, add_embl_names = T,
                                 vsearch_path = ifelse("vsearch_path" %in% config_table[,V1], config_table[V1=="vsearch_path", V2], "vsearch"))
        species[,seqID:=paste0("seq", 1:nrow(species))]
        samps = names(species)[!names(species) %in% c("OTU", "seqID")]
        new_species = merge(species, seq_results, by = "seqID", all.x=T)
        new_species = new_species[,lapply(.SD, sum), by=ModelID, .SDcols = samps]
        new_species[is.na(ModelID), ModelID:="Other"]
        setnames(new_species, "ModelID", "OTU")
        mod_list = seq_results[,unique(ModelID)]
        species = new_species
      }else stop("Model source option not found")
    } else if(config_table[V1=="file1_type", V2==get_text("database_choices")[2]]){ ## GG OTUs
      #Nothing to do
      if(config_table[V1=="ref_choices", V2==get_text("source_choices")[1]]){ #KEGG
        mod_list = species[,OTU]
      } else if(config_table[V1=="ref_choices", V2 %in% get_text("source_choices")[2:3]]){ ## AGORA or embl_gems
        if(config_table[V1=="ref_choices", V2 == get_text("source_choices")[2]]){
          cat("Mapping greengenes OTUs to AGORA...\n")
          species = otus_to_db(species, target_db = "AGORA", data_prefix = paste0(config_table[V1=="data_prefix", V2], "OTU_model_mapping_files/"))
          if(nrow(species[!is.na(OTU)]) == 0) stop("No OTUs mapped to reference models - did you select the correct reference format?")
          cat(paste0("Mapped to ", nrow(species[!is.na(OTU)]), " AGORA species"))
        } else { #embl_gems - should rename this function
          cat("Mapping greengenes OTUs to RefSeq...\n")
          species = otus_to_db(species[!is.na(OTU)], target_db = "RefSeq", data_prefix = paste0(config_table[V1=="data_prefix", V2], "OTU_model_mapping_files/"))
          if(nrow(species[!is.na(OTU)]) == 0) stop("No OTUs mapped to reference models - did you select the correct reference format?")
          cat(paste0("Mapped to ", nrow(species[!is.na(OTU)]), " RefSeq species"))
        }
        mod_list = species[!is.na(OTU),OTU]
      } else stop('Model option not implemented')
    } else if(config_table[V1=="file1_type", V2==get_text("database_choices")[3]]){ # SILVA
      if(config_table[V1=="ref_choices", V2 %in% get_text("source_choices")[2:3]]){ # AGORA
        if(config_table[V1=="ref_choices", V2 == get_text("source_choices")[2]]){
          species = otus_to_db(species, database = "SILVA", target_db = "AGORA", data_prefix = paste0(config_table[V1=="data_prefix", V2], "OTU_model_mapping_files/"))
        } else { #Embl_gems
          species = otus_to_db(species, database = "SILVA", target_db = "RefSeq", data_prefix = paste0(config_table[V1=="data_prefix", V2], "OTU_model_mapping_files/"))
        }
        mod_list = species[!is.na(OTU), OTU]
        if(nrow(species[!is.na(OTU)]) == 0) stop("No OTUs mapped to reference models - did you select the correct reference format?")
        cat(paste0("Mapped to ", nrow(species[!is.na(OTU)]), " species"))
      } else stop("This combination of taxa format and reaction source is not implemented. Please choose a different option.")
    }
    # } else {
    #   stop("This combination of taxa format and reaction source is not implemented. Please choose a different option.")
    # }
    ### Now build network from mod_list
    if(config_table[V1=="ref_choices", V2==get_text("source_choices")[1]]){ #KEGG
      if(config_table[V1=="file1_type", !V2 %in% get_text("database_choices")[c(1, 2, 4, 5)]]){
        stop("SILVA to KEGG not currently implemented")
      }
      if(config_table[V1=="file1_type", V2 %in% get_text("database_choices")[4:5]]){
        if(config_table[V1=="ref_choices", V2 != get_text("source_choices")[1]]){
          stop("This combination of taxa format and reaction source is not implemented. Please choose a different option.")
        }
        #Get network from metagenome KOs
        if(config_table[V1=="file1_type", V2==get_text("database_choices")[4]]){
          species = species[rowSums(species[,names(species) != "KO", with=F]) != 0]
          network_template = fread(paste0(config_table[V1=="kegg_prefix", V2], "/network_template.txt")) ##generate_network_template_kegg(kegg_paths[1], all_kegg = kegg_paths[2:3], write_out = F)
          network = generate_genomic_network(species[,unique(KO)], keggSource = "KeggTemplate", degree_filter = 0, rxn_table = network_template, return_mats = F)
        } else {
          #Humann2 species-KO table
          mod_list = species[,unique(OTU)]
          network_template = fread(paste0(config_table[V1=="kegg_prefix", V2], "/network_template.txt")) ##generate_network_template_kegg(kegg_paths[1], all_kegg = kegg_paths[2:3], write_out = F)
          network = rbindlist(lapply(mod_list, function(x){ #Generate network based on KOs identified in each species
            net1 = generate_genomic_network(species[OTU==x, unique(Gene)], keggSource = "KeggTemplate", degree_filter = 0, rxn_table = network_template, return_mats = F)
            net1[,OTU:=x]
            return(net1)
          }))
          network[,normalized_copy_number:=1]
        }
      } else { ##database_choices 1 or 3
        network = get_kegg_network(mod_list, net_path = paste0(config_table[V1=="data_prefix", V2], "/picrustGG/RxnNetworks/"))
      }
    } else if(config_table[V1=="ref_choices", V2==get_text("source_choices")[2]]){ #AGORA
      network = build_species_networks_w_agora(mod_list, agora_path = paste0(config_table[V1=="data_prefix", V2], "/AGORA/RxnNetworks/"))
    } else if (config_table[V1=="ref_choices", V2==get_text("source_choices")[3]]){ #embl_gems
      network = build_species_networks_w_agora(mod_list, agora_path = paste0(config_table[V1=="data_prefix", V2], "/embl_gems/RxnNetworks/"))
    }else stop('Invalid model format specified')
    if(!(config_table[V1=="file1_type", V2==get_text("database_choices")[4]])){
      #Anything other than generic KOs
      species = species[OTU %in% network[,OTU]]
    }
  } else { ##Option for simulation data
    mod_list = species[,unique(OTU)]
    network = build_species_networks_w_agora(mod_list, agora_path = paste0(config_table[V1=="data_prefix", V2], "AGORA/RxnNetworks/"))
  }
  if(!is.null(netAdd)){
    #if(config_table[V1=="netAdd", V2!=F]){
      #netAdd = fread(config_table[V1=="netAdd", V2])
    cat("Adding and/or removing genes/reactions from the network\n")
    print(netAdd)    
    network = add_to_network(network, netAdd, data_path = config_table[V1=="data_prefix", V2], kegg_path = config_table[V1=="kegg_prefix", V2])
    #}
  }
  # if(config_table[V1=="gapfill", V2 != F]){
  #   #Do stuff
  # }
  if(config_table[V1=="file1_type", V2!=get_text("database_choices")[4]]) network = network[OTU %in% species[,OTU]]
  network = filter_currency_metabolites(network, degree_filter = degree_filt)
  #Get reversible status
  if(config_table[V1=="file1_type", V2==get_text("database_choices")[4]]){
    network = get_non_rev_rxns(network, all_rxns = T, by_species = F)
  } else {
    if("Rev" %in% names(network)){ #Already have rev info (agora and embl)
      setnames(network, "Rev", "Reversible")
    } else {
      network = get_non_rev_rxns(network, all_rxns = T, by_species = T)
    }
  }
  return(list(network, species))
}

#' Updated version of getting all single-species CMP scores for every compound and taxon
#'
#' @import data.table
#' @param species_table OTU abundance table (wide format)
#' @param network Species-specific network table, product of build_network functions
#' @param normalize Whether to normalize rows when making the network EMM
#' @param relAbund Whether to use relative abundance normalization
#' @param manual_agora Whether the metabolite are already specified using AGORA/BiGG (and therefore should not be mapped back to KEGG)
#' @param humann2 Whether the species data is long-form humann2 gene-species abundances
#' @param leave_rxns Return individual abundance/direction scores for each species and rxn
#' @param kos_only Whether to call non-species-specific version of this function instead
#' @param remove_rev Whether to remove reversible reactions before calculating anything
#' @return data.table of cmp scores for each taxon and compound
#' @examples
#' get_species_cmp_scores(species_data, network)
#' @export
get_species_cmp_scores = function(species_table, network, normalize = T, relAbund = T, manual_agora = F, humann2 = F, leave_rxns = F, kos_only = F, remove_rev = T){
  if(kos_only){
    spec_cmps = get_cmp_scores_kos(species_table, network, normalize = normalize, leave_rxns = leave_rxns)
    return(spec_cmps)
  } else {
    network[is.na(stoichReac), stoichReac:=0] #solve NA problem
    network[is.na(stoichProd), stoichProd:=0]
    spec_list = species_table[,unique(OTU)]
    species_table[,OTU:=as.character(OTU)]
    network[,OTU:=as.character(OTU)]
    if(!humann2) species_table = melt(species_table, id.var = "OTU", variable.name = "Sample") else {
      if(all(c("Gene", "CountContributedByOTU") %in% names(species_table))){
        setnames(species_table, c("Gene", "CountContributedByOTU"), c("KO", "value"))
      }
    }
    species_table[,Sample:=as.character(Sample)]
    #Convert species to relative abundance if requested
    if(relAbund){
      species_table[,value:=as.double(value)]
      species_table[,value:=1000*value/sum(value), by=Sample]
      species_table[is.nan(value), value:=0] #Just in case of all-0 samples (although this is bad for other reasons)
    }
    if(length(intersect(spec_list, network[,unique(OTU)]))==0) stop("All taxa missing network information, is this the correct network model?")
    if(!all(spec_list %in% network[,unique(OTU)])) warning("Some taxa missing network information")
    if(remove_rev){ #Remove reversible reactions
      network = network[Reversible==0]
    }
    network_reacs = network[,list(OTU, KO, Reac, stoichReac, normalized_copy_number)]
    network_prods = network[,list(OTU, KO, Prod, stoichProd, normalized_copy_number)]
    network_reacs[,stoichReac:=-1*stoichReac]
    setnames(network_reacs, c("Reac", "stoichReac"), c("compound", "stoich"))
    setnames(network_prods, c("Prod", "stoichProd"), c("compound", "stoich"))
    #Remove multiple-encoded things
    setkey(network_reacs, NULL)
    setkey(network_prods, NULL)
    network_reacs = unique(network_reacs)
    network_prods = unique(network_prods)
    # network_reacs[,stoich:=stoich*normalized_copy_number] #Add in copy num/16S normalization factor
    # network_prods[,stoich:=stoich*normalized_copy_number] #Add in copy num/16S normalization factor
    if(normalize){
      network_reacs[,stoich:=as.double(stoich)]
      network_prods[,stoich:=as.double(stoich)]
      network_reacs[,stoich:=stoich/abs(sum(stoich)), by=list(OTU, compound)]
      network_prods[,stoich:=stoich/sum(stoich), by=list(OTU, compound)]
    }
    net2 = rbind(network_reacs, network_prods, fill = T)
    net2 = net2[!is.na(compound)] #remove NAs
    net2[,stoich:=stoich*normalized_copy_number] #Add in copy num/16S normalization factor
    #network[,stoichProd:=stoichProd*normalized_copy_number] #I don't think this behaves exactly the saem as multiplying later
    
    if(!humann2) spec_cmps = merge(species_table, net2, by = "OTU", allow.cartesian = T) else {
      spec_cmps = merge(species_table, net2, by = c("OTU", "KO"), allow.cartesian = T)
    }
    spec_cmps[,CMP:=value*stoich]
    spec_cmps[,SpecRxn:=paste0(OTU, "_", KO)]
    all_comps = spec_cmps[,unique(compound)]
    #Option to get abundance scores for each species and rxn
    if(!leave_rxns){
      if(length(intersect(all_comps, kegg_mapping[,KEGG])) < 2 & manual_agora==F){ #If compounds are not KEGG IDs
        #Convert AGORA IDs to KEGG IDs
        spec_cmps[,KEGG:=agora_kegg_mets(compound)]
        spec_cmps = spec_cmps[!is.na(KEGG) & grepl("[e]", compound, fixed = T)] #Going to go for just external stuff
        spec_cmps[,compound:=NULL]
        setnames(spec_cmps, "KEGG", "compound")
      }
      spec_cmps = spec_cmps[,sum(CMP), by = list(OTU, Sample, compound, value)]
      # spec_cmps = spec_cmps[,list(sum(CMP), length(unique(KO[CMP > 0])), length(unique(OTU[CMP > 0])), length(unique(SpecRxn[CMP > 0])), 
      #                             length(unique(KO[CMP < 0])), length(unique(OTU[CMP < 0])), length(unique(SpecRxn[CMP < 0]))), by=list(OTU, Sample, compound, value)]
      #spec_cmps[abs(CMP) < 10e-16 & abs(CMP) > 0, CMP:=0]
      #spec_cmps[stoich/value < 10e-]
      #Get rid of tiny values from stoich matrix errors combining synth/deg
      #setnames(spec_cmps, c("OTU", "V1", "V2", "V3", "V4", "V5", "V6", "V7"), c("Species", "CMP", "NumSynthGenes", "NumSynthSpecies", "NumSynthSpecGenes", "NumDegGenes", "NumDegSpecies", "NumDegSpecGenes"))
      setnames(spec_cmps, c("OTU", "V1"), c("Species", "CMP"))
      spec_cmps[abs(CMP)/value < 10e-15, CMP:=0]
      spec_cmps[,value:=NULL]
      #Value might be different if copy number was previously incorporated
      spec_cmps = spec_cmps[,sum(CMP), by=list(Species, Sample, compound)]
      #spec_cmps = spec_cmps[,list(sum(CMP), sum(NumSynthGenes), sum(NumSynthSpecies), sum(NumSynthSpecGenes), sum(NumDegGenes), sum(NumDegSpecies), sum(NumDegSpecGenes)), by=list(Species, Sample, compound)]
      #setnames(spec_cmps, c("V1", "V2", "V3", "V4", "V5", "V6", "V7"), c("CMP", "NumSynthGenes", "NumSynthSpecies", "NumSynthSpecGenes", "NumDegGenes", "NumDegSpecies", "NumDegSpecGenes"))
      setnames(spec_cmps, "V1", "CMP")
    } else {
      if(length(intersect(all_comps, kegg_mapping[,KEGG])) < 2 & manual_agora==F){ #If compounds are not KEGG IDs
        #Convert AGORA IDs to KEGG IDs
        spec_cmps[,KEGG:=agora_kegg_mets(compound)]
        spec_cmps = spec_cmps[!is.na(KEGG) & grepl("[e]", compound, fixed = T)] #Going to go for just external stuff
        spec_cmps[,compound:=NULL]
        spec_cmps = spec_cmps[,sum(CMP), by=list(Sample, KEGG, KO, OTU, SpecRxn)]
        setnames(spec_cmps, c("V1", "KEGG"), c("CMP", "compound"))
      }
      
      setnames(spec_cmps, "OTU", "Species")
      bad_specRxn = spec_cmps[,length(CMP[CMP != 0]), by=SpecRxn][V1==0, SpecRxn]
      spec_cmps = spec_cmps[!SpecRxn %in% bad_specRxn]
      spec_cmps[,SpecRxn:=NULL]
    }
      # if(!leave_rxns){
      #   spec_cmps = spec_cmps[,list(sum(CMP), length(unique(KO[CMP > 0])), length(unique(OTU[CMP > 0])), length(unique(SpecRxn[CMP > 0])), 
      #                               length(unique(KO[CMP < 0])), length(unique(OTU[CMP < 0])), length(unique(SpecRxn[CMP < 0]))),  by=list(Species, KEGG, Sample)] 
      #   #Add in the other stats here!
      #   #separate internal/external?
      #   setnames(spec_cmps, c("KEGG", "V1"), c("compound", "CMP"))
      # } else {
      #   setnames(spec_cmps, "KEGG", "compound")
      #   #Remove all-0 rxns
      #   bad_specRxn = spec_cmps[,length(CMP[CMP != 0]), by=SpecRxn][V1==0, SpecRxn]
      #   spec_cmps = spec_cmps[!SpecRxn %in% bad_specRxn]
      #   spec_cmps[,SpecRxn:=NULL]
      # }
  }
  return(spec_cmps)
}


#' Get summary of CMP score basis across all samples
#'
#' @param species_table OTU abundance table (wide format)
#' @param network Species-specific network table, product of build_network functions
#' @param normalize Whether to normalize rows when making the network EMM
#' @param relAbund Whether to use relative abundance normalization
#' @param manual_agora Whether the metabolite are already specified using AGORA/BiGG (and therefore should not be mapped back to KEGG)
#' @param humann2 Whether the species data is long-form humann2 gene-species abundances
#' @param kos_only Whether to call non-species-specific version of this function instead
#' @param remove_rev Whether to remove reversible reactions before calculating anything
#' @param met_subset Subset of metabolites to keep
#' @param contrib_sizes Contribution values to order by, if applicable
#' @param return_var_shares Whether to also return contributions table with summary information
#'
#' @return List of tables of stats on relevant nonzero species and reactions in the network for each compound
#' @export
#'
#' @examples
#' get_cmp_summary(species, network)
get_cmp_summary = function(species_table, network, normalize = T, relAbund = T, manual_agora = F, humann2 = F, kos_only = F, remove_rev = T, met_subset = NULL, contrib_sizes = NULL, return_var_shares = T){
  spec_rxn_cmps = get_species_cmp_scores(species_table, network, normalize = normalize, relAbund = relAbund, manual_agora = manual_agora, humann2 = humann2, leave_rxns = T, kos_only = kos_only, remove_rev = remove_rev)
  if(!kos_only) spec_rxn_cmps[,SpecRxn:=paste0(Species, "/", KO)]
  if(!is.null(met_subset)){
    spec_rxn_cmps = spec_rxn_cmps[compound %in% met_subset]
  }
  if(!is.null(contrib_sizes)){
    spec_rxn_cmps = merge(spec_rxn_cmps, contrib_sizes[,list(compound, Species, VarShare)], by = c("compound", "Species"), all = T)
  } else {
    spec_rxn_cmps[,VarShare:=NA] #use same ordering code either way
  }
  if(!kos_only){
    spec_rxn_summary = spec_rxn_cmps[Species != "Residual" & !is.na(SpecRxn)][order(abs(VarShare), decreasing = T),list(length(unique(KO[CMP > 0])), length(unique(Species[CMP > 0])), length(unique(SpecRxn[CMP > 0])), 
                                           length(unique(KO[CMP < 0])), length(unique(Species[CMP < 0])), length(unique(SpecRxn[CMP < 0])), 
                                           paste0(unique(KO[CMP > 0]), collapse = " "), paste0(unique(Species[CMP > 0]), collapse = " "), 
                                           paste0(unique(KO[CMP < 0]), collapse = " "), paste0(unique(Species[CMP < 0]), collapse = " "),
                                           paste0(unique(SpecRxn[CMP > 0]), collapse = " "),
                                           paste0(unique(SpecRxn[CMP < 0]), collapse = " ")
                                           ), by = compound]
    setnames(spec_rxn_summary, c("compound", "NumSynthGenes", "NumSynthSpecies", "NumSynthSpecGenes", "NumDegGenes", "NumDegSpecies", "NumDegSpecGenes", "SynthGenes", "SynthSpec", "DegGenes", "DegSpec", "SynthSpecGenes", "DegSpecGenes"))
    spec_rxn_summary[,TopSynthSpecGenes:=sapply(1:nrow(spec_rxn_summary), function(x){
      if(spec_rxn_summary[x,NumSynthSpecGenes==0]){
        return("")
      } else {
        foo = unlist(strsplit(spec_rxn_summary[x,as.character(SynthSpecGenes)], split = " "))
        if(spec_rxn_summary[x,NumSynthSpecGenes] > 5){
          return(paste0(foo[1:5], collapse = " "))
        } else {
          return(paste0(foo[1:spec_rxn_summary[x,NumSynthSpecGenes]], collapse = " "))
        }
      }
      })]
    spec_rxn_summary[,TopDegSpecGenes:=sapply(1:nrow(spec_rxn_summary), function(x){
      if(spec_rxn_summary[x,NumDegSpecGenes==0]){
        return("")
      } else {
        foo = unlist(strsplit(spec_rxn_summary[x,DegSpecGenes], split = " "))
        if(spec_rxn_summary[x,NumDegSpecGenes] > 5){
          return(paste0(foo[1:5], collapse = " "))
        } else {
          return(paste0(foo[1:spec_rxn_summary[x,NumDegSpecGenes]], collapse = " "))
        }
      }
    })]
    if(return_var_shares & !is.null(contrib_sizes)){
      species_level_summary = spec_rxn_cmps[Species != "Residual" & !is.na(SpecRxn), list(length(unique(KO[CMP > 0])),  
                                                                                          length(unique(KO[CMP < 0])),  
                                                                                          paste0(unique(KO[CMP > 0]), collapse = " "), 
                                                                                          paste0(unique(KO[CMP < 0]), collapse = " ")), by = list(compound, Species)]
      setnames(species_level_summary, c("compound", "Species", "NumSynthGenes", "NumDegGenes","SynthGenes",  "DegGenes"))
      return(list(CompLevelSummary = spec_rxn_summary, SpeciesLevelSummary = species_level_summary))
    }
  } else{
    spec_rxn_summary = spec_rxn_cmps[Species != "Residual" & !is.na(KO),list(length(unique(KO[CMP > 0])),  
                                           length(unique(KO[CMP < 0])),  
                                           paste0(unique(KO[CMP > 0]), collapse = " "), 
                                           paste0(unique(KO[CMP < 0]), collapse = " ")), by = compound]
    setnames(spec_rxn_summary, c("compound", "NumSynthGenes", "NumDegGenes","SynthGenes",  "DegGenes"))
    spec_rxn_summary[,TopSynthGenes:=sapply(1:nrow(spec_rxn_summary), function(x){
      if(spec_rxn_summary[x,NumSynthGenes==0]){
        return("")
      } else {
        foo = unlist(strsplit(spec_rxn_summary[x,SynthGenes], split = " "))
        if(spec_rxn_summary[x,NumSynthGenes] > 5){
          return(paste0(foo[1:5], collapse = " "))
        } else {
          return(paste0(foo[1:spec_rxn_summary[x,NumSynthGenes]], collapse = " "))
        }
      }
    })]
    spec_rxn_summary[,TopDegGenes:=sapply(1:nrow(spec_rxn_summary), function(x){
      if(spec_rxn_summary[x,NumDegGenes==0]){
        return("")
      } else {
      foo = unlist(strsplit(spec_rxn_summary[x,DegGenes], split = " "))
      if(spec_rxn_summary[x,NumDegGenes] > 5){
        return(paste0(foo[1:5], collapse = " "))
      } else {
        return(paste0(foo[1:spec_rxn_summary[x,NumDegGenes]], collapse = " "))
      }
      }
    })]
    spec_rxn_summary[,TopSynthSpecGenes:=TopSynthGenes]
    spec_rxn_summary[,TopDegSpecGenes:=TopDegGenes]
    ##Add other columns as NA
    for(colname in c("NumSynthSpecies", "NumSynthSpecGenes",  "NumDegSpecies", "NumDegSpecGenes", "SynthSpec", "DegSpec", "SynthSpecGenes", "DegSpecGenes")){
      spec_rxn_summary[,eval(colname):=""]
    }
  }
  return(list(CompLevelSummary = spec_rxn_summary))
  
}

#' Updated version of getting all sample-level CMP scores from a KO abundance table
#'
#' @import data.table
#' @param ko_table KO abundance table (wide format)
#' @param network Species-specific network table, product of build_network functions
#' @param normalize Whether to normalize rows when making the network EMM
#' @param remove_rev Whether to remove rev rxns before calculating/normalizing
#' @param leave_rxns Whether to leave individual KO contributions or aggregate to metabolite level
#' @return data.table of cmp scores for each taxon and compound
#' @examples
#' get_cmp_scores_kos(ko_data, network)
#' @export
get_cmp_scores_kos = function(ko_table, network, normalize = T, relAbund = T, remove_rev = T, leave_rxns = F){
  network[is.na(stoichReac), stoichReac:=0] #solve NA problem
  network[is.na(stoichProd), stoichProd:=0]
  #network[,stoichReac:=stoichReac*normalized_copy_number] #Add in copy num/16S normalization factor
  #network[,stoichProd:=stoichProd*normalized_copy_number]
  ko_table_melt = melt(ko_table, id.var = "KO", variable.name = "Sample")
  ko_table_melt[,Sample:=as.character(Sample)]
  if(relAbund){
    ko_table_melt[,value:=as.double(value)]
    ko_table_melt[,value:=value/sum(value)*1000, by=Sample]
  }
  if(remove_rev){ #Remove reversible reactions
    #network = get_non_rev_rxns(network, all_rxns = T, by_species = F)
    network = network[Reversible==0]
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
  if(!leave_rxns){
    spec_cmps = spec_cmps[,sum(CMP), by = list(Sample, compound)]
    setnames(spec_cmps, "V1", "CMP")
    # spec_cmps = spec_cmps[,list(sum(CMP), length(unique(KO[CMP > 0])), length(unique(KO[CMP < 0]))), by=list(Sample, compound)]
    # setnames(spec_cmps, c("V1", "V2", "V3"), c("CMP", "NumSynthGenes", "NumDegGenes"))
    # spec_cmps[,NumSynthSpecies:=""]
    # spec_cmps[,NumDegSpecies:=""]
    # spec_cmps[,NumSynthSpecGenes:=""]
    # spec_cmps[,NumDegSpecGenes:=""]
  }
  #all_comps = spec_cmps[,unique(compound)]
  # if(length(intersect(all_comps, kegg_mapping[,KEGG])) < 2){ #If compounds are not KEGG IDs - this will not happen with KOs
  #   #Convert AGORA IDs to KEGG IDs
  #   spec_cmps[,KEGG:=agora_kegg_mets(compound)]
  #   spec_cmps = spec_cmps[!is.na(KEGG)]
  #   spec_cmps = spec_cmps[,list(sum(CMP), length(unique(KO[CMP > 0])), length(unique(KO[CMP < 0]))), by=list( KEGG, Sample)]
  #   setnames(spec_cmps, c("KEGG", "V1", "V2", "V3"), c("compound", "CMP", "NumSynthGenes", "NumDegGenes"))
  # }
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
#' @param kegg_path Path to KEGG database
#' @param data_path Path to reference data with AGORA-KEGG mappings
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
      # if("OTU" %in% table_type){
      #   copyNums = fread(ifelse(target_format=="Cobra", paste0(data_path, "blastDB/agora_NCBItax_processed_nodups.txt"), paste0("gunzip -c ", data_path, "picrustGenomeData/16S_13_5_precalculated.tab.gz")))
      #   if(target_format == "Cobra"){
      #     copyNums = unique(copyNums[,list(AGORA_ID, CopyNum)])
      #     specID = "AGORA_ID"
      #   } else {
      #     setnames(copyNums, c("OTU", "CopyNum"))
      #     specID = "OTU"
      #   }
      #   print(copyNums)
      #   addTable = merge(addTable, copyNums, by.x = "OTU", by.y = specID, all.x = T, all.y=F)
      #   addTable[is.na(CopyNum), CopyNum:=1]
      #   addTable[,normalized_copy_number:=1/CopyNum]
      #   addTable[,CopyNum:=NULL]
      # } else {
        addTable[,normalized_copy_number:=1]
      #}
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
    req_params = c("file1_type", "ref_choices")
  } else {
    req_params = c("file1", "file2", "file1_type", "ref_choices", "data_prefix")
    # if(config_table[V1=="file1_type", V2==get_text("database_choices")[4]]){
    #   req_params[req_params=="file1"] = "metagenome"
    # }
  }
  if(any(!req_params %in% config_table[,V1])){
    missing_param = req_params[!req_params %in% config_table[,V1]]
    stop(paste0("Required parameters missing from configuration file: ", missing_param, "\n"))
  } 
  all_params = c(req_params,  "metType", "netAdd", "simThreshold", "kegg_prefix", "vsearch_path", "compare_only") #Move to package sysdata?  "metagenome", "metagenome_format",
  config_table[V2=="", V2:=FALSE]
  if(!"kegg_prefix" %in% config_table[,V1]){
    config_table = rbind(config_table, data.table(V1 = "kegg_prefix", V2 = paste0(config_table[V1=="data_prefix", V2], "/KEGGfiles/")))
  }
  if(length(all_params[!all_params %in% config_table[,V1]]) > 0){
    config_table = rbind(config_table, data.table(V1 = all_params[!all_params %in% config_table[,V1]], V2 = FALSE))
  }
  #if non-species metagenome is provided, set compare_only flag
  if(config_table[V1=="file1_type", V2==get_text("database_choices")[4]]){
    config_table[V1=="compare_only", V2:="TRUE"]    
  }
  return(config_table)
}

#' Run a MIMOSA 2 analysis
#'
#' @import data.table
#' @param config_table Data.table of input files and settings for MIMOSA analysis OR path to such a table
#' @param species Optionally provide already-read-in species data
#' @param mets Optionally provide already-read-in metabolite data
#' @param make_plots Whether to generate plots for each metabolite and extended output (as provided via the web server)
#' @return Scaling model and variance contribution results
#' @examples
#' run_mimosa2(config_table, species, mets)
#' @export
run_mimosa2 = function(config_table, species = "", mets = "", make_plots = F, save_plots = F){
  #process arguments
  #Read config table if it is filename
    if(typeof(config_table)=="character"){
      config_table = fread(config_table, header = F)
    }
    if(!identical(species, "") & !identical(mets, "")){
      config_table = check_config_table(config_table, app = T)
      data_inputs = list(species = species, mets = mets)
    } else {
      config_table = check_config_table(config_table, app = F)
      file_list = as.list(config_table[V1 %in% c("file1", "file2", "netAddFile"), V2])
      names(file_list) = config_table[V1 %in% c("file1", "file2", "netAddFile"), V1]
      data_inputs = read_mimosa2_files(file_list, config_table, app = F)
      species = data_inputs$species
      mets = data_inputs$mets
    }
    if(!"manualAGORA" %in% config_table[,V1]){
      network_results = build_metabolic_model(species, config_table, degree_filt = 0, netAdd = data_inputs$netAdd)
      network = network_results[[1]]
      species = network_results[[2]]
    } else {
      network_results = build_metabolic_model(species, config_table, manual_agora = T, degree_filt = 0, netAdd = data_inputs$netAdd)
      network = network_results[[1]]
      species = network_results[[2]]
    }
    # if(config_table[V1=="file1_type", !V2 %in% get_text("database_choices")[4:5]]){
    #   #If we are doing a comparison of the species network and the metagenome network
    #   #Metagenome data
    #   #Implement doing stuff with this later
    #   metagenome_network = build_metabolic_model(data_inputs$metagenome, config_table, netAdd = data_inputs$netAdd)
    #   # species2 = metagenome_data[[1]]
    #   # metagenome_network = metagenome_data[[2]]
    #   #Metagenome data
    # }
    
    if(config_table[V1=="metType", V2 ==get_text("met_type_choices")[2]]){ #Assume it is KEGG unless otherwise specified
      mets = map_to_kegg(mets)
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
    if(config_table[V1=="file1_type", V2==get_text("database_choices")[4]]){
      no_spec_param = T
      humann2_param = F
      rel_abund_param = T
    } else if(config_table[V1=="file1_type", V2==get_text("database_choices")[5]]){
      no_spec_param = F
      humann2_param = T
      rel_abund_param = F
      cat("Humann2 format\n")
    } else {
      no_spec_param = F
      humann2_param = F
      rel_abund_param = T
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
    } else if("logTransform" %in% config_table[,V1]){
      if(config_table[V1=="logTransform", V2==T]){
        met_transform = "logplus"
        cat(paste0("Will transform metabolite values, transform is logplus"))
      } else met_transform = ""
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
      indiv_cmps = get_species_cmp_scores(species, network, normalize = !rxn_param, leave_rxns = rxn_param, manual_agora = agora_param, kos_only = no_spec_param, humann2 = humann2_param, relAbund = rel_abund_param)
      if(score_transform != ""){
        indiv_cmps = transform_cmps(indiv_cmps, score_transform)
      }
      indiv_cmps = indiv_cmps[compound %in% mets[,compound]]
      cmp_mods = fit_cmp_mods(indiv_cmps, mets_melt, rank_based = rank_based, rank_type = rank_type)
    }
    #indiv_cmps = add_residuals(indiv_cmps, cmp_mods[[1]], cmp_mods[[2]])
    if(!compare_only & !no_spec_param){ #Option to skip contributions
      if(!humann2_param){
        spec_dat = melt(species, id.var = "OTU", variable.name = "Sample")[,list(value/sum(value), OTU), by=Sample] #convert to relative abundance
        bad_spec = spec_dat[,list(length(V1[V1 != 0])/length(V1), max(V1)), by=OTU]
        bad_spec = bad_spec[V1 < 0.2 & V2 < 0.1, OTU] #Never higher than 10% and absent in at least 90% of samples
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
    #Rxns, taxa summary
    cmp_summary = get_cmp_summary(species, network, normalize = !rxn_param, manual_agora = F, kos_only = no_spec_param, humann2 = humann2_param, 
                                  met_subset = cmp_mods[[1]][!is.na(Rsq) & Rsq != 0,compound], contrib_sizes = var_shares)
    cmp_mods[[1]] = merge(cmp_mods[[1]], cmp_summary$CompLevelSummary, by = "compound", all.x = T)
    #Add species/rxn info
    
    if(length(cmp_summary) > 1) var_shares = merge(var_shares, cmp_summary$SpeciesLevelSummary, by = c("compound", "Species"), all.x = T)
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
          if(is.na(met_names(x))){
            met_id = x
          } else { met_id = met_names(x)}
          plot_contributions(var_shares, met_id, metIDcol = "MetaboliteName", color_palette = contrib_color_palette, include_residual = F, merge_threshold = 0.01) + theme(plot.background = element_blank())
        })
        #Contribution Legend
        leg_dat = data.table(V1 = factor(names(contrib_color_palette), levels = c(all_contrib_taxa, "Other"))) #, "Residual"
        setnames(leg_dat, "V1", "Contributing Taxa")
        legend_plot = ggplot(leg_dat, aes(fill = `Contributing Taxa`, x=`Contributing Taxa`)) + geom_bar() + scale_fill_manual(values = contrib_color_palette, name = "Contributing Taxa")# + theme(legend.text = element_text(size = 10))
        contrib_legend = tryCatch(get_legend(legend_plot), error = function(){ return(NULL)}) 
      } else {
        met_contrib_plots = NULL
        contrib_legend = NULL
      }
    }
    if(config_table[V1=="ref_choices", V2 != get_text("source_choices")[1]]){
      network[,KEGGReac:=agora_kegg_mets(Reac)]
      network[,KEGGProd:=agora_kegg_mets(Prod)]
    } 
    analysis_summary = get_analysis_summary(input_species = data_inputs[[1]], species = species, mets = mets, network = network, indiv_cmps = indiv_cmps, cmp_mods = cmp_mods, var_shares = var_shares, config_table = config_table, pval_threshold = signifThreshold)
    if(save_plots & make_plots){
      #Save plots
      if(!dir.exists("mimosa2plots")){
        dir.create(path = "mimosa2plots", showWarnings = T)
      }
      for(i in 1:length(CMP_plots)){
        print(paste0("mimosa2plots/", names(CMP_plots)[i], ".png"))
        if(!identical(CMP_plots[[i]], NA)){
          save_plot(CMP_plots[[i]], file = paste0("mimosa2plots/", names(CMP_plots)[i], ".png"), base_width = 2, base_height = 2)
        } 
      }
      if(!config_table[V1 == "compare_only", V2==T]){
        for(i in 1:length(met_contrib_plots)){
          print(comp_list[i])
          if(!is.null(met_contrib_plots[[i]])){
            save_plot(met_contrib_plots[[i]] + guides(fill = F), file = paste0("mimosa2plots/", comp_list[i], "_contribs.png"), base_width = 2, base_height = 2)
          }
        }
        if(!is.null(contrib_legend)) save_plot(contrib_legend, file = paste0("mimosa2plots/contribLegend.png"), dpi=120, base_width = 5, base_height = 4)
      } 
    }
    if(make_plots){
      return(list(varShares = var_shares, modelData = cmp_mods[[1]], 
                  networkData = network, newSpecies = species, CMPScores = indiv_cmps[CMP != 0], 
                  analysisSummary = analysis_summary, configs = config_table, CMPplots = CMP_plots, 
                  metContribPlots = met_contrib_plots, plotLegend = contrib_legend))
    } else {
      return(list(varShares = var_shares, modelData = cmp_mods[[1]], networkData = network, newSpecies = species, 
                  CMPScores = indiv_cmps[CMP != 0], analysisSummary = analysis_summary, configs = config_table))
    }
    

}

#' Apply a transformation to metabolite data
#'
#' @param met_dat Dataset of metabolite concentrations
#' @param met_transform transformation to apply (current options: zscore, log1plus, sqrt)
#'
#' @return Dataset with transformed values (rawValue is pre-transofrmation, value is now transformed value)
#' @export
#'
#' @examples
#' transform_mets(mets_melt, met_transform = "zscore")
transform_mets = function(met_dat, met_transform){
  if(met_transform == "zscore"){
    met_dat[,scaledValue:=scale(value), by=compound]
  } else if(met_transform == "logplus"){
    if(any(met_dat[,value < 0])){
      met_dat[,scaledValue:=value - min(value, na.rm = T), by=compound]
    } else met_dat[,scaledValue:=value]
    met_dat[,scaledValue:=log1p(scaledValue)]
  } else if(met_transform == "sqrt"){
    if(any(met_dat[,value < 0])){
      met_dat[,scaledValue:=value - min(value, na.rm = T), by=compound]
    } else met_dat[,scaledValue:=value]
    met_dat[,scaledValue:=sqrt(scaledValue)]
  } else {
    warning("Transformation option not supported, returning original data")
  }
  if("scaledValue" %in% names(met_dat)){
    setnames(met_dat, c("value", "scaledValue"), c("rawValue", "value"))
  }
  return(met_dat)
}

#' Apply a transformation to CMP data
#'
#' @param cmp_dat Dataset of species-specific metabolic potential scores
#' @param score_transform transformation to apply (current options: zscore, log1plus, sqrt)
#'
#' @return Dataset with transformed values (rawValue is pre-transofrmation, value is now transformed value)
#' @export
#'
#' @examples
#' transform_cmps(indiv_cmps, met_transform = "zscore")
transform_cmps = function(cmp_dat, score_transform){
  if(score_transform == "zscore"){
    #Do it so that they add up to the z scores
    tot_cmps = cmp_dat[,sum(CMP), by=list(compound, Sample)]
    cmp_mean_sd = tot_cmps[,list(mean(V1), sd(V1)), by=compound]
    cmp_dat = merge(cmp_dat, cmp_mean_sd, by = "compound")
    cmp_dat[,scaledValue:=(CMP-V1)/V2, by=compound]
    #cmp_dat[,scaledValue:=scale(CMP), by=list(compound, Species)]
  } else if(score_transform == "logplus"){
    if(any(cmp_dat[,CMP < 0])){
      cmp_dat[,scaledValue:=CMP - min(CMP, na.rm = T), by=compound]
    } else cmp_dat[,scaledValue:=CMP]
    cmp_dat[,scaledValue:=log1p(scaledValue)]
  } else if(score_transform == "sqrt"){
    if(any(cmp_dat[,CMP < 0])){
      cmp_dat[,scaledValue:=CMP - min(CMP, na.rm = T), by=compound]
    } else cmp_dat[,scaledValue:=CMP]
    cmp_dat[,scaledValue:=sqrt(scaledValue)]
  } else {
    warning("Transformation option not supported, returning original data")
  }
  if("scaledValue" %in% names(cmp_dat)){
    setnames(cmp_dat, c("CMP", "scaledValue"), c("rawCMP", "CMP"))
  }
  return(cmp_dat)
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
#' @param add_embl_names Whether to add RefSeq/embl_gems names to the table
#' @param otu_tab Whether to return an OTU table, or just the matched sequences
#' @return Table of alignment results (original sequence, hit ID)
#' @examples
#' map_seqvar(seqs)
#' @export
map_seqvar = function(seqs, repSeqDir = "data/AGORA/", repSeqFile = "agora_NCBI_16S.udb", method = "vsearch", vsearch_path = "vsearch",
                      file_prefix = "seqtemp", seqID = 0.99, add_agora_names = T, add_embl_names = F, otu_tab = F){
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
    if(!file.exists(paste0(repSeqDir, file_prefix, "vsearch_results.txt"))){
      stop("Error: Vsearch failed to align sequences to reference DB. Is your input file formatted correctly?")
    }
    results = fread(paste0(repSeqDir, file_prefix, "vsearch_results.txt"), header = F)
    setnames(results, paste0("V", 1:6), c("seqID", "dbID", "matchPerc", "alnlen", "mism", "gapopens"))
    if(add_embl_names){
      results[,dbID:=gsub("_l.*", "", dbID)]
    }

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
    if(add_agora_names & !add_embl_names){
      #seq_data = fread(paste0(repSeqDir, "agora_NCBItax_processed_nodups.txt"))
      seq_data = fread(paste0(repSeqDir, "AGORA_full_genome_info.txt"))
      seq_matches = merge(seq_matches, seq_data, by.x = "dbID", by.y = "ModelAGORA", all.x = T)
    } else if(add_embl_names){
      seq_data = fread(paste0(repSeqDir, "model_list_processed.txt"))
      seq_matches = merge(seq_matches, seq_data, by.x = "dbID", by.y = "assembly_accession", all.x = T)
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


#' Convert metabolite name table to KEGG metabolite table
#'
#' @import data.table
#' @param network Community metabolic network model
#' @return A network with the direction of some reactions filtered
#' @examples
#' refine_rev_rxns(network)
#' @export
refine_rev_rxns = function(network){
  if(!"Rev" %in% names(network)){
    network2 = get_non_rev_rxns(network, all_rxns = T)
  }
}
