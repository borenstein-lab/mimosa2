### MIMOSA linear model w/residual, then get contributions


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

#melted data table of species and their cmps
add_residuals = function(species_cmps, model_dat, resid_dat){
  resid_dat[,Species:="Residual"]
  setnames(resid_dat, "Resid", "newValue")
  for(x in all_comps){
    species_cmps[compound==x, newValue:=value*model_dat[compound==x, Slope]]
    species_cmps = rbind(species_cmps, resid_dat[compound==x], fill = T)
  }
  species_cmps = species_cmps[!is.na(newValue)]
  return(species_cmps)
}


calculate_var_shares = function(species_contribution_table, valueVar = "newValue"){ #generic, table of values for each speices and sample and compound
  spec_list = species_contribution_table[,unique(Species)]
  spec_table_wide = dcast(species_contribution_table, Sample+compound~Species, value.var = valueVar, fill = 0)
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


