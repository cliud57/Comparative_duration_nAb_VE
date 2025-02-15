###general settings-----
library(dplyr)
library(brms)
library(rstan)
library(emmeans)
library(ggpubr)
library(ggplot2)
library(ggforce)
library(purrr)
library(scales)

setwd("C:/Users/cliud57/SPH Dropbox/CHANG LIU/Waning_review/nAb_VE/submission/202407_nature_series/codes&data")

theme_set(
  
  theme_bw() +
    
    theme(axis.text.x = element_text(size = 36,
                                     colour = "black"),
          axis.text.y = element_text(size = 36,
                                     colour = "black"),
          panel.grid.major = element_blank(),
          panel.background = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "black", linewidth = 1),
          axis.line = element_line(colour = "black", linewidth = 0),
          panel.spacing = unit(3.5, "cm"),
          plot.margin = unit(c(0.8,0.8,0.8,0.8), "cm"),
          strip.background = element_blank(),
          axis.title.x = element_text(size = 36, face = 'bold',colour = "black"),
          axis.title.y = element_text(size = 36, face = 'bold',colour = "black"),
          plot.title = element_text(size = 40, face = 'bold',colour = "black"),
          legend.text = element_text(size = 36, colour = "black"),
          legend.title = element_text(size = 36, colour = "black"),
          plot.title.position = 'plot',
          legend.position = "top")
)

my_palette <- c("#00468BFF","#ED0000FF","#42B540FF",
                "#0099B4FF","#925E9FFF", "#FDAF91FF",
                "#AD002AFF", "#ADB6B6FF","#1B1919FF")

set.seed(225)
random_time <- sample(14:180,100, replace=FALSE) %>%
  sort

set.seed(225)
random_time2 <- sample(7:140,100, replace=FALSE) %>%
  sort

type=c("nab", "VE")
Immunity.condition <- c("2*mRNA", "3*mRNA")
Variant.condition <- c("Delta", "Omicron", "Alpha", "Beta", "Gamma")
Assay.condition <-  c("MNT", "PRNT", "FRNT")
###functions-----
GMT_calculation <- function(time, beta, beta_ll, beta_ul, intercept, intercept_ll, intercept_ul){
  GMT=exp(beta*time+intercept)
  GMT_ll=exp(beta_ll*time+intercept_ll)
  GMT_ul=exp(beta_ul*time+intercept_ul)
  
  return(c(GMT, GMT_ll, GMT_ul))
}

logGMT_reduction <- function(starting_time, time, beta, beta_ll, beta_ul, intercept, intercept_ll, intercept_ul){
  logGMT_reduction=100*((starting_time-time)/(starting_time+intercept/beta))
  logGMT_ll_reduction=100*((starting_time-time)/(starting_time+intercept_ll/beta_ul))
  logGMT_ul_reduction=100*((starting_time-time)/(starting_time+intercept_ul/beta_ll))
  
  return(c(logGMT_reduction, logGMT_ll_reduction, logGMT_ul_reduction))
}
estimated_data <- function(model, estimated_input, type){
  
  set.seed(225)
  
  output=tibble(posterior_epred(model, estimated_input, re_formula = NA))
  
  organized_output=if(type=="nab"){
    data.frame(mean_logGMT=apply(output,2,median),
               sd_logGMT=apply(output,2,sd)) %>%
      mutate(logGMT=mean_logGMT,
             logGMT_ll=mean_logGMT-1.96*sd_logGMT,
             logGMT_ul=mean_logGMT+1.96*sd_logGMT) %>%
      cbind(estimated_input) %>%
      mutate_if(names(.) %in% c("Immunity.type","Tested.variant"), factor)
  } else{
    data.frame(mean_logRR=apply(output,2,median),
               sd_logRR=apply(output,2,sd)) %>%
      mutate(VE=100*(1-exp(mean_logRR)),
             VE_ll=100*(1-exp(mean_logRR+1.96*sd_logRR)),
             VE_ul=100*(1-exp(mean_logRR-1.96*sd_logRR))) %>%
      cbind(estimated_input) %>%
      mutate_if(names(.) %in% c("Vaccination.type","Variant", "Severity"), factor)   
  }
  
  return(organized_output)
}
###data input-----
pre_nab <- read.csv("pre_nab_all_uploading.csv") %>% 
  mutate(Immunity.type=factor(Immunity.type, levels=c("2*mRNA", "2*Inactivated", "Primary non-replicating viral vector",
                                                      "Heterologous primary series","Mild-infection convalescents")),
         Tested.variant=factor(Tested.variant, levels=c("Ancestral", "Alpha", "Beta",
                                                        "Gamma", "Delta", "Omicron")),
         Age.group=factor(Age.group, levels=c("Adults", "Elderly", "All")),
         Assay=factor(Assay, levels=c("MNT", "PRNT", "FRNT")))

post_nab <- read.csv("post_nab_all_uploading.csv") %>% 
  mutate(Immunity.type=factor(Immunity.type, levels=c("3*mRNA", "3*Inactivated", "Non-replicating viral vector booster",
                                                      "Heterologous booster", "2*mRNA convalescents", "2*Inactivated convalescents", "2*Viral vector convalescents")),
         Tested.variant=factor(Tested.variant, levels=c("Ancestral", "Alpha", "Beta",
                                                        "Gamma", "Delta", "Omicron")),
         Age.group=factor(Age.group, levels=c("Adults", "Elderly")),
         Assay=factor(Assay, levels=c("MNT","PRNT", "FRNT")))

pre_VE <- read.csv("pre_VE_all_uploading.csv") %>% 
  mutate(Vaccination.type=factor(Vaccination.type, levels=c("2*mRNA", "Primary series viral vector", "Heterologous primary series")),
         Severity=factor(Severity, levels=c("Mild infection", "Severe", "Death")),
         Variant=factor(Variant, levels=c("Delta","Omicron", "Alpha", "Beta", "Gamma", "Mixture of variants")),
         Exposed_age=factor(Exposed_age, levels=c("Adults", "Elderly", "Children","All")))

post_VE <- read.csv("post_VE_all_uploading.csv") %>% 
  mutate(Vaccination.type=factor(Vaccination.type, levels=c("3*mRNA", "Heterologous booster")),
         Severity=factor(Severity, levels=c("Mild infection", "Severe")),
         Variant=factor(Variant, levels=c("Omicron", "Delta")),
         Exposed_age=factor(Exposed_age, levels=c("Adults", "Elderly","All")))




###data management-----
###extract Omicron BA.1/1.1/2
pre_nab<- pre_nab
post_nab <- post_nab %>% 
  filter(!(PMID=="38396065" & Group.ID %in% c("A2", "A3", "A4", "B2", "B3", "B4", "C2", "C3", "C4", "D2", "D3", "D4")),
         !(PMID=="37400857" & Group.ID=="A1"))
pre_VE <- pre_VE %>% 
  filter(PMID!="36795625",
         !(PMID=="36103455" & Group.ID=="B")) %>% 
  mutate(Severity_simplified=ifelse(Severity=="Mild infection", "Mild infection", "Severe/death"))
post_VE <- post_VE %>% 
  filter(!(PMID=="37839947" & Group.ID=="D"),
         !(PMID=="36103455" & Group.ID=="D"))
pre_VE_mild <- pre_VE %>%  filter(Severity_simplified=="Mild infection")
pre_VE_severe <- pre_VE %>%  filter(Severity_simplified=="Severe/death")
post_VE_mild <- post_VE %>%  filter(Severity=="Mild infection")
post_VE_severe <- post_VE%>%  filter(Severity=="Severe")
###model building and comparison----
###nab-----
models_nab <- list(model_1=bf(logGMT ~ Time*Immunity.type*Tested.variant+Age.group+Assay+(1|PMID)), 
                   model_2=bf(logGMT ~ Time*Immunity.type+Tested.variant+Age.group+Assay+(1|PMID)), 
                   model_3=bf(logGMT ~ Time*Tested.variant+Immunity.type+Age.group+Assay+(1|PMID)), 
                   model_4=bf(logGMT ~ Time+Tested.variant+Immunity.type+Age.group+Assay+(1|PMID)),
                   model_5=bf(logGMT ~ Time:Immunity.type:Tested.variant+Time+Age.group+Assay+(1|PMID)), 
                   model_6=bf(logGMT ~ Time:Immunity.type+Time+Tested.variant+Age.group+Assay+(1|PMID)), 
                   model_7=bf(logGMT ~ Time:Tested.variant+Time+Immunity.type+Age.group+Assay+(1|PMID)))

pre_nab_results_list <- lapply(models_nab,function(model){
  brm(formula=model,
      data=pre_nab,
      prior=c(prior(normal(0,1000),class=b, coef = "Time"),
              prior(normal(0,1000), class=Intercept)),
      seed=225,
      warmup = 1000,
      iter = 2000,
      chains = 4,
      save_pars = save_pars(all = TRUE),
      family = gaussian(),
      control = list(adapt_delta = 0.95,
                     max_treedepth=12))
})  ###1,5 not converged
post_nab_results_list <- lapply(models_nab,function(model){
  brm(formula=model,
      data=post_nab,
      prior=c(prior(normal(0,1000),class=b, coef = "Time"),
              prior(normal(0,1000), class=Intercept)),
      seed=225,
      warmup = 1000,
      iter = 2000,
      chains = 4,
      save_pars = save_pars(all = TRUE),
      family = gaussian(),
      control = list(adapt_delta = 0.95,
                     max_treedepth=12))
}) ###1,5 not converged

###VE-----
models_VE <- list(
  model_1VE=bf(logRR ~ Time*Vaccination.type+Variant+Exposed_age+(1|PMID)),
  model_2VE=bf(logRR ~ Time*Variant+ Vaccination.type+Exposed_age+(1|PMID)),
  model_3VE=bf(logRR ~ Time+Variant+ Vaccination.type+Exposed_age+(1|PMID)),
  model_4VE=bf(logRR ~ Time:Vaccination.type+Time+Variant+Exposed_age+(1|PMID)),
  model_5VE=bf(logRR ~ Time:Variant+Time+ Vaccination.type+Exposed_age+(1|PMID)),
  model_6VE=bf(logRR ~ Time*Variant*Vaccination.type+Exposed_age+(1|PMID)),
  model_7VE=bf(logRR ~ Time:Variant:Vaccination.type+Time+Exposed_age+(1|PMID)))

pre_VE_results_list_mild<- lapply(models_VE,function(model){
  brm(formula=model,
      data=pre_VE_mild,
      prior=c(prior(normal(0,1000),class=b, coef = "Time"),
              prior(normal(0,1000), class=Intercept)),
      seed =225,
      warmup = 1000,
      iter = 2000,
      chains = 4,
      save_pars = save_pars(all = TRUE),
      family = gaussian(),
      control = list(adapt_delta = 0.95,
                     max_treedepth=12))
}) ###6,7 not converged
pre_VE_results_list_severe<- lapply(models_VE,function(model){
  brm(formula=model,
      data=pre_VE_severe,
      prior=c(prior(normal(0,1000),class=b, coef = "Time"),
              prior(normal(0,1000), class=Intercept)),
      seed =225,
      warmup = 1000,
      iter = 2000,
      chains = 4,
      save_pars = save_pars(all = TRUE),
      family = gaussian(),
      control = list(adapt_delta = 0.95,
                     max_treedepth=12))
}) ###6,7 not converged
post_VE_results_list_mild <- lapply(models_VE,function(model){
  brm(formula=model,
      data=post_VE_mild,
      prior=c(prior(normal(0,1000),class=b, coef = "Time"),
              prior(normal(0,1000), class=Intercept)),
      warmup = 1000,
      iter = 2000,
      chains = 4,
      seed=225,
      save_pars = save_pars(all = TRUE),
      family = gaussian(),
      control = list(adapt_delta = 0.99,
                     max_treedepth=12))
}) ###6,7 not converged
post_VE_results_list_severe <- lapply(models_VE,function(model){
  brm(formula=model,
      data=post_VE_severe,
      prior=c(prior(normal(0,1000),class=b, coef = "Time"),
              prior(normal(0,1000), class=Intercept)),
      warmup = 1000,
      iter = 2000,
      chains = 4,
      seed=225,
      save_pars = save_pars(all = TRUE),
      family = gaussian(),
      control = list(adapt_delta = 0.99,
                     max_treedepth=12))
})  ###6,7 not converged


###comparison-----
suffix <- c("_pre_nab", "_post_nab",
            "_pre_VE", "_post_VE",
            "pre_VE_mild", "pre_VE_severe", "post_VE_mild", "post_VE_severe",
            "pre_nab_mRNA", "pre_nab_viral_vector", 
            "pre_VE_mild_mRNA", "pre_VE_mild_viral_vector", "pre_VE_severe_mRNA", "pre_VE_severe_viral_vector")

model_comparison <- function(list, observed_data, suffix){
  names <- c(paste0("model_", seq(1,length(list), 1), suffix))
  
  for(i in 1:length(list)){
    assign(names[i], list[[i]])
  }
  
  ###looic
  set.seed(225)
  looic=lapply(list, loo)
  looic_result=tibble(name=names,
                      looic_number=sapply(looic, function(i)i[["estimates"]]["looic","Estimate"]),
                      se_looic=sapply(looic, function(i)i[["estimates"]]["looic","SE"]),
                      elpd_loo=sapply(looic, function(i)i[["estimates"]]["elpd_loo","Estimate"]))
  comparison=loo_compare(looic)
  ###bayesian_r2
  bayesian_r2 <- lapply(list, bayes_R2) 
  bayesian_r2_result= tibble(name=names, 
                             estimate=sapply(bayesian_r2, function(i)i[,"Estimate"]),
                             Q2.5=sapply(bayesian_r2, function(i)i[,"Q2.5"]),
                             Q97.5=sapply(bayesian_r2, function(i)i[,"Q97.5"])) 
  
  ###rmse
  if(suffix=="_pre_nab"|suffix=="_post_nab"){
    observed=observed_data$logGMT
  }else{observed=observed_data$logRR}
  
  rmse_result=lapply(list, function(i){
    value=posterior_epred(i)
    value_median=sapply(value, median)
    difference=observed-value_median
    
    rmse=sqrt(mean(difference^2))
  }) %>% do.call("rbind", .) 
  
  final=list(looic_value=looic_result,
             looic_comparison=comparison,
             r2_value=bayesian_r2_result,
             rmse_value=rmse_result)
  
  return(final)
}

model_comparison_pre_nab <- model_comparison(pre_nab_results_list[c(2,3,4,6,7)], pre_nab, suffix[1])  ###converged ones
model_comparison_post_nab <- model_comparison(post_nab_results_list[c(2,3,4,6,7)], post_nab, suffix[2])    ###converged ones
model_comparison_pre_VE_mild <- model_comparison(pre_VE_results_list_mild[c(1:5)], pre_VE_mild, suffix[3])  ###converged ones
model_comparison_pre_VE_severe <- model_comparison(pre_VE_results_list_severe[c(1:5)], pre_VE_severe, suffix[3])  ###converged ones
model_comparison_post_VE_mild <- model_comparison(post_VE_results_list_mild[c(1:5)], post_VE_mild, suffix[5])   ###converged ones
model_comparison_post_VE_severe <- model_comparison(post_VE_results_list_severe[c(1:5)], post_VE_severe, suffix[5])  ###converged ones
###figure 1-----
booster_status= c("pre_nab", "post_nab")
figure1_function <- function(model,booster_status){
  if(booster_status=="pre_nab"){
    time=c(14,30,90,180)
    
    input = tibble(Time=rep(time,3*6),
                   Immunity.type=c(rep("2*mRNA",6*length(time)),
                                   rep("2*Inactivated",6*length(time)),
                                   # rep("Heterologous primary series",3*6*length(random_time)),      ###remove insignificant type
                                   # rep("Mild-infection convalescents",3*6*length(random_time)),
                                   rep("Primary non-replicating viral vector",6*length(time))),
                   Assay= rep(c(rep("MNT",6*length(time))),3),
                   Tested.variant=rep(c(rep("Ancestral",length(time)),
                                        rep("Alpha",length(time)),
                                        rep("Beta",length(time)),
                                        rep("Delta",length(time)),
                                        rep("Gamma",length(time)),
                                        rep("Omicron",length(time))),3),
                   Age.group=rep("Adults",3*6*length(time)))
    
    dat = as_tibble(posterior_epred(model, input, re_formula = NA))
  }
  else {
    time=c(7,30,90,180)
    
    input = tibble(Time=rep(time,7*6),
                   Immunity.type=c(rep("3*mRNA",6*length(time)),
                                   rep("3*Inactivated",6*length(time)),
                                   rep("Heterologous booster",6*length(time)),  
                                   rep("Non-replicating viral vector booster",6*length(time)),
                                   rep("2*mRNA convalescents",6*length(time)),
                                   rep("2*Inactivated convalescents",6*length(time)),
                                   rep("2*Viral vector convalescents",6*length(time))),
                   Assay= rep(c(rep("MNT",6*length(time))),7),
                   Tested.variant=rep(c(rep("Ancestral",length(time)),
                                        rep("Alpha",length(time)),
                                        rep("Beta",length(time)),
                                        rep("Delta",length(time)),
                                        rep("Gamma",length(time)),
                                        rep("Omicron",length(time))),7),
                   Age.group=rep("Adults",7*6*length(time)))
    
    dat = as_tibble(posterior_epred(model, input, re_formula = NA))
  } 
  
  median = dat %>% summarise(across(everything(),~median(.))) %>% as.numeric()
  ll = dat %>% summarise(across(everything(),~quantile(., 0.025))) %>% as.numeric()
  ul = dat %>% summarise(across(everything(),~quantile(., 0.975))) %>% as.numeric()
  
  fold_change_dat = as_tibble(t(dat)) %>% 
    mutate(group = rep(1:(ncol(dat) / (4*6)), each = (4*6))) %>% 
    group_by(group) %>% 
    mutate(across(everything(), ~ first(.)-.))
  
  fold_median = exp(as.numeric(apply(fold_change_dat, 1, median)))
  fold_ll= exp(as.numeric(apply(fold_change_dat, 1, function(x)quantile(x,0.025))))
  fold_ul= exp(as.numeric(apply(fold_change_dat, 1, function(x)quantile(x,0.975))))
  
  plot_dat =cbind(input,
                  median, ll, ul,
                  fold_median, fold_ll, fold_ul) %>% 
    mutate(Immunity.abbre=case_when(
      Immunity.type %in% c("2*mRNA", "3*mRNA") ~ "mRNA",
      Immunity.type %in% c("2*Inactivated", "3*Inactivated") ~ "Inactivated",
      Immunity.type %in% c("Primary non-replicating viral vector", "Non-replicating viral vector booster") ~ "Viral\nvector",
      Immunity.type=="Heterologous booster" ~ "Hetero-\nlogous",
      Immunity.type=="2*mRNA convalescents" ~ "Infected\n+\nmRNA",
      Immunity.type=="2*Inactivated convalescents" ~ "Infected\n+\nInactivated",
      Immunity.type=="2*Viral vector convalescents" ~ "Infected\n+\nViral\nvector")) %>% 
    mutate(Tested.variant=ifelse(Tested.variant=="Omicron", "Omicron BA.1/1.1/2", Tested.variant)) %>% 
    mutate(Time=factor(Time, levels=c("7", "14", "30", "90", "180")), 
           Tested.variant=factor(Tested.variant, levels=c("Ancestral", "Alpha", "Beta","Delta", "Gamma", "Omicron BA.1/1.1/2")),
           Immunity.abbre=factor(Immunity.abbre, levels=c("mRNA", "Inactivated", "Viral\nvector", "Hetero-\nlogous",
                                                          "Infected\n+\nmRNA","Infected\n+\nInactivated","Infected\n+\nViral\nvector")))
  
  assign(paste0(booster_status,"_fig1_data"), plot_dat, envir = globalenv())
  
  y_min=round(min(plot_dat$fold_ll),1)
  
  title_number=if(booster_status=="pre_nab"){"A"}else{"B"}
  
  annotation = if(booster_status=="pre_nab"){
    tibble(Immunity.abbre= factor(c("mRNA","Inactivated","Viral\nvector")),
           fold_median=0.8,
           Tested.variant=factor("Ancestral"),
           Time=factor(time[1]))}
  else{
    tibble(Immunity.abbre= factor(c("mRNA","Inactivated","Viral\nvector", "Hetero-\nlogous",
                                    "Infected\n+\nmRNA","Infected\n+\nInactivated","Infected\n+\nViral\nvector")),
           fold_median=0.8,
           Tested.variant=factor("Ancestral"),
           Time=factor(time[1]))
  } 
  
  waning_plot <- ggplot(data=plot_dat,
                        aes(x=Immunity.abbre,y=fold_median, color=Time))+
    geom_point(position=position_dodge(0.7), size=10,
               show.legend = TRUE)+
    geom_errorbar(aes(ymin=fold_ll,
                      ymax=fold_ul),
                  width=0.5,
                  size=2,
                  position=position_dodge(0.7))+
    geom_hline(yintercept = 1, linetype = "dashed", color = "black", linewidth=2)+
    scale_color_manual(values = c("14"=my_palette[1],
                                  "30"=my_palette[2],
                                  "90"=my_palette[3],
                                  "180"=my_palette[4],
                                  "7"=my_palette[6]),
                       drop=FALSE)+
    # coord_cartesian(ylim=c(0.005,10.24),
    #                 expand = FALSE)+
    scale_x_discrete(expand = c(0,0.1))+
    labs(title=title_number,
         x="Immunization type",
         color="Days since complete immunizations")
  
  waning_plot = if(booster_status=="pre_nab"){
    waning_plot+
      scale_y_continuous(limits = c(135,0.5),
                         expand=c(0,0.1),
                         breaks=c(128,32,8,2,1, y_min),
                         labels=c("128", "32", "8", "2", "1", as.character(y_min)),
                         trans=trans_reverser('log2'),
                         name="Fold drop of GMT")+
      facet_wrap(~Tested.variant, nrow=1, scales="free")+
      geom_text(data=annotation, label="Ref", color="black", size=16,  hjust=0.7)+
      theme(strip.text = element_text(size = 36, face="bold"),
            axis.text.x = element_text(size = 36))}
  else{
    waning_plot+
      scale_y_continuous(limits = c(130,0.2),
                         expand=c(0,0.1),
                         breaks=c(128,32,8,2,1, y_min),
                         labels=c("128", "32", "8", "2", "1", as.character(y_min)),
                         trans=trans_reverser('log2'),
                         name="Fold drop of GMT")+
      facet_wrap(~Tested.variant, 
                 nrow=2,
                 scales = "free")+
      geom_text(data=annotation, label="Ref", color="black", size=16,  hjust=0.7)+
      theme(strip.text = element_text(size = 36, face="bold"),
            axis.text.x = element_text(size = 36))
  }
  
  return(waning_plot)
}

pre_nab_plot <- figure1_function(pre_nab_results_list[[2]], booster_status[1])
post_nab_plot<- figure1_function(post_nab_results_list[[4]], booster_status[2])


figure_1 <- ggarrange(pre_nab_plot,
                      post_nab_plot,
                      nrow=2,
                      common.legend = TRUE,
                      heights=c(1,2))+
  theme(plot.background = element_rect(fill = 'white', colour = 'white'))


ggsave("Figure_1.pdf", figure_1, width=45, height=35, limitsize = FALSE, device = "pdf" )

###figure 2-----
figure2_function <- function(model_mild, model_severe, Variant.condition, type){
  if(type=="_pre_VE"){
    time=c(14, 30, 90, 180)
    
    mild_input= tibble(Time=rep(time,3*6),
                       Vaccination.type=c(rep("2*mRNA", 6*length(time)),
                                          rep("Primary series viral vector",6*length(time)),
                                          rep("Heterologous primary series",6*length(time))),
                       Variant=rep(c(rep("Delta",length(time)),
                                     rep("Alpha",length(time)),
                                     rep("Beta",length(time)),
                                     rep("Gamma",length(time)),
                                     rep("Omicron",length(time)),
                                     rep("Mixture of variants",length(time))), 3),
                       Exposed_age=rep("Adults",6*3*length(time)))
    severe_input=tibble(Time=rep(time,3*4),
                        Vaccination.type=c(rep("2*mRNA", 4*length(time)),
                                           rep("Primary series viral vector",4*length(time)),
                                           rep("Heterologous primary series",4*length(time))),
                        Variant=rep(c(rep("Delta",length(time)),
                                      rep("Alpha",length(time)),
                                      rep("Omicron",length(time)),
                                      rep("Mixture of variants",length(time))), 3),
                        Exposed_age=rep("Adults",4*3*length(time)))
    
    dat_mild = as_tibble(posterior_epred(model_mild, mild_input, re_formula = NA))
    dat_severe = as_tibble(posterior_epred(model_severe, severe_input, re_formula = NA))
  }
  else{
    time=c(7, 30, 90, 140)
    
    input = tibble(Time=rep(time,4),
                   Vaccination.type=c(rep("3*mRNA",2*length(time)),
                                      rep("Heterologous booster",2*length(time))),
                   Variant=rep(c(rep("Delta",length(time)),
                                 rep("Omicron",length(time))),2),
                   Exposed_age=rep("Adults",2*2*length(time)))
    
    dat_mild = as_tibble(posterior_epred(model_mild, input, re_formula = NA))
    dat_severe = as_tibble(posterior_epred(model_severe, input, re_formula = NA))
  }
  
  median_mild = dat_mild %>% summarise(across(everything(),~median(.))) %>% as.numeric()
  ll_mild = dat_mild %>% summarise(across(everything(),~quantile(., 0.025))) %>% as.numeric()
  ul_mild = dat_mild %>% summarise(across(everything(),~quantile(., 0.975))) %>% as.numeric()
  
  median_severe = dat_severe %>% summarise(across(everything(),~median(.))) %>% as.numeric()
  ll_severe = dat_severe %>% summarise(across(everything(),~quantile(., 0.025))) %>% as.numeric()
  ul_severe = dat_severe %>% summarise(across(everything(),~quantile(., 0.975))) %>% as.numeric()
  
  if(type=="_pre_VE"){
    plot_dat_mild = cbind(mild_input,median_mild, ll_mild, ul_mild) 
    plot_dat_severe = cbind(severe_input,median_severe, ll_severe, ul_severe)
  }
  else{
    plot_dat_mild = cbind(input,median_mild, ll_mild, ul_mild) 
    plot_dat_severe = cbind(input,median_severe, ll_severe, ul_severe)
  }
  Variant.level=c("Delta" , "Omicron BA.1/1.1/2", "Alpha", "Beta", "Gamma")
  
  plot_dat_mild = plot_dat_mild %>% 
    mutate(Vaccine.abbre = case_when(
      Vaccination.type %in% c("2*mRNA", "3*mRNA") ~ "mRNA",
      Vaccination.type=="Primary series viral vector" ~ "Viral\nvector",
      Vaccination.type %in% c("Heterologous primary series","Heterologous booster") ~"Hetero-\nlogous")) %>% 
    mutate(Variant=ifelse(Variant=="Mixture of variants", "Mixed variants", 
                          ifelse(Variant=="Omicron", "Omicron BA.1/1.1/2", Variant))) %>% 
    mutate(Variant=factor(Variant, levels=c(Variant.level, "Mixed variants")),
           Time=factor(Time, levels=c("7", "14", "30", "90", "140", "180")),
           Vaccine.abbre=factor(Vaccine.abbre, levels=c("mRNA", "Viral\nvector", "Hetero-\nlogous"))) %>% 
    mutate(VE_median=(1-exp(median_mild))*100,
           VE_ll=(1-exp(ul_mild))*100,
           VE_ul=(1-exp(ll_mild))*100) %>% 
    mutate(VE_median_plot=ifelse(VE_median>0, VE_median, 0),
           VE_ll_plot=ifelse(VE_ll>0, VE_ll, 0),
           VE_ul_plot=ifelse(VE_ul>0, VE_ul, 0))
  
  plot_dat_severe = plot_dat_severe %>% 
    mutate(Vaccine.abbre = case_when(
      Vaccination.type %in% c("2*mRNA", "3*mRNA") ~ "mRNA",
      Vaccination.type=="Primary series viral vector" ~ "Viral\nvector",
      Vaccination.type %in% c("Heterologous primary series","Heterologous booster") ~"Hetero-\nlogous")) %>% 
    mutate(Variant=ifelse(Variant=="Mixture of variants", "Mixed variants", 
                          ifelse(Variant=="Omicron", "Omicron BA.1/1.1/2", Variant))) %>% 
    mutate(Variant=factor(Variant, levels=c(Variant.level, "Mixed variants")),
           Time=factor(Time, levels=c("7", "14", "30", "90", "140", "180")),
           Vaccine.abbre=factor(Vaccine.abbre, levels=c("mRNA", "Viral\nvector", "Hetero-\nlogous")))%>% 
    mutate(VE_median=(1-exp(median_severe))*100,
           VE_ll=(1-exp(ul_severe))*100,
           VE_ul=(1-exp(ll_severe))*100) %>% 
    mutate(VE_median_plot=ifelse(VE_median>0, VE_median, 0),
           VE_ll_plot=ifelse(VE_ll>0, VE_ll, 0),
           VE_ul_plot=ifelse(VE_ul>0, VE_ul, 0))
  
  y_intercept_omicron_mild = if(type=="_pre_VE"){
    plot_dat_mild$VE_median_plot[ plot_dat_mild$Vaccine.abbre=="mRNA"&
                                    plot_dat_mild$Variant=="Omicron BA.1/1.1/2" &
                                    plot_dat_mild$Time=="14"]
  }else{
    plot_dat_mild$VE_median_plot[ plot_dat_mild$Vaccine.abbre=="mRNA"&
                                    plot_dat_mild$Variant=="Omicron BA.1/1.1/2" &
                                    plot_dat_mild$Time=="7"]}
  y_intercept_omicron_severe =if(type=="_pre_VE"){
    plot_dat_severe$VE_median_plot[plot_dat_severe$Vaccine.abbre=="mRNA"&
                                     plot_dat_severe$Variant=="Omicron BA.1/1.1/2" &
                                     plot_dat_severe$Time=="14"]
  }else{
    plot_dat_severe$VE_median_plot[plot_dat_severe$Vaccine.abbre=="mRNA"&
                                     plot_dat_severe$Variant=="Omicron BA.1/1.1/2" &
                                     plot_dat_severe$Time=="7"]
  }
  
  title_number=if(type=="_pre_VE"){c("A", "B")}else{c("C","D")}
  # status=if(type=="_pre_VE"){"Primary-immunized"}else{"Booster-immunized"}
  # severe_status=if(type=="_pre_VE"){"severe/fatal infections"}else{"severe infections"}
  
  plot_mild=ggplot(data=plot_dat_mild, 
                   aes(x=Vaccine.abbre,y=VE_median_plot, color=Time))+
    geom_point(position = position_dodge(width=0.5), size=10,
               show.legend = TRUE)+
    geom_errorbar(aes(ymin=VE_ll_plot,
                      ymax=VE_ul_plot),
                  width=0.5,
                  size=2,
                  position = position_dodge(width=0.5))+
    geom_hline(yintercept = y_intercept_omicron_mild, linetype = "dashed", color = "black", linewidth=2)+
    scale_color_manual(values = c("14"=my_palette[1],
                                  "30"=my_palette[2],
                                  "90"=my_palette[3],
                                  "180"=my_palette[4],
                                  "140"=my_palette[5],
                                  "7"=my_palette[6]),
                       drop=FALSE)+
    guides(color = guide_legend(nrow = 1))+
    coord_cartesian(ylim=c(-2,100),
                    expand = TRUE)+
    scale_x_discrete(expand = c(0,0.1))+
    scale_y_continuous(limits = c(-500,100),
                       expand=c(0,2),
                       breaks = seq(0,100,25),
                       labels=as.expression(c(bquote(""<=.(0)),"25","50","75","100")),
                       name="Vaccine effectiveness (%)")+
    labs(title = title_number[1],
      x="Immunization type",
      color="Days since complete immunizations")
  
  plot_severe=ggplot(data=plot_dat_severe, 
                     aes(x=Vaccine.abbre,y=VE_median_plot, color=Time))+
    geom_point(position = position_dodge(width=0.5), size=10,
               show.legend = TRUE)+
    geom_errorbar(aes(ymin=VE_ll_plot,
                      ymax=VE_ul_plot),
                  width=0.5,
                  size=2,
                  position = position_dodge(width=0.5))+
    geom_hline(yintercept = y_intercept_omicron_severe, linetype = "dashed", color = "black", linewidth=2)+
    scale_color_manual(values = c("14"=my_palette[1],
                                  "30"=my_palette[2],
                                  "90"=my_palette[3],
                                  "180"=my_palette[4],
                                  "140"=my_palette[5],
                                  "7"=my_palette[6]),
                       drop=FALSE)+
    guides(color = guide_legend(nrow = 1))+
    coord_cartesian(ylim=c(-2,100),
                    expand = TRUE)+
    scale_x_discrete(expand = c(0,0.1))+
    scale_y_continuous(limits = c(-500,100),
                       expand=c(0,2),
                       breaks = seq(0,100,25),
                       labels=as.expression(c(bquote(""<=.(0)),"25","50","75","100")),
                       name="Vaccine effectiveness (%)")+
    labs( title = title_number[2],
      x="Immunization type",
      color="Days since complete immunizations")+
    facet_wrap(~Variant, nrow=1)+
    theme(strip.text = element_text(size = 36, face="bold"))
  
  if(type=="_pre_VE"){
    plot_mild =  plot_mild +
      facet_wrap(~Variant, nrow=2,
                 scales = "free")+
      theme(strip.text = element_text(size = 36, face="bold"))
    
    combined_plot = ggarrange(plot_mild,
                              plot_severe,
                              nrow=2,
                              common.legend = TRUE,
                              legend="top",
                              heights=c(1.5,1))+
      theme(plot.background = element_rect(fill = 'white', colour = 'white'))
  }
  else{
    plot_mild =  plot_mild +
      facet_wrap(~Variant, nrow=1, scales="free")+
      theme(strip.text = element_text(size = 36, face="bold"))
    
    combined_plot = ggarrange(plot_mild,
                              plot_severe,
                              nrow=1,
                              legend ="none",
                              widths=c(1,1))+
      theme(plot.background = element_rect(fill = 'white', colour = 'white'))
  }
  
  combined_plot
}

pre_VE_plot <- figure2_function(pre_VE_results_list_mild[[3]],
                                pre_VE_results_list_severe[[3]],
                                Variant.condition,
                                "_pre_VE")
post_VE_plot <- figure2_function(post_VE_results_list_mild[[3]],
                                 post_VE_results_list_severe[[3]],
                                 Variant.condition,
                                 "_post_VE")

figure_2 <- ggarrange(pre_VE_plot,
                      post_VE_plot,
                      nrow=2,
                      common.legend = TRUE,
                      legend="top",
                      heights=c(2.7,1))+
  theme(plot.background = element_rect(fill = 'white', colour = 'white'))

ggsave("Figure_2.pdf", figure_2, width=30, height=32, limitsize = FALSE, device = "pdf" )

###figure 3-----
figure3_function <- function(pre_nab_input, post_nab_input,    
                             pre_VE_input_mild,post_VE_input_mild,   
                             pre_nab_model, post_nab_model,
                             pre_VE_mild_model,post_VE_mild_model){
  
  pre_nab_estimated_data <- estimated_data(pre_nab_model, pre_nab_input, type[1])
  post_nab_estimated_data <- estimated_data(post_nab_model, post_nab_input, type[1])
  pre_VE_estimated_data_mild <- estimated_data(pre_VE_mild_model, pre_VE_input_mild, type[2])
  post_VE_estimated_data_mild <- estimated_data(post_VE_mild_model, post_VE_input_mild, type[2])
  
  pre_boost_mRNA_delta_plot<- ggplot()+
    geom_line(data=pre_nab_estimated_data %>%
                filter(Immunity.type=="2*mRNA",
                       Tested.variant=="Delta",
                       Assay=="MNT"),
              aes(x=Time,
                  y=logGMT, color="color1"),size=2)+
    geom_line(data=pre_VE_estimated_data_mild %>%
                filter(Vaccination.type=="2*mRNA",
                       Variant=="Delta"),
              aes(x=Time,
                  y=VE*(log(64)/100)+log(10),color="color2"),linewidth=2)+
    geom_ribbon(data=pre_nab_estimated_data %>%
                  filter(Immunity.type=="2*mRNA",
                         Tested.variant=="Delta",
                         Assay=="MNT"),
                aes(x=Time,
                    ymin=logGMT_ll, ymax=logGMT_ul,
                    fill="color1"), alpha=0.1, color=NA)+
    geom_ribbon(data=pre_VE_estimated_data_mild %>%
                  filter(Vaccination.type=="2*mRNA",
                         Variant=="Delta"),
                aes(x=Time,
                    ymin=VE_ll*(log(64)/100)+log(10), ymax=VE_ul*(log(64)/100)+log(10),
                    fill="color2"), alpha=0.1, color=NA)+
    coord_cartesian(xlim=c(14,180),
                    ylim=c(log(10),log(640)),
                    expand=FALSE,
                    clip="on")+
    scale_x_continuous(limits=c(14,180),
                       expand = c(0,1),
                       breaks=c(14,30,60,90,120,150,180),
                       name="Days since last immunization")+
    scale_y_continuous(limits = c(-300,1500),
                       expand=c(0,0.02),
                       breaks=c(log(10),log(20),log(40),log(80),log(160),log(320),log(640)),
                       labels=c("10","20","40","80","160","320","640"),
                       name=paste0("<span style='font-size: 36pt'>Delta</span><br><br><span style='font-size: 32pt'>Geometric mean titer</span><br>"),
                       sec.axis = sec_axis(~(.-log(10))/(log(64)/100), breaks = seq(0,100,20), name = "Vaccine effectiveness (%)"))+
    scale_color_manual(labels=c("Model estimated nAb GMT",
                                "Model estimated VE"),
                       values = c("color1"=my_palette[7],
                                  "color2"=my_palette[4]))+
    scale_fill_manual(values = c("color1"=my_palette[7],
                                 "color2"=my_palette[4]))+
    guides(fill="none")+
    labs(color="",
         title="A. Primary-immunized nAb GMT and VE against mild\n    delta infections")+
    theme(axis.title.y = ggtext::element_markdown())+
    theme(plot.margin = margin(2,1.5,0.8,0.8,"cm")) 
  
  
  post_boost_mRNA_delta_plot<- ggplot()+
    geom_line(data=post_nab_estimated_data %>%
                filter(Immunity.type=="3*mRNA",
                       Tested.variant=="Delta",
                       Assay=="MNT"),
              aes(x=Time,
                  y=logGMT, color="color1"),size=2)+
    geom_line(data=post_VE_estimated_data_mild %>%
                filter(Vaccination.type=="3*mRNA",
                       Variant=="Delta"),
              aes(x=Time,
                  y=VE*(log(64)/100)+log(10),color="color2"),linewidth=2)+
    geom_ribbon(data=post_nab_estimated_data %>%
                  filter(Immunity.type=="3*mRNA",
                         Tested.variant=="Delta",
                         Assay=="MNT"),
                aes(x=Time,
                    ymin=logGMT_ll, ymax=logGMT_ul,
                    fill="color1"), alpha=0.1, color=NA)+
    geom_ribbon(data=post_VE_estimated_data_mild %>%
                  filter(Vaccination.type=="3*mRNA",
                         Variant=="Delta"),
                aes(x=Time,
                    ymin=VE_ll*(log(64)/100)+log(10), ymax=VE_ul*(log(64)/100)+log(10),
                    fill="color2"), alpha=0.1, color=NA)+
    coord_cartesian(xlim=c(7,140),
                    ylim=c(log(10),log(640)),
                    expand=FALSE,
                    clip="on")+
    scale_x_continuous(limits=c(7,140),
                       expand = c(0,1),
                       breaks=c(7,30,60,90,120,140),
                       name="Days since last immunization")+
    scale_y_continuous(limits = c(-300,3000),
                       expand=c(0,0.02),
                       breaks=c(log(10),log(20),log(40),log(80),log(160),log(320),log(640)),
                       labels=c("10","20","40","80","160","320","640"),
                       name="Geometric mean titer",
                       sec.axis = sec_axis(~(.-log(10))/(log(64)/100), breaks = seq(0,100,20), name = "Vaccine effectiveness (%)"))+
    scale_color_manual(labels=c("Model estimated nAb GMT",
                                "Model etimated VE"),
                       values = c("color1"=my_palette[7],
                                  "color2"=my_palette[4]))+
    scale_fill_manual(values = c("color1"=my_palette[7],
                                 "color2"=my_palette[4]))+
    guides(fill="none")+
    labs(color="",
         title="B. Booster-immunized nAb GMT and VE against mild\n    delta infections")+
    theme(plot.margin = margin(2,0.8,0.8,1.5,"cm"),
          axis.title.y = element_text(vjust = 0.3)) 
  
  
  pre_boost_mRNA_omicron_plot<- ggplot()+
    geom_line(data=pre_nab_estimated_data %>%
                filter(Immunity.type=="2*mRNA",
                       Tested.variant=="Omicron",
                       Assay=="MNT"),
              aes(x=Time,
                  y=logGMT, color="color1"),size=2)+
    geom_line(data=pre_VE_estimated_data_mild %>%
                filter(Vaccination.type=="2*mRNA",
                       Variant=="Omicron"),
              aes(x=Time,
                  y=VE*(log(64)/100)+log(10),color="color2"),linewidth=2)+
    geom_ribbon(data=pre_nab_estimated_data %>%
                  filter(Immunity.type=="2*mRNA",
                         Tested.variant=="Omicron",
                         Assay=="MNT"),
                aes(x=Time,
                    ymin=logGMT_ll, ymax=logGMT_ul,
                    fill="color1"), alpha=0.1, color=NA)+
    geom_ribbon(data=pre_VE_estimated_data_mild %>%
                  filter(Vaccination.type=="2*mRNA",
                         Variant=="Omicron"),
                aes(x=Time,
                    ymin=VE_ll*(log(64)/100)+log(10), ymax=VE_ul*(log(64)/100)+log(10),
                    fill="color2"), alpha=0.1, color=NA)+
    coord_cartesian(xlim=c(14,180),
                    ylim=c(log(10),log(640)),
                    expand=FALSE,
                    clip="on")+
    scale_x_continuous(limits=c(14,180),
                       expand = c(0,1),
                       breaks=c(14,30,60,90,120,150,180),
                       name="Days since last immunization")+
    scale_y_continuous(limits = c(-700,1500),
                       expand=c(0,0.02),
                       breaks=c(log(10),log(20),log(40),log(80),log(160),log(320),log(640)),
                       labels=c("10","20","40","80","160","320","640"),
                       name=paste0("<span style='font-size: 36pt'>Omicron BA.1&#47;1.1&#47;2</span><br><br><span style='font-size: 32pt'>Geometric mean titer</span><br>"),
                       sec.axis = sec_axis(~(.-log(10))/(log(64)/100), breaks = seq(0,100,20), name = "Vaccine effectiveness (%)"))+
    scale_color_manual(labels=c("Model estimated nAb GMT",
                                "Model estimated VE"),
                       values = c("color1"=my_palette[7],
                                  "color2"=my_palette[4]))+
    scale_fill_manual(values = c("color1"=my_palette[7],
                                 "color2"=my_palette[4]))+
    guides(fill="none")+
    labs(color="",
         title="C. Primary-immunized nAb GMT and VE against mild\n    BA.1/1.1/2 infections")+
    theme(axis.title.y = ggtext::element_markdown())+
    theme(plot.margin = margin(0.8,1.5,0.8,0.8,"cm")) 
  
  
  post_boost_mRNA_omicron_plot<- ggplot()+
    geom_line(data=post_nab_estimated_data %>%
                filter(Immunity.type=="3*mRNA",
                       Tested.variant=="Omicron",
                       Assay=="MNT"),
              aes(x=Time,
                  y=logGMT, color="color1"),size=2)+
    geom_line(data=post_VE_estimated_data_mild %>%
                filter(Vaccination.type=="3*mRNA",
                       Variant=="Omicron"),
              aes(x=Time,
                  y=VE*(log(64)/100)+log(10),color="color2"),linewidth=2)+
    geom_ribbon(data=post_nab_estimated_data %>%
                  filter(Immunity.type=="3*mRNA",
                         Tested.variant=="Omicron",
                         Assay=="MNT"),
                aes(x=Time,
                    ymin=logGMT_ll, ymax=logGMT_ul,
                    fill="color1"), alpha=0.1, color=NA)+
    geom_ribbon(data=post_VE_estimated_data_mild %>%
                  filter(Vaccination.type=="3*mRNA",
                         Variant=="Omicron"),
                aes(x=Time,
                    ymin=VE_ll*(log(64)/100)+log(10), ymax=VE_ul*(log(64)/100)+log(10),
                    fill="color2"), alpha=0.1, color=NA)+
    coord_cartesian(xlim=c(7,140),
                    ylim=c(log(10),log(640)),
                    expand=FALSE,
                    clip="on")+
    scale_x_continuous(limits=c(7,140),
                       expand = c(0,1),
                       breaks=c(7,30,60,90,120,140),
                       name="Days since last immunization")+
    scale_y_continuous(limits = c(-300,1500),
                       expand=c(0,0.02),
                       breaks=c(log(10),log(20),log(40),log(80),log(160),log(320),log(640)),
                       labels=c("10","20","40","80","160","320","640"),
                       name="Geometric mean titer",
                       sec.axis = sec_axis(~(.-log(10))/(log(64)/100), breaks = seq(0,100,20), name = "Vaccine effectiveness (%)"))+
    scale_color_manual(labels=c("Model estimated nAb GMT",
                                "Model estimated VE"),
                       values = c("color1"=my_palette[7],
                                  "color2"=my_palette[4]))+
    scale_fill_manual(values = c("color1"=my_palette[7],
                                 "color2"=my_palette[4]))+
    guides(fill="none")+
    labs(color="",
         title="D. Booster-immunized nAb GMT and VE against mild\n    BA.1/1.1/2 infections")+
    theme(plot.margin = margin(0.8,0.8,0.8,1.5,"cm"),
          axis.title.y = element_text(vjust = 0.3)) 
  
  
  combined_plot <- ggarrange(pre_boost_mRNA_delta_plot,
                             post_boost_mRNA_delta_plot,
                             pre_boost_mRNA_omicron_plot,
                             post_boost_mRNA_omicron_plot,
                             nrow = 2,
                             ncol=2,
                             align="hv",
                             common.legend = TRUE,
                             widths = c(0.505,0.495),
                             labels = c("Primary immunization", "Booster immunization"),
                             font.label = list(size=36),hjust=c(-1.3,-1.3),vjust=c(1.5,1.5))+
    theme(plot.background = element_rect(fill = 'white', colour = 'white'))
  
  combined_plot
} 

pre_nab_fig3_input <- data.frame(Time=rep(random_time,15),
                                 Immunity.type=rep(c(rep("2*mRNA",100),
                                                     rep("2*Inactivated",100),
                                                     rep("Primary non-replicating viral vector",100)),5),
                                 Tested.variant=c(rep("Delta",300),
                                                  rep("Omicron",300),
                                                  rep("Alpha",300),
                                                  rep("Beta",300),
                                                  rep("Gamma",300)),
                                 Age.group=rep("Adults",1500),
                                 Assay=rep("MNT",1500))
post_nab_fig3_input <- data.frame(Time=rep(random_time2,6),
                                  Immunity.type=rep(c(rep("3*mRNA",100),
                                                      rep("3*Inactivated",100),
                                                      rep("Heterologous booster",100)),2),
                                  Tested.variant=c(rep("Delta",300),
                                                   rep("Omicron",300)),
                                  Age.group=rep("Adults",600),
                                  Assay=rep("MNT",600))
pre_VE_fig3_mild_input <- tibble(Time=rep(random_time,18),
                                 Vaccination.type=c(rep("2*mRNA", 6*length(random_time)),
                                                    rep("Primary series viral vector",6*length(random_time)),
                                                    rep("Heterologous primary series",6*length(random_time))),
                                 Variant=rep(c(rep("Delta",length(random_time)),
                                               rep("Alpha",length(random_time)),
                                               rep("Beta",length(random_time)),
                                               rep("Gamma",length(random_time)),
                                               rep("Omicron",length(random_time)),
                                               rep("Mixture of variants",length(random_time))), 3),
                                 Exposed_age=rep("Adults",6*3*length(random_time)))
post_VE_fig3_subset_input <- tibble(Time=rep(random_time2,4),
                                    Vaccination.type=c(rep("3*mRNA",2*length(random_time2)),
                                                       rep("Heterologous booster",2*length(random_time2))),
                                    Variant=rep(c(rep("Delta",length(random_time2)),
                                                  rep("Omicron",length(random_time2))),2),
                                    Exposed_age=rep("Adults",2*2*length(random_time2)))



figure_3 <-figure3_function  (pre_nab_fig3_input,
                              post_nab_fig3_input,
                              pre_VE_fig3_mild_input,
                              post_VE_fig3_subset_input,
                              pre_nab_results_list[[2]] ,
                              post_nab_results_list[[4]] ,
                              pre_VE_results_list_mild[[3]],
                              post_VE_results_list_mild[[3]])



ggsave("Figure_3.pdf",figure_3, width=35, height=25, limitsize = FALSE,device="pdf")

###figure 4-----
figure4_function <- function(mild_correlation_data,
                             severe_correlation_data,         
                             mild_quantile_list,
                             severe_quantile_list,
                             Variant.condition){
  
  correlation_list = list(mild=do.call("rbind", mild_correlation_data),
                          severe=do.call("rbind", severe_correlation_data))
  correlation_quantile_list=list(mild=mild_quantile_list,
                                 severe=severe_quantile_list)
  
  corA= sprintf("%.3f",cor(correlation_list[[1]]$GMT, correlation_list[[1]]$VE, method = "spearman"))
  corB= sprintf("%.3f",cor(correlation_list[[2]]$GMT, correlation_list[[2]]$VE, method = "spearman"))
  
  annotationA =data.frame(GMT=160, VE=20)
  annotationB =data.frame(GMT=160, VE=20)
  
  plotA = ggplot()+
    geom_ribbon(data= correlation_quantile_list[[1]][[Variant.condition[1]]][[2]],
                aes(y=VE_median, x=GMT_median,
                    xmin=GMT_ll,xmax=GMT_ul,fill="color4"),alpha=0.3, color=NA)+
    geom_ribbon(data= correlation_quantile_list[[1]][[Variant.condition[2]]][[2]],
                aes(y=VE_median,x=GMT_median,
                    xmin=GMT_ll,xmax=GMT_ul,fill="color5"),alpha=0.3, color=NA)+
    geom_ribbon(data= correlation_quantile_list[[1]][[Variant.condition[3]]][[2]],
                aes(y=VE_median,x=GMT_median,
                    xmin=GMT_ll,xmax=GMT_ul,fill="color1"),alpha=0.3, color=NA)+
    geom_ribbon(data=correlation_quantile_list[[1]][[Variant.condition[4]]][[2]],
                aes(y=VE_median,x=GMT_median,
                    xmin=GMT_ll,xmax=GMT_ul,fill="color3"),alpha=0.3, color=NA)+
    geom_ribbon(data=correlation_quantile_list[[1]][[Variant.condition[5]]][[2]],
                aes(y=VE_median,x=GMT_median,
                    xmin=GMT_ll,xmax=GMT_ul,fill="color2"),alpha=0.3, color=NA)+
    geom_line(data=correlation_quantile_list[[1]][[Variant.condition[1]]][[1]],
              aes(x=GMT, y=VE,color="color4"),alpha=1, linewidth=7)+
    geom_line(data=correlation_quantile_list[[1]][[Variant.condition[2]]][[1]],
              aes(x=GMT, y=VE, color="color5"),alpha=1, linewidth=7)+
    geom_line(data=correlation_quantile_list[[1]][[Variant.condition[3]]][[1]],
              aes(x=GMT, y=VE, color="color1"),alpha=1, linewidth=7)+
    geom_line(data=correlation_quantile_list[[1]][[Variant.condition[4]]][[1]],
              aes(x=GMT, y=VE, color="color3"),alpha=1, linewidth=7)+
    geom_line(data=correlation_quantile_list[[1]][[Variant.condition[5]]][[1]],
              aes(x=GMT, y=VE, color="color2"),alpha=1, linewidth=7)+
    coord_cartesian(xlim=c(5,640),
                    ylim=c(0,105),
                    expand = FALSE,
                    clip="on")+
    scale_x_continuous(limits = c(0,1280),
                       breaks=c(10,20,40,80,160,320,640),
                       trans="log2",
                       name="Geometric mean titer")+
    scale_y_continuous(limits=c(-80000,105),
                       breaks = seq(0,100,by=20),
                       name="Vaccine effectiveness(%)")+
    scale_color_manual(values = c("color1"=my_palette[1],
                                  "color2"=my_palette[2],
                                  "color3"=my_palette[3],
                                  "color4"=my_palette[4],
                                  "color5"=my_palette[5]),
                       labels =c("color4"="Delta",
                                 "color5"="Omicron BA.1/1.1/2",
                                 "color1"="Alpha",
                                 "color3"="Beta",
                                 "color2"="Gamma"))+
    scale_fill_manual(values = c("color1"=my_palette[1],
                                 "color2"=my_palette[2],
                                 "color3"=my_palette[3],
                                 "color4"=my_palette[4],
                                 "color5"=my_palette[5]),
                      guide="none")+
    theme(plot.margin =margin(2,1,0.8,0.8,"cm") )+
    geom_text(data=annotationA,
              aes(x=GMT, y=VE),
              label=paste('rho',"==",corA),
              parse=TRUE,
              color="black", size=16, hjust=0.7)+
    labs(title="A",
         color="Variant")
  
  plotB = ggplot()+
    geom_ribbon(data= correlation_quantile_list[[2]][[Variant.condition[1]]][[2]],
                aes(y=VE_median, x=GMT_median,
                    xmin=GMT_ll,xmax=GMT_ul,fill="color4"),alpha=0.3, color=NA)+
    geom_ribbon(data= correlation_quantile_list[[2]][[Variant.condition[2]]][[2]],
                aes(y=VE_median,x=GMT_median,
                    xmin=GMT_ll,xmax=GMT_ul,fill="color5"),alpha=0.3, color=NA)+
    geom_ribbon(data= correlation_quantile_list[[2]][[Variant.condition[3]]][[2]],
                aes(y=VE_median,x=GMT_median,
                    xmin=GMT_ll,xmax=GMT_ul,fill="color1"),alpha=0.3, color=NA)+
    geom_line(data=correlation_quantile_list[[2]][[Variant.condition[1]]][[1]],
              aes(x=GMT, y=VE,color="color4"),alpha=1, linewidth=7)+
    geom_line(data=correlation_quantile_list[[2]][[Variant.condition[2]]][[1]],
              aes(x=GMT, y=VE, color="color5"),alpha=1, linewidth=7)+
    geom_line(data=correlation_quantile_list[[2]][[Variant.condition[3]]][[1]],
              aes(x=GMT, y=VE, color="color1"),alpha=1, linewidth=7)+
    coord_cartesian(xlim=c(5,640),
                    ylim=c(0,105),
                    expand = FALSE,
                    clip="on")+
    scale_x_continuous(limits = c(0,1280),
                       breaks=c(10,20,40,80,160,320,640),
                       trans="log2",
                       name="Geometric mean titer")+
    scale_y_continuous(limits=c(-80000,105),
                       breaks = seq(0,100,by=20),
                       name="Vaccine effectiveness(%)")+
    scale_color_manual(values = c("color1"=my_palette[1],
                                  "color4"=my_palette[4],
                                  "color5"=my_palette[5]),
                       labels =c("color4"="Delta",
                                 "color5"="Omicron BA.1/1.1/2",
                                 "color1"="Alpha"))+
    scale_fill_manual(values = c("color1"=my_palette[1],
                                 "color4"=my_palette[4],
                                 "color5"=my_palette[5]),
                      guide="none")+
    theme(plot.margin =margin(2,1,0.8,0.8,"cm") )+
    geom_text(data=annotationB,
              aes(x=GMT, y=VE),
              label=paste('rho',"==",corB),
              parse=TRUE,
              color="black", size=16, hjust=0.7)+
    labs(title="B",
         color="Variant")
  
  
  combined_plot <- ggarrange(plotA, plotB,
                             nrow = 1,
                             align="v",
                             common.legend = TRUE,
                             widths = c(1,1),
                             labels = c("Mild", "Severe/death"),
                             font.label = list(size=36),hjust=c(-6.3,-1.7),vjust=c(1.5,1.5))+
    theme(plot.background = element_rect(fill = 'white', colour = 'white'),
          panel.spacing = unit(5, "cm"))
}

correlation_data_function<- function(estimated_nab_data,
                                     estimated_VE_data,
                                     Immunity.condition,
                                     Variant.condition,
                                     Assay.condition,
                                     time_seq){
  
  correlation_list <- replicate(1000,
                                data.frame(time=time_seq,
                                           logGMT=numeric() %>% 
                                             rep(NA,100),
                                           logRR=numeric() %>% 
                                             rep(NA,100),
                                           logGMT_change=numeric() %>% 
                                             rep(NA,100),
                                           logRR_change=numeric() %>% 
                                             rep(NA,100),
                                           logGMT_change_percent=numeric() %>% 
                                             rep(NA,100),
                                           logRR_change_percent=numeric() %>% 
                                             rep(NA,100)),
                                simplify=FALSE)
  
  min_time = min(time_seq)
  
  for(i in seq_len(length(correlation_list))){
    nab_data=estimated_nab_data %>%
      filter(Immunity.type==Immunity.condition,
             Tested.variant==Variant.condition,
             Assay==Assay.condition)
    VE_data=estimated_VE_data %>%
      filter(Vaccination.type==Immunity.condition,
             Variant==Variant.condition)
    
    for (j in seq_len(nrow(correlation_list[[i]]))){
      correlation_list[[i]]$logGMT[j]=rnorm(1, nab_data$mean_logGMT[j], nab_data$sd_logGMT[j])
      
      correlation_list[[i]]$logRR[j]=rnorm(1, VE_data$mean_logRR[j], VE_data$sd_logRR[j])
      
    }
    
    logGMT_ref=correlation_list[[i]]$logGMT[correlation_list[[i]]$time==min_time]
    logRR_ref=correlation_list[[i]]$logRR[correlation_list[[i]]$time==min_time]
    
    correlation_list[[i]]$logGMT_change=correlation_list[[i]]$logGMT-logGMT_ref
    correlation_list[[i]]$logRR_change=correlation_list[[i]]$logRR-logRR_ref
    correlation_list[[i]]$logGMT_change_percent=(-correlation_list[[i]]$logGMT+logGMT_ref)/logGMT_ref
    correlation_list[[i]]$logRR_change_percent=(-correlation_list[[i]]$logRR+logRR_ref)/logRR_ref
  }
  
  combined_correlation_list <- map_dfr(correlation_list, ~.) %>%
    mutate(GMT=exp(logGMT),
           VE=100*(1-exp(logRR)))
  
  return(combined_correlation_list)
  
}

correlation_data_quantile_subset <- function(estimated_correaltion_data, 
                                             estimated_nab_data,
                                             estimated_VE_data,
                                             Immunity.condition,
                                             Variant.condition,
                                             Assay.condition){
  
  estimated_results <- tibble(GMT=exp(estimated_nab_data$logGMT[estimated_nab_data$Immunity.type==Immunity.condition 
                                                                & estimated_nab_data$Tested.variant==Variant.condition
                                                                &estimated_nab_data$Assay==Assay.condition]),
                              VE=estimated_VE_data$VE[estimated_VE_data$Vaccination.type==Immunity.condition 
                                                      & estimated_VE_data$Variant==Variant.condition]) 
  grouped_correlated_data <- estimated_correaltion_data %>% 
    group_by(time) %>% 
    mutate(GMT_ll=exp(quantile(logGMT, 0.025)),
           GMT_ul=exp(quantile(logGMT, 0.975)),
           VE_ll=(1-exp(quantile(logRR, 0.975)))*100,
           VE_ul=(1-exp(quantile(logRR, 0.025)))*100,
           logGMT_change_ll=quantile(logGMT_change, 0.025),
           logGMT_change_ul=quantile(logGMT_change, 0.975),
           logGMT_change_percent_ll=quantile(logGMT_change_percent, 0.025),
           logGMT_change_percent_ul=quantile(logGMT_change_percent, 0.975),
           logRR_change_ll=quantile(logRR_change, 0.025),
           logRR_change_ul=quantile(logRR_change, 0.975),
           logRR_change_percent_ll=quantile(logRR_change_percent, 0.025),
           logRR_change_percent_ul=quantile(logRR_change_percent, 0.975)) %>% 
    ungroup() %>%
    slice(1:100)%>%
    select(-c(logGMT,logRR,logGMT_change,logRR_change,logGMT_change_percent,logRR_change_percent, GMT, VE)) %>%
    cbind(GMT_median=estimated_results$GMT,
          VE_median=estimated_results$VE) %>% 
    mutate(logGMT_change_median=log(GMT_median)-log(GMT_median[time==min(time)]),
           logRR_change_median=log(1-VE_median/100)-log(1-VE_median[time==min(time)]/100))
  
  return(list(estimated_results,grouped_correlated_data))
  
}
###data preparation
pre_VE_fig4_severe_input <- tibble(Time=rep(random_time,12),
                                   Vaccination.type=c(rep("2*mRNA", 4*length(random_time)),
                                                      rep("Primary series viral vector",4*length(random_time)),
                                                      rep("Heterologous primary series",4*length(random_time))),
                                   Variant=rep(c(rep("Delta",length(random_time)),
                                                 rep("Alpha",length(random_time)),
                                                 rep("Omicron",length(random_time)),
                                                 rep("Mixture of variants",length(random_time))), 3),
                                   Exposed_age=rep("Adults",4*3*length(random_time)))


pre_nab_estimated_fig4_data <- estimated_data(pre_nab_results_list[[2]], pre_nab_fig3_input , type[1])   ###same input as fig.3
pre_VE_estimated_fig4_data_mild <- estimated_data(pre_VE_results_list_mild[[3]], pre_VE_fig3_mild_input, type[2]) ###same input as fig.3
pre_VE_estimated_fig4_data_severe <- estimated_data(pre_VE_results_list_severe[[3]], pre_VE_fig4_severe_input, type[2])


pre_mild_correlation_list <-lapply(Variant.condition, function(v){
  set.seed(225)
  dat = correlation_data_function(pre_nab_estimated_fig4_data, pre_VE_estimated_fig4_data_mild,
                                  Immunity.condition[1],v,Assay.condition[1], random_time)
  
  return(dat)})
names(pre_mild_correlation_list) = Variant.condition

pre_severe_correlation_list<- lapply(Variant.condition[1:3], function(v){
  set.seed(225)
  dat = correlation_data_function(pre_nab_estimated_fig4_data, pre_VE_estimated_fig4_data_severe,
                                  Immunity.condition[1],v,Assay.condition[1], random_time)
  
  return(dat)})
names(pre_severe_correlation_list)= Variant.condition[1:3]

pre_mild_correlation_quantile_list <- lapply(Variant.condition, function(v){
  correlation_data_quantile_subset (pre_mild_correlation_list[[v]], 
                                    pre_nab_estimated_fig4_data, 
                                    pre_VE_estimated_fig4_data_mild,
                                    Immunity.condition[1],
                                    v,
                                    Assay.condition[1])
})
names(pre_mild_correlation_quantile_list)=Variant.condition 

pre_severe_correlation_quantile_list <- lapply(Variant.condition[1:3], function(v){
  correlation_data_quantile_subset (pre_severe_correlation_list [[v]], 
                                    pre_nab_estimated_fig4_data, 
                                    pre_VE_estimated_fig4_data_severe,
                                    Immunity.condition[1],
                                    v,
                                    Assay.condition[1])
})
names(pre_severe_correlation_quantile_list)=Variant.condition[1:3] 
###all set



figure_4 <- figure4_function(pre_mild_correlation_list,
                             pre_severe_correlation_list,          
                             pre_mild_correlation_quantile_list,
                             pre_severe_correlation_quantile_list,
                             Variant.condition)

ggsave("Figure_4.pdf",figure_4, width=25, height=10, limitsize = FALSE,device="pdf")

