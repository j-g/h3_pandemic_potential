################################################################################
# Functions to run all of analysis
################################################################################

# Load libraries
library(mgcv)
library(cenGAM)
library(plyr)
library(dplyr)
library(rstan)

################################################################################
# Load MCMC output from titer model manuscript to take draws for susceptibility
# model
################################################################################

# Read in MCMC data
load("../data/titer_model_MCMC.RData")

# Get distributions
df<-extract(fitchi10)
s0<-df$strain0
pt<-df$peak_titer
tb<-df$t_buildup
tw<-df$t_wane
amag<-df$as_mag
atim<-df$as_time
tmax<-46
params<-data.frame(s0, pt, tb, tw, amag, atim, tmax)

################################################################################
# Coudeville params for susceptability as function of titer
# https://doi.org/10.1186/1471-2288-10-18
################################################################################

alpha<-2.844
beta<-1.299
alpha_min<-2.23
beta_min<-1
alpha_max<-3.36
beta_max<-1.69

# Assume posterior distribution is Gaussian and reverse engineer standard
# deviations
sd_a<-(alpha_max-alpha_min)/3.92
sd_b<-(beta_max-beta_min)/3.92

# Group params together
coude_params<-c(alpha, beta, sd_a, sd_b)

################################################################################
# Model functions
################################################################################

#Take average value per bin starting at 0.
boxcar_avg<-function(df, bin_width, group, value)
{
  df$group_mid<-round_any(df[,group]-(bin_width/2), bin_width)+(bin_width/2)
  df <-setNames(aggregate(df[,value], by=list(df[,"group_mid"]), FUN=mean), c(group, value))
  return(df)
}

#Turn Log Titer into susceptibility using Coudeville relation
protection_func<-function(LogTiter, alpha, beta)
{
  titer<-10*(2^(LogTiter-1))
  x<-1/(1+exp(beta*(log(titer)-alpha)))
  return(x)
}

# Create susc from random Coudeville params (alpha, beta) and random susc model
# params. Bin width and max age can vary but =5 and 70 for most contact surveys
# StrainYear = year of strain under question, nonHAI =F if not considering broad
# immunity, mag=magnitude of non HAI effects and StudyYear=year in question.
# output is dataframe of Age, susceptibility, Time since Circulation (TSC)
predict_susc<-function(params, alpha, beta, bin_width, max_age, StrainYear, nonHAI, mag, StudyYear)
{
  # Sequence of ages in bins of 1 year
  Age<-seq(from=0, to=max_age+bin_width, by=1)
  
  TSC<-StudyYear-StrainYear
  YoB<-StudyYear-Age
  
  # Time Since First Exposure (TSFE) as defined in titer model manuscript
  TSFE<-Age-TSC
  TSFE[YoB<=1968]<-StrainYear-1968
  
  newd<-data.frame(Age, TSC, TSFE)
  
  # Predict susc from params
  s0<-params$s0
  pt<-params$pt
  tb<-params$tb
  tw<-params$tw
  amag<-params$amag
  atim<-params$atim
  tmax<-params$tmax
  
  # Antigenic seiority term
  ant_sen_c<-sqrt(2*pi)*atim*amag/tmax*(0.5-pnorm(tmax, 0, atim))
  Y_as1<-amag*exp(-0.5*(TSFE/atim)^2)+ant_sen_c
  Y_as2<-(amag+ant_sen_c)+(TSFE/3)
  # Smoothing width for pre birth titers. Experimentally these flatline since
  # the HAI titer is 0. We allow them tio transition to a linear decay.
  wid<-2
  Y_as<-((tanh(TSFE/wid)+1)*0.5)*Y_as1+(1-(tanh(TSFE/wid)+1)*0.5)*Y_as2
  
  # Time since ciculation term
  Y_tsc<-s0+pt*(1-exp(-TSC/tb))-(TSC/tw)
  
  # Modelled titers
  newd$LogTiter<-Y_as+Y_tsc
  
  
  # Take weighted average of susceptibility estimates over distribution of
  # individual variation in the titer model
  
  # Standard deviation of the log titer model ~ individual variation
  sd<-1.325759
  
  # Break distribution into 100 slices and calculate weighted average over 4sd's
  n_int<-100
  dlt<-4*sd/n_int
  
  # Add up contributions to get weighted susc
  susc<-rep(0, nrow(newd))
  for(lt in seq(-2*sd, 2*sd, dlt))
  {
    susc<-susc+dlt*dnorm(lt, mean=0, sd=sd)*sapply( (newd$LogTiter+lt), alpha=alpha, beta=beta, protection_func)
  }
  
  # If nonHAI=T, add in non specific component
  if(nonHAI)
  {
    phi<-1-mag*(Age/70)
    
    protect<-phi*(1-susc)+(1-phi)*0.5
    
    susc<-1-protect
    
    susc[susc<0]<-0
  }
  
  newd$susc<-susc
  
  #Average yearly alues over (typically 5 year) age bins
  newd<-boxcar_avg(newd, bin_width, "Age", "susc")
  newd$TSC<-TSC
  
  return(newd)
}

#Calculate R0 given next generaion matrix
calculate_R0 <- function(cont, susc)
{
  #next gen matrix
  M<-cont*susc
  e<-eigen(M)
  
  return(max(Re(e$values)))
}

# Takes list of countries (countries), params (params, coude_params, bin_width, max_age)
# vector of strain years (StrainYears), options for non HAI protection for adults
# (nonHAI and degree through mag) and StudyYear. Outputs estimates of R_eff for
# N_rep repetitions for each config.
rand_sweep<-function(countries, params, coude_params, bin_width, max_age, N_rep, StrainYears, nonHAI, mag, StudyYear)
{
  alpha<-coude_params[1]
  beta<-coude_params[2]
  sd_a<-coude_params[3]
  sd_b<-coude_params[4]
  
  # Declare data structures
  N<-N_rep*length(StrainYears)*length(countries)
  R0<-rep(0, N)
  country_name<-rep("a", N)
  rep<-rep(0, N)
  k<-rep(0, N)
  med_age<-rep(0, N)
  strain_yr<-rep(0,N)
  
  j<-1
  for(i in 1:N_rep)
  {
    t0<-Sys.time()
    cat(i, "\t")
    
    # Generate random params
    # Coudeville susc
    alpha_rand<-rnorm(1, mean = alpha, sd = sd_a)
    beta_rand<-rnorm(1, mean = beta, sd = sd_b)
    
    # Titer model
    isamp<-sample(1:nrow(params), 1)
    psamp<-params[isamp,]
    
    for(StrainYear in StrainYears)
    {
      # Predict susc
      susc<-predict_susc(psamp, alpha_rand, beta_rand, bin_width, max_age, StrainYear, nonHAI, mag, StudyYear)
      susc<-susc$susc
      
      # Loop over countries
      for(Country in countries)
      {
        R0[j]<-calculate_R0(Country$cont, susc)
        country_name[j]<-Country$country
        rep[j]<-i
        k[j]<-Country$k
        med_age[j]<-Country$med_age
        strain_yr[j]<-StrainYear
        j<-j+1
      } 
    }
    t1<-Sys.time()
    cat(difftime(t1, t0), "\t\t")
  }
  df_out<-data.frame(R0, country_name, k, med_age, rep, strain_yr)
  
  return(df_out)
}

# Same as rand_sweep but for average parameter values
avg_sweep<-function(countries, params, coude_params, bin_width, max_age, StrainYears, nonHAI, mag, StudyYear)
{
  alpha<-coude_params[1]
  beta<-coude_params[2]
  
  N<-length(StrainYears)*length(countries)
  
  R0<-rep(0, N)
  country_name<-rep("a", N)
  rep<-rep(0, N)
  k<-rep(0, N)
  med_age<-rep(0, N)
  strain_yr<-rep(0,N)
  
  j<-1
  
  pars<-as.data.frame(t(colMeans(params)))
  
  for(StrainYear in StrainYears)
  {
    susc<-predict_susc(pars, alpha, beta, bin_width, max_age, StrainYear, nonHAI, mag, StudyYear)
    susc<-susc$susc
    
    for(Country in countries)
    {
      R0[j]<-calculate_R0(Country$cont, susc)
      country_name[j]<-Country$country
      k[j]<-Country$k
      med_age[j]<-Country$med_age
      strain_yr[j]<-StrainYear
      j<-j+1
    } 
  }
  
  df_out<-data.frame(R0, country_name, k, med_age, rep, strain_yr)
  
  return(df_out)
}

# Make new country, combining behav and demog of input countries. Baseline
# connectivity depends on normalization convention: norm="pand" for pandemic.
mix<-function(country_behav, country_demog, norm)
{
  cont<-(country_demog$rand)*(country_behav$behav)
  
  # Choice to normalize baseline connectivity so that R_eff fits an expected
  # pandemic limit, or so that connectivity of the behavioural country is
  # conserved.
  if(norm=="pand")
  {
    #Normalize to pandemic limit
    e<-eigen(cont)
    R0<-(max(Re(e$values)))
    cont<-2.5*cont/R0
  }
  else
  {
    # Conserve connectivity
    connect<-connectivity(cont, country_demog$dem_bin)
    cont<-country_behav$connectivity*cont/connect
  }
  
  L<-list("cont"=cont,
          "country"=paste(country_behav$country, country_demog$country, sep="_"),
          "k"=country_behav$k,
          "dem_bin"=country_demog$dem_bin,
          "med_age"=country_demog$med_age)
  return(L)
}

# Return list of all combinations of behaviour and demography of the list
# 'countries'
list_of_mixes<-function(countries, norm)
{
  L<-list()
  i<-1
  for(country_behav in countries)
  {
    for(country_demog in countries)
    {
      L[[i]]<-mix(country_behav, country_demog, norm)
      i<-i+1
    }
  }
  return(L)
}