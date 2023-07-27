################################################################################
# Functions to pre process contact matrices
################################################################################

#Correct for reciprocity
reciprocity<-function(dem, cont)
{
  age_cats<-rownames(cont)
  dem_bin<-bin_demo(dem, age_cats)
  cont2<-cont
  
  n<-nrow(dem_bin)
  
  for(i in 1:n)
  {
    for(j in 1:n)
    {
      cont2[i,j]<-0.5*(cont[i,j]+(cont[j,i]*dem_bin$Value[i]/dem_bin$Value[j]) )
    }
  }
  
  return(cont2)
}

# Takes demography distribution and caps it at a certain max age
cap_demo<-function(dem, max_age)
{
  n<-nrow(dem)
  
  # ages 0,1,2,3,4,....max_age-1
  dem_cap<-dem[1:max_age,]
  
  # oldest category
  Age<-paste0(max_age, "+")
  Value<-sum(dem[(max_age+1):n,]$Value)
  
  dem_cap<-rbind(dem_cap, data.frame(Age, Value))
  rownames(dem_cap)<-dem_cap$Age
  return(dem_cap)
}

# Bin demographic data to age_cats
bin_demo<-function(dem, age_cats)
{
  n_age<-length(age_cats)
  n<-nrow(dem)
  
  dem_binned<-NULL
  for(i in 1:(n_age))
  {
    if(!grepl("+", age_cats[i], fixed=TRUE))
    {
      min_age<-as.numeric(gsub("-.*", "", age_cats[i]))
      max_age<-as.numeric(gsub(".*-", "", age_cats[i]))
      Value<-sum(dem[(min_age+1):(max_age+1),]$Value)
    }
    else
    {
      min_age<-as.numeric(gsub("\\+.*", "", age_cats[i]))
      Value<-sum(dem[(min_age+1):n,]$Value)
    }
    
    Age<-age_cats[[i]]
    dem_binned<-rbind(dem_binned, data.frame(Age, Value))
  }
  rownames(dem_binned)<-age_cats
  
  return(dem_binned)
}

#Form numeric sequence from min to max.
#if input is "min-max" output is min, min+1, min+2 ..... max
#if input is "max+" output is max
age_seq<-function(age_cat)
{
  if(grepl("+", age_cat, fixed=TRUE))
  {
    return(c(as.numeric(gsub("\\+.*", "", age_cat))) )
  }
  else
  {
    lo<-as.numeric(gsub("-.*","", age_cat))
    hi<-as.numeric(gsub(".*-","", age_cat))
    return(seq(from=lo, to=hi, by=1))
  }
}

# In a bin of width ... do yearly breakdown
find_sub_conts<-function(age_cat1, age_cat2, dem_cap, cont)
{
  # Age_cat1 is contacts, age_cat 2 is participant
  # Contacts between categories
  c_12<-cont[age_cat1, age_cat2]
  
  ages1_i<-age_seq(age_cat1)+1
  ages2_i<-age_seq(age_cat2)+1
  
  N<-length(ages1_i)
  M<-length(ages2_i)
  
  m<-matrix(c_12, nrow=N, ncol=M)
  
  # prob distribution of ages of contacts
  sub_dem<-dem_cap$Value[ages1_i]
  sub_dem<-sub_dem/sum(sub_dem)
  
  # Spread contacts out between these
  m<-m*sub_dem
  
  return(m)
}

# Build matrix of contacts, where age categories are by year.
build_yearly_matrix<-function(cont, dem_raw)
{
  age_cats<-colnames(cont)
  n<-length(age_cats)
  max_age<-as.numeric(gsub("\\+.*", "", age_cats[n]))
  
  # Get age bins, ages="0,1,2,3...max-1, max+"
  ages<-seq(0, (max_age-1), 1)
  ages<-c(ages, age_cats[n])
  N<-length(ages)
  
  # Empty matrix to receive values
  m<-matrix(0, nrow = N, ncol = N)
  rownames(m)<-ages
  colnames(m)<-ages
  
  # Cap demography at last age category, so that it has same bins as m
  dem_cap<-cap_demo(dem_raw, max_age)
  
  # Break down contact matrix into yearly resolution
  for(i in age_cats)
  {
    ages1<-age_seq(i)+1
    for(j in age_cats)
    {
      ages2<-age_seq(j)+1
      m[ages1, ages2]<-find_sub_conts(i, j, dem_cap, cont)
    }
  }
  
  return(m)
}

# If max category is e.g. 50+, modify to achieve e.g. 70+
increase_max<-function(cont, dem_raw, new_max)
{
  age_cats<-colnames(cont)
  n<-length(age_cats)
  pre_max<-as.numeric(gsub("\\+.*", "", age_cats[n]))
  
  ######################
  # Add new participants
  ######################
  # Add extra column, if e.g. 50+ behave in a certain way then we assume that
  # 50-69 and 70+ behave the same too
  cont<-cbind(cont, cont[,n, drop = F])
  
  ##################
  # Add new contacts
  ##################
  N<-nrow(dem_raw)
  # Number of people between pre_max-max
  N_lo<-sum(dem_raw$Value[(pre_max+1):(new_max)])
  # Number of people max+
  N_hi<-sum(dem_raw$Value[(new_max+1):N])
  
  # Fraction to split between rows
  tot<-N_lo+N_hi
  f_lo<-N_lo/tot
  f_hi<-N_hi/tot
  
  # Split and add to matrix
  row_to_split<-cont[n, , drop = F]
  row1<-row_to_split*f_lo
  row2<-row_to_split*f_hi
  cont[n,]<-row1
  cont<-rbind(cont, row2)
  
  ######################
  # Rename cols and rows
  ######################
  cat1<-paste(pre_max, new_max-1, sep="-")
  cat2<-paste0(new_max, sep="+")
  
  rownames(cont)[n]<-cat1
  rownames(cont)[n+1]<-cat2
  colnames(cont)[n]<-cat1
  colnames(cont)[n+1]<-cat2
  
  return(cont)
}


# Reduce max if e.g. 80+ and we want 70+ (works on yearly matrix)
reduce_max<-function(cont, dem_raw, maximum)
{
  age_cats<-colnames(cont)
  n<-length(age_cats)
  pre_max<-as.numeric(gsub("\\+.*", "", age_cats[n]))
  dem_cap<-cap_demo(dem_raw, pre_max)
  
  i_old<-pre_max+1
  i_max<-maximum+1
  
  #####################################################
  # Add up contacts in new maximum+category for new row
  #####################################################
  x<-colSums(cont[i_max:i_old,])
  # replace age rows with new max+ row
  cont<-cont[-i_max:-i_old,]
  cont<-rbind(cont,x)
  
  ###############################
  # Combine last columns into one
  ###############################
  
  # Age distribution for ages between i_max-i_old
  dem_sub<-dem_cap$Value[i_max:i_old]
  dem_sub<-dem_sub/sum(dem_sub)
  dem_sub<-as.matrix(dem_sub)
  
  # Sub matrix for last columns
  m_sub<-cont[,i_max:i_old]
  
  # Take weighted sum to get average contacts
  last_col<-m_sub %*% dem_sub
  
  # Replace last column
  cont<- cont[,-i_max:-i_old]
  cont<- cbind(cont,last_col)
  
  colnames(cont)[i_max]<-paste0(maximum, "+")
  rownames(cont)[i_max]<-paste0(maximum, "+")
  
  return(cont)
}


# Sort yearly matrix into bins
bin_matrix<-function(cont_yr, bins, dem_raw)
{
  N<-length(bins)
  max_age<-as.numeric(gsub("\\+.*", "", bins[N]))
  dem_cap<-cap_demo(dem_raw, max_age)
  
  # Initialize empty matrix
  bin_mat<-matrix(0, nrow = N, ncol = N)
  rownames(bin_mat)<-bins
  colnames(bin_mat)<-bins
  
  #Do things by index, not row/col name. Since cont_yr goes 0, 1, 2, 3... the indices are shifted by 1.
  for(i in bins[1:N])
  {
    ages1<-age_seq(i)+1
    
    for(j in bins[1:N])
    {
      ages2<-age_seq(j)+1
      
      subm<-cont_yr[ages1, ages2, drop = F]
      
      subm<-colSums(subm)
      
      sub_dem<-dem_cap$Value[ages2]
      sub_dem<-sub_dem/sum(sub_dem)
      
      bin_mat[i,j]<-as.numeric(sub_dem %*% subm)
    }
  }
  return(bin_mat)
}


# Process matrix to have desired max age and binwidth
process_matrix<-function(dem_raw, cont, max_age, bin_width)
{
  #######################################
  # (1) Set max age bin of contact matrix
  #######################################
  
  # Find current max age bin
  age_cats<-colnames(cont)
  n<-length(age_cats)
  old_max<-as.numeric(gsub("\\+.*", "", age_cats[n]))
  
  # If not=max_age, modify the matrix
  if(old_max<max_age)
  {
    cont<-increase_max(cont, dem_raw, max_age)
    age_cats<-colnames(cont)
    n<-length(age_cats)
  }
  
  ##################################
  # (2) Form a yearly contact matrix
  ##################################
  
  cont_yr<-build_yearly_matrix(cont, dem_raw)
  if(old_max>max_age)
  {
    cont_yr<-reduce_max(cont_yr, dem_raw, max_age)
  }
  
  #######################
  # (3) Form new age bins
  #######################
  lower_bin<-seq(from=0, to=max_age-bin_width, by=bin_width)
  upper_bin<-seq(from=bin_width-1, to=max_age-1, by=bin_width)
  bins<-paste(lower_bin, upper_bin, sep="-")
  bins<-c(bins, paste0(max_age, "+"))
  
  ###################################
  # (4) Bin yearly data into new bins
  ###################################
  cont_bin<-bin_matrix(cont_yr, bins, dem_raw)
  
  return(cont_bin)
}



################################################################################
# Functions to get metrics from matrices
################################################################################

# Connectivity as defined by Arregui
connectivity<-function(cont, dem)
{
  dem<-dem$Value
  N<-length(dem)
  
  out<-0
  for(i in 1:N)
  {
    for(j in 1:N)
    {
      out<-out+cont[i,j]*dem[i]
    }
  }
  
  out<-out/sum(dem)
  
  return(out)
}

# Get random mix
rand_mix<-function(dem)
{
  x<-dem$Value
  x<-x/sum(x)
  N<-length(x)
  
  m<-matrix(x,N,N)
  
  rownames(m)<-rownames(dem)
  colnames(m)<-rownames(dem)
  
  return(m)
}

# Get behavioural matrix
get_behav<-function(cont, dem)
{
  k<-connectivity(cont, dem)
  cont<-cont/k
  
  rand<-rand_mix(dem)
  
  beh<-cont/rand
  
  return(beh)
}

# Derive Y parameter
Y_fun<-function(behav)
{
  N<-nrow(behav)
  
  age_cats<-colnames(behav)
  upper<-as.numeric(gsub(".*-", "", age_cats))
  y_max<-max(which(upper<20))
  
  out<-0
  
  for(i in 1:y_max)
  {
    for(j in 1:N)
    {
      out<-out+behav[i,j]
    }
  }
  
  return(out/sum(behav))
}

median_age<-function(demo_raw)
{
  demo_raw$Age<-as.numeric(as.character(demo_raw$Age))
  
  demo_raw$cumsum<-cumsum(demo_raw$Value)
  tot<-sum(demo_raw$Value)
  demo_raw$mid_diff<-demo_raw$cumsum-0.5*tot
  i<-findInterval(0, demo_raw$mid_diff)
  demo_raw<-demo_raw[i:(i+1),]
  demo_raw$mid_diff<-abs(demo_raw$mid_diff)
  
  med_age<-demo_raw$Age[1]+(demo_raw$mid_diff[1]/(sum(demo_raw$mid_diff)))
  
  return(med_age)
}

################################################################################
# Put it all together
################################################################################


process_country<-function(country_dat, maximum, bin_width)
{
  country<-country_dat$country
  cont_raw<-country_dat$cont_raw
  dem<-country_dat$dem
  
  # Adjust for reciprocity
  cont_raw<-reciprocity(dem, cont_raw)
  
  cont<-process_matrix(dem, cont_raw, maximum, bin_width)
  dem_bin<-bin_demo(dem, colnames(cont))
  
  #Normalize such that R0=2.5 for a fully susceptible population
  e<-eigen(cont)
  R0<-(max(Re(e$values)))
  scale<-2.5/R0
  cont<-scale*cont
  
  connect<-connectivity(cont, dem_bin)
  behav<-get_behav(cont, dem_bin)
  
  avg_age<-median_age(dem)
  k<-Y_fun(behav)
  rand<-rand_mix(dem_bin)
  
  return(list("country"=country,
              "cont_raw"=cont_raw,
              "cont"=cont,
              "connectivity"=connect,
              "dem"=dem,
              "dem_bin"=dem_bin,
              "behav"=behav,
              "rand"=rand,
              "med_age"=avg_age,
              "k"=k
  ) )
}