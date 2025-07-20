library(glue)
library(dplyr)
library(tidyr)
library(fst)
library(parallel)
library(foreach)
library(doParallel)

#functions----------------------------------------------------------------------

check29Feb<-function(dir)
{
  f<-list.files(dir,full.names = T)[1]
  ch29 <- system(paste0("cdo showdate ", f, " | grep -c '02-29'"), intern = TRUE)
  has29<-ifelse(ch29=="0",F,T)
  return(has29)
}

correct_csvFile<-function(file.name)
{
  df     <-read.csv(file.name,header = F)
  df[,2] <- paste(df[,3],df[,2],sep = ",")
  df     <- df[,-3]
  colnames(df)<-c("date","coords","values")
  df     <- df %>% pivot_wider(names_from = date,values_from = values) %>% as.data.frame()
  
  coords <- df[,1]
  coord  <<-strsplit(coords,",") %>% lapply(.,as.numeric) %>% do.call(rbind,.) 
  df     <- df[,-1]
  df     <- t(df)
  
  return(df)
}

#start parallel-----------------------------------------------------------------

my.cluster <- parallel::makeCluster(30,type = "PSOCK")
print(my.cluster)
doParallel::registerDoParallel(cl = my.cluster)

#inputs-------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
variable   <- args[1]
varType    <- args[2]

#pre-processing ERA5------------------------------------------------------------

  dir.era5<-glue("/scratch/aniruddha_s_hy.iitr/era5/data/{variable}") 
  var.era5<-glue("era5_day_{variable}")
  system(glue("cdo remapbil,/scratch/aniruddha_s_hy.iitr/era5/grid.txt {dir.era5}/{var.era5}.nc {dir.era5}/{var.era5}_1.nc")) #regridding to 1 degree
   
  system(glue("cdo ydaymean {dir.era5}/{var.era5}.nc {dir.era5}/{var.era5}_average.nc"))
  system(glue("cdo ydaymean {dir.era5}/{var.era5}_1.nc {dir.era5}/{var.era5}_1_average.nc")) #calculating daily climatology (averages)
  {
    if(varType == "additive")system(glue("cdo sub {dir.era5}/{var.era5}_1.nc {dir.era5}/{var.era5}_1_average.nc {dir.era5}/{var.era5}_1_anomaly.nc")) #calculating anomaly (additive)
    else                     system(glue("cdo div {dir.era5}/{var.era5}_1.nc {dir.era5}/{var.era5}_1_average.nc {dir.era5}/{var.era5}_1_anomaly.nc")) #calculating anomaly (multiplicative)
  }
  system(glue("cdo -outputtab,date,lat:6,lon:6,value:8 {dir.era5}/{var.era5}_1_anomaly.nc | grep -v '#' | tr -s ' ' | sed -e 's/ /,/g;s/^.//;s/.$//' >> {dir.era5}/{var.era5}_1_anomaly.csv"))
  era5.ano <- correct_csvFile(glue("{dir.era5}/{var.era5}_1_anomaly.csv"))  #preparing .csv into our required format
  era5.ano <- as.data.frame(era5.ano)
  write_fst(era5.ano, sprintf("%s/%s_1_anomaly.fst",dir.era5,var.era5),compress = 0,uniform_encoding = TRUE) #saving as .fst
  
#pre-processing CMIP6-----------------------------------------------------------

model.name<- c("ACCESS_CM2","ACCESS_ESM1_5","CanESM5","CESM2_WACCM","CMCC_CM2_SR5","CMCC_ESM2","FGOALS_g3","GFDL_ESM4","INM_CM4_8","INM_CM5_0","MRI_ESM2_0","MIROC6") 
cmip6.scenarios<- c("historical","ssp126","ssp245","ssp370","ssp585")

combo <- expand.grid(model.name, cmip6.scenarios)
colnames(combo)<-c("model","scenario")

foreach(i = 1:nrow(combo))%dopar%
{
  library(dplyr)
  library(tidyr)
  library(fst)
  
  model<-combo$model[i]
  dir.create(paste0("/scratch/aniruddha_s_hy.iitr/cmip6/repo/", model))
  
    scenario<-combo$scenario[i]
    
    dir.create(paste0("/scratch/aniruddha_s_hy.iitr/cmip6/repo/", model, "/", scenario))
    dir.create(paste0("/scratch/aniruddha_s_hy.iitr/cmip6/repo/", model, "/", scenario, "/", variable))
    
    dir.cmip6.raw <- paste0("/scratch/aniruddha_s_hy.iitr/cmip6/global/", model, "/", scenario, "/", variable)
    dir.cmip6     <- paste0("/scratch/aniruddha_s_hy.iitr/cmip6/repo/", model, "/", scenario, "/", variable)
    has29         <- check29Feb(dir.cmip6.raw)
    
    {#merging input files and converting into 1 degree
    if (scenario == "historical") system(paste0("cdo -remapbil,", dir.era5, "/", var.era5, "_1.nc -seldate,1960-01-01,2014-12-31 -mergetime ", dir.cmip6.raw, "/*.nc ", dir.cmip6, "/", variable, "_with29F.nc")) 
    else                          system(paste0("cdo -remapbil,", dir.era5, "/", var.era5, "_1.nc -seldate,2015-01-01,2100-12-31 -mergetime ", dir.cmip6.raw, "/*.nc ", dir.cmip6, "/", variable, "_with29F.nc"))
    }
     
    {#removing 29th Feb  
    if (has29) system(paste0("cdo -delete,dom=29feb ", dir.cmip6, "/", variable, "_with29F.nc ", dir.cmip6, "/", variable, "_1.nc")) 
    else       file.rename(from = paste0(dir.cmip6, "/", variable, "_with29F.nc"), to = paste0(dir.cmip6, "/", variable, "_1.nc"))
    }
    
    {
    system(paste0("cdo ydaymean ", dir.cmip6, "/", variable, "_1.nc ", dir.cmip6, "/", variable, "_1_average.nc"))#calculating daily climatology (averages)
    if (varType == "additive") system(paste0("cdo sub ", dir.cmip6, "/", variable, "_1.nc ", dir.cmip6, "/", variable, "_1_average.nc ", dir.cmip6, "/", variable, "_1_anomaly.nc")) #calculating anomaly (additive)
    else                       system(paste0("cdo div ", dir.cmip6, "/", variable, "_1.nc ", dir.cmip6, "/", variable, "_1_average.nc ", dir.cmip6, "/", variable, "_1_anomaly.nc")) #calculating anomaly (multiplicative)
    }
    
    system(paste0("cdo -outputtab,date,lat:6,lon:6,value:8 ", dir.cmip6, "/", variable, "_1_anomaly.nc | grep -v '#' | tr -s ' ' | sed -e 's/ /,/g;s/^.//;s/.$//' >> ", dir.cmip6, "/", variable, "_1_anomaly.csv")) #.nc to .csv
    
    ano<- correct_csvFile(paste0(dir.cmip6, "/", variable, "_1_anomaly.csv")) #preparing .csv into our required format
    ano<-as.data.frame(ano)
    
    file.remove(paste0(dir.cmip6, "/", variable, "_1_anomaly.nc"))
    file.remove(paste0(dir.cmip6,"/",variable,"_with29F.nc"))
    file.remove(paste0(dir.cmip6,"/",variable,"_1_anomaly.csv"))
    write_fst(ano, sprintf("%s/%s_1_anomaly.fst",dir.cmip6,variable),compress = 0,uniform_encoding = TRUE) #saving as .fst
    
    NULL
}


#end parallel-------------------------------------------------------------------
foreach::registerDoSEQ()
parallel::stopCluster(cl = my.cluster) 


