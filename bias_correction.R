library(dplyr)
library(tidyr)
library(glue)
library(magrittr)
library(raster)
library(ncdf4)
library(withr)

library(parallel)
library(foreach)
library(doParallel)

date.fut      <-data.frame(ts=seq(as.Date("2015-01-01"),as.Date("2100-12-31"),1))
date.fut$year <-as.numeric(format(date.fut$ts,"%Y"))                 
date.fut$month<-as.numeric(format(date.fut$ts,"%m"))                 
date.fut$day  <-as.numeric(format(date.fut$ts,"%d"))                 
pos.fut       <-which(date.fut$month==2 & date.fut$day==29)

date.his      <-data.frame(ts=seq(as.Date("1960-01-01"),as.Date("2014-12-31"),1))
date.his$year <-as.numeric(format(date.his$ts,"%Y"))                 
date.his$month<-as.numeric(format(date.his$ts,"%m"))                 
date.his$day  <-as.numeric(format(date.his$ts,"%d"))                 
pos.his       <-which(date.his$month==2 & date.his$day==29)

#functions----------------------------------------------------------------------


calculate_run_mean<-function(grid,pad=4)
{
  n.avg    <- pad*2+1
  grid.len <- length(grid)
  run.mean <- c()
  for(t in 1:grid.len)
  {
    if( t <= pad)
    {
      if(t==1) c(grid[5:2],grid[1:5]) %>% mean() -> run.mean[1]
      if(t==2) c(grid[4:1],grid[2:6]) %>% mean() -> run.mean[2]
      if(t==3) c(grid[4:1],grid[3:7]) %>% mean() -> run.mean[3]
      if(t==4) c(grid[4:1],grid[4:8]) %>% mean() -> run.mean[4]
    }
    else if( t > (grid.len-pad))
    {
      if(t==(grid.len-3)) c(grid[(grid.len-7):(grid.len)],grid[(grid.len-4)])              %>% mean() -> run.mean[grid.len-3]
      if(t==(grid.len-2)) c(grid[(grid.len-6):(grid.len)],grid[(grid.len-4):(grid.len-3)]) %>% mean() -> run.mean[grid.len-2]
      if(t==(grid.len-1)) c(grid[(grid.len-5):(grid.len)],grid[(grid.len-4):(grid.len-2)]) %>% mean() -> run.mean[grid.len-1]
      if(t==grid.len)     c(grid[(grid.len-4):(grid.len)],grid[(grid.len-4):(grid.len-1)]) %>% mean() -> run.mean[grid.len]
      
    }
    else
    {
      run.mean[t] <- c(grid[(t-pad):(t+pad)]) %>% mean()
    }
  }
  return(run.mean)
}

prep_anomaly_data<-function(file.name)
{
  df     <-read.csv(file.name,header = F)
  df[,2] <- paste(df[,3],df[,2],sep = ",")
  df     <- df[,-3]
  colnames(df)<-c("date","coords","values")
  df     <- df %>% pivot_wider(names_from = date,values_from = values) %>% as.data.frame()
  
  coords <- df[,1]
  coord  <<-strsplit(coords,",") %>% lapply(.,as.numeric) %>% do.call(rbind,.) 
  df     <- df[,-1]
  
  return(df)
}

get_surrounding_doy <- function(day_of_year) 
{
  base_date <- as.Date("2002-01-01")  # Reference start of the year
  # Convert DOY to actual date
  date <- base_date + (day_of_year - 1)
  # Generate sequence from 5 days before to 5 days after
  seq_dates <- seq(date - 15, date + 15, by="day")
  # Convert back to DOY
  doy_values <- as.numeric(format(seq_dates, "%j")) # %j gives day of year
  return(doy_values)
}

doY_to_date<-function(doY)
{
  date<-as.Date(doY-1,origin = "2001-01-01")
  return(paste0("-",date %>% format(.,"%m"),"-",date %>% format(.,"%d")))
}

prepare_nc_file<-function(ncfname,date,fillvalue,name,
                          longname,units,nc.data,model,scenario   )
{
  # define dimensions
  londim  <- ncdim_def("lon","degrees_east", as.double(seq(59.5,98.5,1))) 
  latdim  <- ncdim_def("lat","degrees_north",as.double(seq( 4.5,39.5,1)))
  timedim <- ncdim_def("time","days since 1850-01-01",as.double(date))
  
  # define variables
  
  var_def   <- ncvar_def(name  = name,
                         units = units,
                         dim   = list(londim,latdim,timedim),
                         missval  = fillvalue,
                         longname = longname,
                         prec     = "double")
  
  # create netCDF file and put arrays
  ncout <- nc_create(ncfname,list(var_def),force_v4=TRUE)
  
  # put variables
  ncvar_put(ncout,var_def,nc.data)
  ncatt_put(ncout,"lon","axis","X") 
  ncatt_put(ncout,"lat","axis","Y")
  ncatt_put(ncout,"time","axis","T")
  
  # add global attributes
  ncatt_put(ncout,0,"title",paste0(model,"_",scenario ))
  ncatt_put(ncout,0,"institution","BC by : DoH IIT Roorkee")
  ncatt_put(ncout,0,"source","https://esgf-ui.ceda.ac.uk/cog/search/cmip6-ceda/")
  
  nc_close(ncout)
}

#input data --------------------------------------------------------------------  

direc.cmip<- "/scratch/aniruddha_s_hy.iitr/cmip6/repo"
direc.era5<- "/scratch/aniruddha_s_hy.iitr/era5/data"
variable  <- "rsds"
varName   <- "rsds"
varUnit   <- "Wm-2"
varLonName<- "Solar Radiation"

#reading observation dataset ---------------------------------------------------

print(glue("ERA5 {variable} : reading observed daily anomolies"))
obs.ano   <- prep_anomaly_data(file.name = glue("{direc.era5}/{variable}/era5_anomaly_{variable}_1.csv"))

#start parallel-----------------------------------------------------------------

my.cluster <- parallel::makeCluster(48,type = "PSOCK")
print(my.cluster)
doParallel::registerDoParallel(cl = my.cluster)

#historical---------------------------------------------------------------------

model.name<-c("ACCESS_CM2","ACCESS_ESM1_5","CanESM5","CESM2_WACCM","CMCC_CM2_SR5","CMCC_ESM2","FGOALS_g3","GFDL_ESM4","INM_CM4_8","INM_CM5_0","MIROC6","MRI_ESM2_0")

for(model in model.name[1:12])
{
  scenario<-"historical"  
  {
    source("/home/aniruddha_s_hy.iitr/src/preprocessing_cmip.R")
    
    print(glue("{variable} CMIP6 {model} historical: preparing daily anomolies......."))
    his.ano   <- prep_anomaly_data(file.name = glue("{direc.cmip}/{model}/{scenario}/{variable}/{variable}_anomaly_{model}_{scenario}_1.csv"))
    
    #calculating adjustment factor function
    print(glue("{variable} CMIP6 {model} : calculating adjustment factor function......."))
    adj.doY<-foreach(grid = 1:nrow(his.ano)) %dopar%
      {
        library(dplyr)
        library(magrittr)
        
        obs.ano.grid <- obs.ano[grid,] %>% t() #reading observational anomaly time series for a grid (1960-2014)
        his.ano.grid <- his.ano[grid,] %>% t() #reading historical    anomaly time series for a grid (1960-2014)
        
        adj.fac.doY<-list() #list to store the adjustment factor function (dataframe) for each doY for a grid
        
        for(doY in 1:365)
        {
          sliding.window<-get_surrounding_doy(doY) #extracting +-15days for a given doY to prepare a sliding window
          
          rep(1:365,55) %>% is_in(sliding.window) %>% which() %>% obs.ano.grid[.,] -> obs.ano.grid.swData #extracting anomalies for the sliding window for each year (31x55)
          rep(1:365,55) %>% is_in(sliding.window) %>% which() %>% his.ano.grid[.,] -> his.ano.grid.swData #extracting anomalies for the sliding window for each year (31x55)
          
          #computing adjustment factor function, by sampling q with 99 values starting from 0.01 to 0.99 by steps of 0.01
          adj.fac.doY[[doY]]<-data.frame(q=seq(0.01,0.99,0.01),
                                         A=lapply(seq(0.01,0.99,0.01), function(q) quantile(obs.ano.grid.swData,q)-quantile(his.ano.grid.swData,q))%>% unlist())
        }
        adj.fac.doY
      }
    
    #bias correction (quantile mapping)
    print(glue("{variable} CMIP6 {model} {scenario}: Bias correction starts .. {(Sys.time())} "))
    his.ano.bc<-foreach(grid = 1:nrow(his.ano),.combine='rbind') %dopar%
      {
        library(magrittr)
        library(dplyr)
        
        his.ano.grid <- his.ano[grid,] %>% t()
        mod.ano.grid <- his.ano[grid,] %>% t() %>% as.data.frame()
        mod.ano.grid$q <- NA
        mod.ano.grid$A <- NA
        
        bc.grid<-c()
        
        for(doY in 1:365)
        {
          sliding.window<-get_surrounding_doy(doY)
          
          dates<-mod.ano.grid %>% rownames()
          mod.ano.grid.doY<-mod.ano.grid[which(dates %in% paste0(1960:2014,doY_to_date(doY))),]
          rep(1:365,55) %>% is_in(sliding.window) %>% which() %>% his.ano.grid[.,] -> his.ano.grid.swData
          
          q<-lapply(1960:2014,function(year) ecdf(his.ano.grid.swData)(mod.ano.grid.doY[paste0(year,doY_to_date(doY)),1])) %>% unlist()
          q[q<0.006]<-0.006
          q[q>0.994]<-0.994
          mod.ano.grid.doY$q<-round(q,2)
          
          adj.df<- adj.doY[[grid]][[doY]]
          
          mod.ano.grid.doY$A<-lapply(1:nrow(mod.ano.grid.doY),function(i) adj.df$A[mod.ano.grid.doY$q[i]*100]) %>% unlist()
          
          bc<-mod.ano.grid.doY[,1]+mod.ano.grid.doY$A
          names(bc)<-rownames(mod.ano.grid.doY)
          
          bc.grid<-c(bc.grid,bc)
        }
        
        bc.grid %>% names() %>% as.Date() %>% order(.,decreasing = F) %>% bc.grid[.] -> bc.grid
        bc.grid
      }
    print(glue("{variable} CMIP6 {model} {scenario}: Bias correction ends .. {(Sys.time())} "))
    
    date <- his.ano.bc %>% colnames() %>% as.Date()
    date <- as.numeric(date-as.Date("1850-01-01"))
    prepare_nc_file(ncfname   = glue("{direc.cmip}/{model}/{scenario}/{variable}/{variable}_day_{model}_{scenario}_ano_bc.nc"),
                    date      = date,
                    fillvalue = -9999,
                    name      = varName,
                    longname  = varLonName,
                    units     = varUnit,
                    nc.data   = his.ano.bc,
                    model     = model,
                    scenario  = scenario)
    rm(his.ano.bc,date)
    print(glue("{variable} CMIP6 {model} {scenario}: Bias corrected anomaly data written as .nc file at 1 degree resolution........ "))
  }
  
  #bias correction of anomalies
  for(scenario in c("ssp126","ssp245","ssp370","ssp585"))
  {
    #pre-processing
    source("/home/aniruddha_s_hy.iitr/src/preprocessing_cmip.R")
    
    #reading future anomaly 
    print(glue("{variable} CMIP6 {model} {scenario}: reading daily anomolies......."))
    mod.ano   <- prep_anomaly_data(file.name = glue("{direc.cmip}/{model}/{scenario}/{variable}/{variable}_anomaly_{model}_{scenario}_1.csv"))
    
    #bias correction (quantile mapping)
    print(glue("{variable} CMIP6 {model} {scenario}: Bias correction starts .. {(Sys.time())} "))
    mod.ano.bc<-foreach(grid = 1:nrow(mod.ano),.combine='rbind') %dopar%
      {
        library(magrittr)
        library(dplyr)
        
        his.ano.grid <- his.ano[grid,] %>% t()
        mod.ano.grid <- mod.ano[grid,] %>% t() %>% as.data.frame()
        mod.ano.grid$q <- NA
        mod.ano.grid$A <- NA
        
        bc.grid<-c()
        
        for(doY in 1:365)
        {
          sliding.window<-get_surrounding_doy(doY)
          
          dates<-mod.ano.grid %>% rownames()
          mod.ano.grid.doY<-mod.ano.grid[which(dates %in% paste0(2015:2100,doY_to_date(doY))),]
          rep(1:365,55) %>% is_in(sliding.window) %>% which() %>% his.ano.grid[.,] -> his.ano.grid.swData
          
          q<-lapply(2015:2100,function(year) ecdf(his.ano.grid.swData)(mod.ano.grid.doY[paste0(year,doY_to_date(doY)),1])) %>% unlist()
          q[q<0.006]<-0.006
          q[q>0.994]<-0.994
          mod.ano.grid.doY$q<-round(q,2)
          
          adj.df<- adj.doY[[grid]][[doY]]
          
          mod.ano.grid.doY$A<-lapply(1:nrow(mod.ano.grid.doY),function(i) adj.df$A[mod.ano.grid.doY$q[i]*100]) %>% unlist()
          
          bc<-mod.ano.grid.doY[,1]+mod.ano.grid.doY$A
          names(bc)<-rownames(mod.ano.grid.doY)
          
          bc.grid<-c(bc.grid,bc)
        }
        
        bc.grid %>% names() %>% as.Date() %>% order(.,decreasing = F) %>% bc.grid[.] -> bc.grid
        bc.grid
      }
    rm(mod.ano)
    print(glue("{variable} CMIP6 {model} {scenario}: Bias correction ends .. {(Sys.time())} "))
    
    date <- mod.ano.bc %>% colnames() %>% as.Date()
    date <- as.numeric(date-as.Date("1850-01-01"))
    
    prepare_nc_file(ncfname   = glue("{direc.cmip}/{model}/{scenario}/{variable}/{variable}_day_{model}_{scenario}_ano_bc.nc"),
                    date      = date,
                    fillvalue = -9999,
                    name      = varName,
                    longname  = varLonName,
                    units     = varUnit,
                    nc.data   = mod.ano.bc,
                    model     = model,
                    scenario  = scenario)
    rm(mod.ano.bc,date)
    print(glue("{variable} CMIP6 {model} {scenario}: Bias corrected anomaly data written as .nc file at 1 degree resolution........ "))
  }
  
  rm(his.ano,adj.doY)
  
  #down scaling 
  for(scenario in c("historical","ssp126","ssp245","ssp370","ssp585"))
  {
    
    #preparing bias corrected data (at 1 degree)
    system(glue("cdo add {direc.cmip}/{model}/{scenario}/{variable}/{variable}_day_{model}_{scenario}_ano_bc.nc -add {direc.cmip}/{model}/{scenario}/{variable}/{variable}_mean_{model}_{scenario}_1.nc -timmean -sub {direc.era5}/{variable}/era5_mean_{variable}_1.nc {direc.cmip}/{model}/historical/{variable}/{variable}_mean_{model}_historical_1.nc {direc.cmip}/{model}/{scenario}/{variable}/{variable}_day_{model}_{scenario}_bc.nc"))
    
    #preparing bias corrected downscaled data
    system(glue("cdo -ydayadd -remapbil,{direc.era5}/{variable}/era5_dayClim_{variable}.nc -ydaysub {direc.cmip}/{model}/{scenario}/{variable}/{variable}_day_{model}_{scenario}_bc.nc {direc.era5}/{variable}/era5_dayClim_{variable}_1.nc {direc.era5}/{variable}/era5_dayClim_{variable}.nc {direc.cmip}/{model}/{scenario}/{variable}/{variable}_day_{model}_{scenario}_bcds_temp.nc"))
    
    #correcting values greater than 100 and lower than 0
    system(glue("cdo remapbil,{direc.era5}/{variable}/era5_dayClim_{variable}.nc -div {direc.cmip}/{model}/{scenario}/{variable}/{variable}_day_{model}_{scenario}_bcds_temp.nc -mul -gtc,0 {direc.cmip}/{model}/{scenario}/{variable}/{variable}_day_{model}_{scenario}_bcds_temp.nc -ltc,500 {direc.cmip}/{model}/{scenario}/{variable}/{variable}_day_{model}_{scenario}_bcds_temp.nc {direc.cmip}/{model}/{scenario}/{variable}/{variable}_day_{model}_{scenario}_bcds.nc"))
    
    print(glue("{variable} CMIP6 {model} {scenario}: Downscaling completed ..................................... "))
    
    #adding 29th Feb
    print(glue("{variable} CMIP6 {model} {scenario}: Adding 29th Feb ..................................... "))
    input.file <- glue("{direc.cmip}/{model}/{scenario}/{variable}/{variable}_day_{model}_{scenario}_bcds.nc")
    files<-""
    
    {
      if(scenario=="historical")
      {
        df  <- date.his
        pos <- pos.his
      }
      else
      {
        df  <- date.fut
        pos <- pos.fut
      }
    }
    
    for(y in df[pos,"year"])
    {
      system(glue("cdo seldate,{y}-02-28 {input.file} {direc.cmip}/{model}/{scenario}/{variable}/feb28_{y}.nc"))
      system(glue("cdo seldate,{y}-03-01 {input.file} {direc.cmip}/{model}/{scenario}/{variable}/mar1_{y}.nc"))
      
      system(glue("cdo -mergetime {direc.cmip}/{model}/{scenario}/{variable}/feb28_{y}.nc {direc.cmip}/{model}/{scenario}/{variable}/mar1_{y}.nc {direc.cmip}/{model}/{scenario}/{variable}/merged_{y}.nc"))
      system(glue("cdo timmean {direc.cmip}/{model}/{scenario}/{variable}/merged_{y}.nc {direc.cmip}/{model}/{scenario}/{variable}/feb29_{y}.nc"))
      
      file.remove(c(glue("{direc.cmip}/{model}/{scenario}/{variable}/feb28_{y}.nc"),
                    glue("{direc.cmip}/{model}/{scenario}/{variable}/mar1_{y}.nc"),
                    glue("{direc.cmip}/{model}/{scenario}/{variable}/merged_{y}.nc")))
      files<-paste0(files," ",glue("{direc.cmip}/{model}/{scenario}/{variable}/feb29_{y}.nc"))
    }
    system(glue("cdo mergetime {input.file} {files} {direc.cmip}/{model}/{scenario}/{variable}/{variable}_day_{model}_{scenario}_bcds_wtp.nc"))
    for(y in df[pos,"year"])
      file.remove(glue("{direc.cmip}/{model}/{scenario}/{variable}/feb29_{y}.nc"))
    
    file.remove(glue("{direc.cmip}/{model}/{scenario}/{variable}/{variable}_day_{model}_{scenario}_bcds.nc"))
    file.rename(from = glue("{direc.cmip}/{model}/{scenario}/{variable}/{variable}_day_{model}_{scenario}_bcds_wtp.nc"),
                to = glue("{direc.cmip}/{model}/{scenario}/{variable}/{variable}_day_{model}_{scenario}_bcds.nc"))
    
    print(glue("{variable} CMIP6 {model} {scenario}: / Bias correction is completed  {(Sys.time())}............................................................."))
  }
  
  print(glue("{variable} CMIP6 {model} : / Bias correction is completed  for the model ....................................................................."))
  print("                                                                                                                                                  ")
  print("--------------------------------------------------------------------------------------------------------------------------------------------------")
  print("                                                                                                                                                  ")
}

#end parallel
foreach::registerDoSEQ()
parallel::stopCluster(cl = my.cluster) 




