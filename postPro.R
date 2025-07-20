library(parallel)
library(foreach)
library(doParallel)


#functions----------------------------------------------------------------------
prepare_nc_file<-function(ncfname,date,fillvalue,name,
                          longname,units,nc.data,model,scenario   ){
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

generate_dates <- function(start_date, end_date ) {
  all_dates <- seq(as.Date(start_date), as.Date(end_date), by = "day")
  filtered_dates <- all_dates[!(format(all_dates, "%m-%d") == "02-29")]
  return(filtered_dates)
}


args <- commandArgs(trailingOnly = TRUE)
variable   <- args[1]
varUnit    <- args[2]
varLonName <- args[3]
varType    <- args[4]

model.name<- c("ACCESS_CM2","ACCESS_ESM1_5","CanESM5","CESM2_WACCM","CMCC_CM2_SR5","CMCC_ESM2","FGOALS_g3","GFDL_ESM4","INM_CM4_8","INM_CM5_0","MRI_ESM2_0","MIROC6") 
cmip6.scenarios<- c("historical","ssp126","ssp245","ssp370","ssp585")
combo <- expand.grid(model.name, cmip6.scenarios)
colnames(combo)<-c("model","scenario")

#start parallel-----------------------------------------------------------------

my.cluster <- parallel::makeCluster(30,type = "PSOCK")
print(my.cluster)
doParallel::registerDoParallel(cl = my.cluster)

#postPro-------------------------------------------------------------------------------

foreach(i = 1:nrow(combo)) %dopar%
{
  suppressWarnings(suppressMessages(library(ncdf4)))
  
  model   <-combo$model[i]
  scenario<-combo$scenario[i]
  
  if(scenario=="historical") date<-generate_dates(start_date = "1960-01-01", end_date = "2014-12-31")
  else                       date<-generate_dates(start_date = "2015-01-01", end_date = "2100-12-31")
  
  date <- as.numeric(date-as.Date("1850-01-01"))
  
  dir.cmip6<-paste0("/scratch/aniruddha_s_hy.iitr/cmip6/repo/", model, "/", scenario, "/", variable)
  dir.era5 <-paste0("/scratch/aniruddha_s_hy.iitr/era5/data/",variable)
  
  #preparing bias corrected anomaly nc file
  anomaly.bc<-readRDS(paste0(dir.cmip6,"/",variable,"_1_anomaly_bc.RDS"))
  prepare_nc_file(ncfname  =paste0(dir.cmip6,"/",variable,"_1_anomaly_bc.nc"),
                  date     =date,
                  fillvalue=-9999,
                  name     =variable,
                  longname =varLonName,
                  units    =varUnit,
                  nc.data  =anomaly.bc,
                  model    =model,
                  scenario =scenario)
  
  if(varType == "additive") 
  {
    system(paste0("cdo sub ", dir.era5, "/era5_day_", variable, "_1_average.nc ", dir.cmip6, "/", variable, "_1_average.nc ", dir.cmip6, "/c_factor.nc"))
    system(paste0("cdo -add ", dir.cmip6, "/", variable, "_1_average.nc ", dir.cmip6, "/c_factor.nc ", dir.cmip6, "/", variable, "_1_bc_step1.nc"))
    system(paste0("cdo add ", dir.cmip6, "/", variable, "_1_anomaly_bc.nc ", dir.cmip6, "/", variable, "_1_bc_step1.nc ", dir.cmip6, "/", variable, "_1_bc.nc"))
    file.remove(paste0(dir.cmip6, "/", variable, "_1_bc_step1.nc"))
    system(paste0("cdo -ydayadd -remapbil,", dir.era5, "/era5_day_", variable, "_average.nc -ydaysub ", dir.cmip6, "/", variable, "_1_bc.nc ", dir.era5, "/era5_day_", variable, "_1_average.nc ", dir.era5, "/era5_day_", variable, "_average.nc ", dir.cmip6, "/", variable, "_bc.nc"))
  } 
  else 
  {
    system(paste0("cdo div ", dir.era5, "/era5_day_", variable, "_1_average.nc ", dir.cmip6, "/", variable, "_1_average.nc ", dir.cmip6, "/c_factor.nc"))
    system(paste0("cdo -mul ", dir.cmip6, "/", variable, "_1_average.nc ", dir.cmip6, "/c_factor.nc ", dir.cmip6, "/", variable, "_1_bc_step1.nc"))
    system(paste0("cdo mul ", dir.cmip6, "/", variable, "_1_anomaly_bc.nc ", dir.cmip6, "/", variable, "_1_bc_step1.nc ", dir.cmip6, "/", variable, "_1_bc.nc"))
    file.remove(paste0(dir.cmip6, "/", variable, "_1_bc_step1.nc"))
    system(paste0("cdo -ydaymul -remapbil,", dir.era5, "/era5_day_", variable, "_average.nc -ydaydiv ", dir.cmip6, "/", variable, "_1_bc.nc ", dir.era5, "/era5_day_", variable, "_1_average.nc ", dir.era5, "/era5_day_", variable, "_average.nc ", dir.cmip6, "/", variable, "_bc.nc"))
  }
  
  NULL
}

#end parallel-------------------------------------------------------------------
foreach::registerDoSEQ()
parallel::stopCluster(cl = my.cluster) 

  
  
  
  
  
  
  
  