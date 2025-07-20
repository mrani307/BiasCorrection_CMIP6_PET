library(pbdMPI)

init()


suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(magrittr)))
suppressWarnings(suppressMessages(library(fst)))
suppressWarnings(suppressMessages(library(ncdf4)))

#functions----------------------------------------------------------------------
get_surrounding_doy <- function(day_of_year) {
  base_date <- as.Date("2002-01-01")  # Reference start of the year
  # Convert DOY to actual date
  date <- base_date + (day_of_year - 1)
  # Generate sequence from 5 days before to 5 days after
  seq_dates <- seq(date - 15, date + 15, by="day")
  # Convert back to DOY
  doy_values <- as.numeric(format(seq_dates, "%j")) # %j gives day of year
  return(doy_values)
}

doY_to_date<-function(doY){
  date<-as.Date(doY-1,origin = "2001-01-01")
  return(paste0("-",date %>% format(.,"%m"),"-",date %>% format(.,"%d")))
}

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

biasCorrection_historical<-function(grid){
  his.ano.grid <- read.fst(paste0(dir.cmip6,"/",variable,"_1_anomaly.fst"),columns = paste0("V",grid) )
  mod.ano.grid <- as.data.frame( read.fst(paste0(dir.cmip6,"/",variable,"_1_anomaly.fst"),columns = paste0("V",grid) ) )
  mod.ano.grid$q <- NA
  mod.ano.grid$A <- NA
  
  rownames(his.ano.grid) <- dates.his
  rownames(mod.ano.grid) <- dates.mod
  
  bc.grid<-c()
  
  for(doY in 1:365)
  {
    sliding.window<-get_surrounding_doy(doY)
    dates<-rownames(mod.ano.grid)
    mod.ano.grid.doY <- mod.ano.grid[which(dates %in% paste0(1960:2014,doY_to_date(doY))),]
    rep(1:365,55) %>% is_in(sliding.window) %>% which() %>% his.ano.grid[.,] -> his.ano.grid.swData
    
    q<-lapply(1960:2014,function(year) ecdf(his.ano.grid.swData)(mod.ano.grid.doY[paste0(year,doY_to_date(doY)),1])) %>% unlist()
    q[q<0.006]<-0.006
    q[q>0.994]<-0.994
    mod.ano.grid.doY$q<-round(q,2)
    
    adj.df<- adjFac[[grid]][[doY]]
    
    mod.ano.grid.doY$A<-lapply(1:nrow(mod.ano.grid.doY),function(i) adj.df$A[mod.ano.grid.doY$q[i]*100]) %>% unlist()
    
    bc<-mod.ano.grid.doY[,1]+mod.ano.grid.doY$A
    names(bc)<-rownames(mod.ano.grid.doY)
    
    bc.grid<-c(bc.grid,bc)
  }
  
  bc.grid %>% names() %>% as.Date() %>% order(.,decreasing = F) %>% bc.grid[.] -> bc.grid
  return(bc.grid)
}

biasCorrection_future<-function(grid)     {
  library(magrittr)
  library(dplyr)
  library(fst)
  
  his.ano.grid <- read.fst(paste0(dir.cmip6.his,"/",variable,"_1_anomaly.fst"),columns = paste0("V",grid) )
  mod.ano.grid <- as.data.frame( read.fst(paste0(dir.cmip6,"/",variable,"_1_anomaly.fst"),columns = paste0("V",grid) ) )
  mod.ano.grid$q <- NA
  mod.ano.grid$A <- NA
  
  rownames(his.ano.grid) <- dates.his
  rownames(mod.ano.grid) <- dates.mod
  
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
    
    adj.df<- adjFac[[grid]][[doY]]
    
    mod.ano.grid.doY$A<-lapply(1:nrow(mod.ano.grid.doY),function(i) adj.df$A[mod.ano.grid.doY$q[i]*100]) %>% unlist()
    
    if(varType=="additive") bc<-mod.ano.grid.doY[,1]+mod.ano.grid.doY$A
    else                    bc<-mod.ano.grid.doY[,1]*mod.ano.grid.doY$A

    names(bc)<-rownames(mod.ano.grid.doY)
    
    bc.grid<-c(bc.grid,bc)
  }
  
  bc.grid %>% names() %>% as.Date() %>% order(.,decreasing = F) %>% bc.grid[.] -> bc.grid
  bc.grid
}

#inputs-------------------------------------------------------------------------
  args <- commandArgs(trailingOnly = TRUE)
  model.name <- args[1]
  variable   <- args[2]  
  varUnit    <- args[3]
  varLonName <- args[4]
  varType    <- args[5]

  adjFac    <- readRDS(paste0("/scratch/aniruddha_s_hy.iitr/cmip6/repo/",model.name,"/historical/",variable,"/adjFac.RDS"))
  
#parallel setup-----------------------------------------------------------------  
  
  rank <- comm.rank()
  size <- comm.size()
  chunks   <- split(1:1440, cut(seq_along(1:1440), size, labels = FALSE))
  my_chunk <- chunks[[rank + 1]]  # +1 because R is 1-based
  
#historical scenario------------------------------------------------------------

  scenario   <- "historical"
  dir.cmip6 <- paste0("/scratch/aniruddha_s_hy.iitr/cmip6/repo/",model.name,"/",scenario,"/",variable)
  
  dates.his<-as.character( generate_dates(start_date = "1960-01-01", end_date = "2014-12-31") )
  dates.mod<-as.character( generate_dates(start_date = "1960-01-01", end_date = "2014-12-31") )

  date <- as.Date(dates.his)
  date <- as.numeric(date-as.Date("1850-01-01"))
  
  result_list   <- lapply(my_chunk, biasCorrection_historical) %>% do.call(rbind,.)
  gathered_data <- gather(result_list, rank.dest = 0)
  
  if (rank == 0) 
  {
    gathered_data <- do.call(rbind,gathered_data)
    saveRDS(gathered_data,paste0(dir.cmip6,"/",variable,"_1_anomaly_bc.RDS"))
  }
    
#future scenarios (ssp126)------------------------------------------------------

  scenario   <- "ssp126"
  
  dir.cmip6     <- paste0("/scratch/aniruddha_s_hy.iitr/cmip6/repo/",model.name,"/",scenario,"/",variable)
  dir.cmip6.his <- paste0("/scratch/aniruddha_s_hy.iitr/cmip6/repo/",model.name,"/historical/",variable)
  
  dates.his<-as.character( generate_dates(start_date = "1960-01-01", end_date = "2014-12-31") )
  dates.mod<-as.character( generate_dates(start_date = "2015-01-01", end_date = "2100-12-31") )
    
  date <- as.Date(dates.mod)
  date <- as.numeric(date-as.Date("1850-01-01"))
  
  result_list   <- lapply(my_chunk, biasCorrection_future) %>% do.call(rbind,.)
  gathered_data <- gather(result_list, rank.dest = 0)
  
  if (rank == 0) 
  {
    gathered_data <- do.call(rbind,gathered_data)
    saveRDS(gathered_data,paste0(dir.cmip6,"/",variable,"_1_anomaly_bc.RDS"))
  }  
  
#future scenarios (ssp245)------------------------------------------------------
  
  scenario   <- "ssp245"
  
  dir.cmip6     <- paste0("/scratch/aniruddha_s_hy.iitr/cmip6/repo/",model.name,"/",scenario,"/",variable)
  dir.cmip6.his <- paste0("/scratch/aniruddha_s_hy.iitr/cmip6/repo/",model.name,"/historical/",variable)
  
  dates.his<-as.character( generate_dates(start_date = "1960-01-01", end_date = "2014-12-31") )
  dates.mod<-as.character( generate_dates(start_date = "2015-01-01", end_date = "2100-12-31") )
  
  date <- as.Date(dates.mod)
  date <- as.numeric(date-as.Date("1850-01-01"))
  
  result_list   <- lapply(my_chunk, biasCorrection_future) %>% do.call(rbind,.)
  gathered_data <- gather(result_list, rank.dest = 0)
  
  if (rank == 0) 
  {
    gathered_data <- do.call(rbind,gathered_data)
    saveRDS(gathered_data,paste0(dir.cmip6,"/",variable,"_1_anomaly_bc.RDS"))
  }  
  
#future scenarios (ssp370)------------------------------------------------------
  
  scenario   <- "ssp370"
  
  dir.cmip6     <- paste0("/scratch/aniruddha_s_hy.iitr/cmip6/repo/",model.name,"/",scenario,"/",variable)
  dir.cmip6.his <- paste0("/scratch/aniruddha_s_hy.iitr/cmip6/repo/",model.name,"/historical/",variable)
  
  dates.his<-as.character( generate_dates(start_date = "1960-01-01", end_date = "2014-12-31") )
  dates.mod<-as.character( generate_dates(start_date = "2015-01-01", end_date = "2100-12-31") )
  
  date <- as.Date(dates.mod)
  date <- as.numeric(date-as.Date("1850-01-01"))
  
  result_list   <- lapply(my_chunk, biasCorrection_future) %>% do.call(rbind,.)
  gathered_data <- gather(result_list, rank.dest = 0)
  
  if (rank == 0) 
  {
    gathered_data <- do.call(rbind,gathered_data)
    saveRDS(gathered_data,paste0(dir.cmip6,"/",variable,"_1_anomaly_bc.RDS"))
  }  
  
#future scenarios (ssp585)------------------------------------------------------
  
  scenario   <- "ssp585"
  
  dir.cmip6     <- paste0("/scratch/aniruddha_s_hy.iitr/cmip6/repo/",model.name,"/",scenario,"/",variable)
  dir.cmip6.his <- paste0("/scratch/aniruddha_s_hy.iitr/cmip6/repo/",model.name,"/historical/",variable)
  
  dates.his<-as.character( generate_dates(start_date = "1960-01-01", end_date = "2014-12-31") )
  dates.mod<-as.character( generate_dates(start_date = "2015-01-01", end_date = "2100-12-31") )
  
  date <- as.Date(dates.mod)
  date <- as.numeric(date-as.Date("1850-01-01"))
  
  result_list   <- lapply(my_chunk, biasCorrection_future) %>% do.call(rbind,.)
  gathered_data <- gather(result_list, rank.dest = 0)
  
  if (rank == 0) 
  {
    gathered_data <- do.call(rbind,gathered_data)
    saveRDS(gathered_data,paste0(dir.cmip6,"/",variable,"_1_anomaly_bc.RDS"))
  }  
    
finalize()
