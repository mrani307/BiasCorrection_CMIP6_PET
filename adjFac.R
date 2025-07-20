library(pbdMPI)

init()

suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(magrittr)))
suppressWarnings(suppressMessages(library(fst)))

#reading inputs---------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
model.name <- args[1]  
variable   <- args[2]
varType    <- args[3]


#function---------------------------------------------------------------------

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

calc_adjFac<-function(grid)
{

  obs.ano.grid <- read.fst(paste0("/scratch/aniruddha_s_hy.iitr/era5/data/",variable,"/era5_day_",variable,"_1_anomaly.fst"),columns=paste0("V",grid))                   #reading observational anomaly time series for a grid (1960-2014)
  his.ano.grid <- read.fst(paste0("/scratch/aniruddha_s_hy.iitr/cmip6/repo/",model.name,"/historical/",variable,"/",variable,"_1_anomaly.fst"),columns=paste0("V",grid)) #reading historical    anomaly time series for a grid (1960-2014)
  
  adj.fac.doY<-list() #list to store the adjustment factor function (dataframe) for each doY for a grid
  
  for(doY in 1:365)
  {
    sliding.window<-get_surrounding_doy(doY) #extracting +-15days for a given doY to prepare a sliding window
    
    rep(1:365,55) %>% is_in(sliding.window) %>% which() %>% obs.ano.grid[.,] -> obs.ano.grid.swData #extracting anomalies for the sliding window for each year (31x55)
    rep(1:365,55) %>% is_in(sliding.window) %>% which() %>% his.ano.grid[.,] -> his.ano.grid.swData #extracting anomalies for the sliding window for each year (31x55)
    
    #computing adjustment factor function, by sampling q with 99 values starting from 0.01 to 0.99 by steps of 0.01
    if(varType == "additive")
      adj.fac.doY[[doY]]<-data.frame(q=seq(0.01,0.99,0.01),
                                     A=lapply(seq(0.01,0.99,0.01), function(q) quantile(obs.ano.grid.swData,q)-quantile(his.ano.grid.swData,q))%>% unlist())
    else
      adj.fac.doY[[doY]]<-data.frame(q=seq(0.01,0.99,0.01),
                                     A=lapply(seq(0.01,0.99,0.01), function(q) quantile(obs.ano.grid.swData,q)/quantile(his.ano.grid.swData,q))%>% unlist())
  }
  return(adj.fac.doY)
}

rank <- comm.rank()
size <- comm.size()
chunks   <- split(1:1440, cut(seq_along(1:1440), size, labels = FALSE))
my_chunk <- chunks[[rank + 1]]  # +1 because R is 1-based
result_list <- lapply(my_chunk, calc_adjFac)

gathered_data <- gather(result_list, rank.dest = 0)
if (rank == 0) 
{
  #reading adjustment factors
  adjFac<-list()
  grid<-1
  for(i in 1:256)
  {
    print(i)
    ls<-gathered_data[[i]]
    for(j in 1:length(ls))
    {
      adjFac[[grid]]<-ls[[j]]
      grid<-grid+1
    }
  }
  
  
  saveRDS(adjFac, file = paste0("/scratch/aniruddha_s_hy.iitr/cmip6/repo/",model.name,"/historical/",variable,"/adjFac.RDS"))
}
finalize()
