library(optparse)
option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL,
              help="dataset file name"),
  make_option(c("-o", "--out"), type="character", default="out.csv",
              help="output file name [default= %default]"),
  make_option(c("-a", "--age"), type="numeric", default=NULL,
              help="age of the allele"),
  make_option(c("-c", "--cores"), type="integer", default=1,
              help="number of cores to us in parallel"),
  make_option(c("-i", "--init"), type="integer", default=2,
              help="number of initial points for simulated annealing"),
  make_option(c("-l", "--loglik"), type="character", default="ph",
              help="whether to use pseudohaploid genotypes (ph) or genotype likelihoods (gl)"),
  make_option(c("--lon"), type="numeric", default=NULL,
              help="longitude of the origin of the allele. Has to be in the range from -10 to 80"),
  make_option(c("--lat"), type="numeric", default=NULL,
              help="latitude of the origin of the allele. Has to be in the range from 30 to 75"),
  make_option(c("--seed"), type="numeric", default=5,
              help="Seed for optimization")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

library(deSolve)
library(ReacTran)
library(RColorBrewer)
library(geosphere)
library(fields)
library(stats4)
library(GEOmap)
library(geomapdata)
library(plsgenomics)
library(paramtest)
library(bbmle)
library(parallel)
library(maps)
library(mapdata)
library(mapplots)
library(rworldmap)
library(scales)
library(FME)
library(dplyr)
source("./stepfunctions.R")
options(scipen = 999)


# Define map boundaries
la1 <- 30
la2 <- 75
lo1 <- -10
lo2 <- 120

# Create topographical data grid
topo  <- createTopo(la1, la2, lo1, lo2)
topoLon <- topo$topoLon
topoLat <- topo$topoLat
topo <- topo$topo

# Load data file
ts_data <- read.csv(opt$file)

# Convert to generations
ts_data$age <- ts_data$age/29

# Round times and scale by 10
ts_data$age <- round(ts_data$age/10)

present_ccr5 <- read.csv("data/November_S1.csv", sep=";")

present_ccr5 <- present_ccr5[,c("Longitude", "Latitude", "X.32.Frequency", "Sample.Size")]

present_ccr5["Sample.Size"] <- present_ccr5["Sample.Size"]*2
colnames(present_ccr5) <- c("longitude", "latitude", "frequency", "chrom_count")
present_ccr5["frequency"] <- as.numeric(gsub(",", ".", present_ccr5$frequency))

present_ccr5 <- present_ccr5[(present_ccr5["longitude"] <= 80 &
                                present_ccr5["longitude"] >= -10 &
                                present_ccr5["latitude"] <= 75 &
                                present_ccr5["latitude"] >= 30),]

ts_data <- ts_data[(ts_data["longitude"] <= 80 &
                      ts_data["longitude"] >= -10 &
                      ts_data["latitude"] <= 75 &
                      ts_data["latitude"] >= 30),]

#ts_data[, "genotype"] <- apply(ts_data[,c("AA", "AD", "DD")], 1, which.max)-1

# Points with potential initial allele
lonseq <- seq(-5, 100, 10)
latseq <- seq(37, 70, 5)
Lonpoints <- rep(lonseq, length(latseq))
Latpoints <- rep(latseq, each=length(lonseq))
ptsremoved <- c(2,3,4,6,8,10,11,13,14,16,18,20,21,22,23,26,28,30,32,33,34,36,38,40,42,43,44,45,46,48,50,52,54,55,56,57,58,60,62,64,65,66,67,68,70,72,74,76,77)
Lonpoints <- Lonpoints[-ptsremoved]
Latpoints <- Latpoints[-ptsremoved]
Lonpoints[25] <- Lonpoints[25]+1
Lonpoints[14] <- Lonpoints[14]+1

# Initialise parameters
Nx <- nrow(topo)
Ny <- ncol(topo)
D = 2.5 # Initial population density
timeP = round(opt$age/290) 
timeSplit = round(2000/290)

x0_fixed <- Nx+1-which(topoLat == closest(topoLat, opt$lat))
y0_fixed <- which(topoLon == closest(topoLon, opt$lon))

if (is.null(opt$lat) == FALSE & is.null(opt$lon) == FALSE){
  print(c(opt$lat, opt$lon))
} else{print("is null")}

######################################## MLE PH ######################################## 

# Function for calculating likelihoods
LL_PH <- function(iter, Nx=Nx, Ny=Ny, s=s, sigmax=sigmax, sigmay=sigmay, D=D, x0=x0, y0=y0, xtran=xtran, ytran=ytran, timeP=timeP, timeSplit=timeSplit, oldnew=NULL){
  
  # Set boundaries for SANN algorithm
  s.low <- 0
  s.high <- 0.1
  sigma.low <- 1
  sigma.high <- 100
  tran.low <- -2.5
  tran.high <- 2.5
  if (s < s.low | s > s.high | sigmax < sigma.low | sigmay < sigma.low |
      sigmax > sigma.high | sigmay > sigma.high |
      xtran < tran.low | ytran < tran.low | xtran > tran.high | ytran > tran.high
      | x0 < 1 | y0 < 1 | x0 > Nx | y0 > Ny | (topo[x0,y0] < 0) == TRUE){
    return(c(99999))
  } else{

    x0 <- round(x0)
    y0 <- round(y0)
    
    if (oldnew == 0){
      out <- FischerSolver(Nx, Ny, s, sigmax, sigmay, D, x0, y0, xtran, ytran, (timeP-timeSplit), oldnew = "old")
      out <- out[,2:ncol(out)]
      
      time_periods <- seq(timeP,timeSplit,-1)
      timevect <- timeP+1 - time_periods
    }

    
    if (oldnew == 1){
      out <- FischerSolver(Nx, Ny, s, sigmax, sigmay, D, x0, y0, xtran, ytran, timeSplit, oldnew="new", initMat=initMat)
      out <- out[,2:ncol(out)]
      
      time_periods <- seq(timeSplit,0,-1)
      timevect <- timeSplit+1 - time_periods
    }

    
    l <- 0

    # Optimizer crashes with high s values
    if (nrow(out) != length(time_periods)){
      print("ERRORRRRR")
      return(c(100000))
    }
    else{

      for (tp in 1:length(timevect)){

        mat_out <- matrix(unlist(out[timevect[tp],]),nrow=Nx,ncol=Ny)
        p <- t(mat_out[nrow(mat_out):1,])

        if (time_periods[tp] != 0){
          current_time <- ts_data[which(ts_data$age == time_periods[tp]),]
  
          # Pseudo haploid calculation
          if (nrow(current_time) != 0){
            
            for (i in 1:nrow(current_time)){

              # Extract spatial points matching observed point
              y_topo <- which(topoLon == closest(topoLon, current_time$longitude[i]))[1]
              x_topo <- which(topoLat == closest(topoLat, current_time$latitude[i]))[1]
  
              k <- current_time$genotype[i]/2
              small_num <- 1e-5
              if (p[y_topo,x_topo] == 0 | round(p[y_topo,x_topo], 5) == 0){l <- l + log(small_num^k*(1-small_num)^(1-k))}
              else if (p[y_topo,x_topo] == 1 | round(p[y_topo,x_topo], 5) == 1){l <- l + log((1-small_num)^k*(small_num)^(1-k))}
              else{l <- l + log(p[y_topo,x_topo]^k*(1-p[y_topo,x_topo])^(1-k))}
            }
          }
        }
        # Likelihood calculation for present day
        if (time_periods[tp] == 0){
          current_time <- ts_data[which(ts_data$age == 0),]
  
          for (i in 1:nrow(current_time)){
            y_topo <- which(topoLon == closest(topoLon, current_time$longitude[i]))[1]
            x_topo <- which(topoLat == closest(topoLat, current_time$latitude[i]))[1]
            k <- current_time$genotype[i]
            small_num <- 1e-5
            if (p[y_topo,x_topo] == 0 | round(p[y_topo,x_topo], 5) == 0){l <- l+log(choose(2, k)*small_num^k*(1-small_num)^(2-k))}
            else if (p[y_topo,x_topo] == 1 | round(p[y_topo,x_topo], 5) == 1){l <- l+log(choose(2, k)*(1-small_num)^k*small_num^(2-k))}
            else {l <- l+log(choose(2, k)*p[y_topo,x_topo]^k*(1-p[y_topo,x_topo])^(2-k))}

          }
        }
      }

      print(round(c(s, sigmax, sigmay, xtran, ytran, -l), 3))

      if (is.nan(-l)){
        return(c(99999))
      } else{
        return(c(-l))}
    }
  }
}

######################################## MLE GL######################################## 

# Function for calculating likelihoods
LL_GL <- function(iter, Nx=Nx, Ny=Ny, s=s, sigmax=sigmax, sigmay=sigmay, D=D, x0=x0, y0=y0, xtran=xtran, ytran=ytran, timeP=timeP, timeSplit=timeSplit, oldnew=NULL){
  
  # Set boundaries for SANN algorithm
  s.low <- 0
  s.high <- 0.1
  sigma.low <- 1
  sigma.high <- 100
  tran.low <- -2.5
  tran.high <- 2.5
  x0 <- floor(x0)
  y0 <- floor(y0)
  if (s < s.low | s > s.high | sigmax < sigma.low | sigmay < sigma.low |
      sigmax > sigma.high | sigmay > sigma.high |
      xtran < tran.low | ytran < tran.low | xtran > tran.high | ytran > tran.high
      | x0 < 1 | y0 < 1 | x0 > Nx | y0 > Ny){
    return(c(99999))
    
  } else if ((topo[x0,y0] < 0) == TRUE){
    return(c(99999))
  } else{
    
    if (oldnew == 0){
      out <- FischerSolver(Nx, Ny, s, sigmax, sigmay, D, x0, y0, xtran, ytran, (timeP-timeSplit), oldnew = "old")
      out <- out[,2:ncol(out)]
      
      time_periods <- seq(timeP,timeSplit,-1)
      timevect <- timeP+1 - time_periods
    }
    
    if (oldnew == 1){
      out <- FischerSolver(Nx, Ny, s, sigmax, sigmay, D, x0, y0, xtran, ytran, timeSplit, oldnew="new", initMat=initMat)
      out <- out[,2:ncol(out)]
      
      time_periods <- seq(timeSplit,0,-1)
      timevect <- timeSplit+1 - time_periods
    }
    
    l <- 0
    # Optimizer crashes with high s values
    if (nrow(out) != length(time_periods)){
      print("ERRORRRRR")
      return(c(100000))
    }
    else{
      
      for (tp in 1:length(timevect)){
        
        mat_out <- matrix(unlist(out[timevect[tp],]),nrow=Nx,ncol=Ny)
        p <- t(mat_out[nrow(mat_out):1,])
        
        if (time_periods[tp] != 0){
          current_time <- ts_data[which(ts_data$age == time_periods[tp]),]
          
          # Genotype likelihood calculations
          if (nrow(current_time) != 0){
            for (i in 1:nrow(current_time)){
                
              # Extract spatial points matching observed point
              y_topo <- which(topoLon == closest(topoLon, current_time$longitude[i]))[1]
              x_topo <- which(topoLat == closest(topoLat, current_time$latitude[i]))[1]
              
              n <- 2
              k <- current_time$genotype[i]
              l <- l+log(choose(n, k)*p[y_topo,x_topo]^k*(1-p[y_topo,x_topo])^(n-k))
              
              }
            }
        }
        # Likelihood calculation for present day
        if (time_periods[tp] == 0){
          for (i in 1:nrow(present_ccr5)){
            y_topo <- which(topoLon == closest(topoLon, present_ccr5$longitude[i]))[1]
            x_topo <- which(topoLat == closest(topoLat, present_ccr5$latitude[i]))[1]
            n <- present_ccr5$chrom_count[i]
            k <- round(present_ccr5$chrom_count[i]*present_ccr5$frequency[i])
            small_num <- 1e-5
            if (p[y_topo,x_topo] == 0 | round(p[y_topo,x_topo], 5) == 0){l <- l+bignchoosek(n,k)+k*log(small_num)+(n-k)*log(1-small_num)}
            else if (p[y_topo,x_topo] == 1 | round(p[y_topo,x_topo], 5) == 1){l <- l+bignchoosek(n,k)+k*log(1-small_num)+(n-k)*log(small_num)}
            else {l <- l+bignchoosek(n,k)+k*log(p[y_topo,x_topo])+(n-k)*log(1-p[y_topo,x_topo])}
          }
        }
      }
      
      print(round(c(s, sigmax, sigmay, xtran, ytran, -l), 3))
      
      if (is.nan(-l)){
        return(c(99999))
      } else{
        return(c(-l))}
    }
  }
}

#########################################################################################
if(opt$loglik == "ph"){
  func <- LL_PH
} else if (opt$loglik == "gl"){
  func <- LL_GL
}

start_time <- Sys.time()
if (is.null(opt$lat) == FALSE & is.null(opt$lon) == FALSE){
  ResOld <- param_optim_fixedInit(init_points = opt$init, seed = opt$seed, num_cores = opt$cores, LL=func, oldnew="old")
} else{ResOld <- param_optim_old(init_points = opt$init, seed = opt$seed, num_cores = opt$cores, LL=func)}

CIres <- getCI(ResOld, "old", LL=func)
bestResOld <- CIres[1][[1]]
oldCIvect <- CIres[2:length(CIres)]

outOld <- FischerSolver(Nx, Ny, coef(bestResOld)[["s"]], coef(bestResOld)[["sigmax"]], coef(bestResOld)[["sigmay"]], D, coef(bestResOld)[["x0"]], coef(bestResOld)[["y0"]], coef(bestResOld)[["xtran"]], coef(bestResOld)[["ytran"]], (timeP-timeSplit), oldnew="old")
outOld <- outOld[,2:ncol(outOld)]
initMat <- matrix(unlist(outOld[nrow(outOld),]),nrow=Nx,ncol=Ny)
initMat[which(topo < 0, arr.ind = TRUE)] <- 0

if (is.null(opt$lat) == FALSE & is.null(opt$lon) == FALSE){
  ResNew <- param_optim_fixedInit(init_points = opt$init, seed = opt$seed, num_cores = opt$cores, LL=func, oldnew="new")
} else{ResNew <- param_optim_new(init_points = opt$init, seed = opt$seed, num_cores = opt$cores, LL=func)}

CIresNew <- getCI(ResNew, "new", LL=func)
bestResNEW <- CIresNew[1][[1]]
newCIvect <- CIresNew[2:length(CIresNew)]
end_time <- Sys.time()
print(c("Run time: ", (end_time-start_time)))

outNew <- FischerSolver(Nx, Ny, coef(bestResNEW)[["s"]], coef(bestResNEW)[["sigmax"]], coef(bestResNEW)[["sigmay"]], D, coef(bestResNEW)[["x0"]], coef(bestResNEW)[["y0"]], coef(bestResNEW)[["xtran"]], coef(bestResNEW)[["ytran"]], timeSplit, oldnew="new", initMat=initMat)
outNew <- outNew[,2:ncol(outNew)]

outDF <- data.frame(. = c(paste("before 2000"), "CI low", "CI high"), s = c(round(coef(bestResOld)[["s"]], 4), oldCIvect[[1]], oldCIvect[[2]]),
                    sigmax=c(round(coef(bestResOld)[["sigmax"]], 3), max(oldCIvect[[3]], 0), oldCIvect[[4]]),
                    sigmay=c(round(coef(bestResOld)[["sigmay"]], 3), max(oldCIvect[[5]], 0), oldCIvect[[6]]),
                    xtran=c(round(coef(bestResOld)[["xtran"]], 3), oldCIvect[[7]], oldCIvect[[8]]),
                    ytran=c(round(coef(bestResOld)[["ytran"]], 3), oldCIvect[[9]], oldCIvect[[10]]),
                    longitude=round(topoLon[coef(bestResOld)[["y0"]]]), latitude=round(topoLat[Nx+1-coef(bestResOld)[["x0"]]]),
                    alleleAge.years=timeP*290)
 
outDF <- rbind(outDF, data.frame(. = c(paste("after 2000"), "CI low", "CI high"), s = c(round(coef(bestResNEW)[["s"]], 4), newCIvect[[1]], newCIvect[[2]]),
                    sigmax=c(round(coef(bestResNEW)[["sigmax"]], 3), max(newCIvect[[3]], 0), newCIvect[[4]]),
                    sigmay=c(round(coef(bestResNEW)[["sigmay"]], 3), max(newCIvect[[5]], 0), newCIvect[[6]]),
                    xtran=c(round(coef(bestResNEW)[["xtran"]], 3), newCIvect[[7]], newCIvect[[8]]),
                    ytran=c(round(coef(bestResNEW)[["ytran"]], 3), newCIvect[[9]], newCIvect[[10]]),
                    longitude=round(topoLon[coef(bestResNEW)[["y0"]]]), latitude=round(topoLat[Nx+1-coef(bestResNEW)[["x0"]]]), alleleAge.years=timeP*290))


write.csv(outDF, opt$out, quote=FALSE)
save.image(paste0(opt$out,".RData"))

