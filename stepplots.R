setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

list.of.packages <- c("GEOmap","geomapdata","ReacTran","geosphere","deSolve","magrittr","maptools","ggplot2","cowplot","animation","optparse")
new.packages <-
  list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])]
if (length(new.packages)) {
  install.packages(new.packages,repos = "http://cran.us.r-project.org")
}


# Import libraries
library(GEOmap)

library(geomapdata)
library(ReacTran)
library(geosphere)
library(deSolve)

library(magrittr)
library(maptools)
library(ggplot2)
library(cowplot)
library(animation)
source("./plotfunctions.R")

######################################## Load data ######################################## 

for (dn in c(1,2,3)){

datatype <- c("HAPI", "permissive", "strict")

ccr5_data <- read.csv(paste0("data/",datatype[dn],"_input.csv"))

ccr5_data <- ccr5_data[(ccr5_data["longitude"] <= 120 &
                      ccr5_data["longitude"] >= -10 &
                      ccr5_data["latitude"] <= 75 &
                      ccr5_data["latitude"] >= 30),]

ccr5_data$longitude <- ccr5_data$longitude+rnorm(nrow(ccr5_data),0,0.3)
ccr5_data$latitude <- ccr5_data$latitude+rnorm(nrow(ccr5_data),0,0.3)

ccr5_data <- ccr5_data[(ccr5_data["genotype"] == 1 | ccr5_data["genotype"] == 2),]

######################################## Producing PDFs ######################################## 
opt=list(file=paste0("output/out_",datatype[dn],".csv"), nimage=7, out=paste0("figure_",datatype[dn]))

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

# Initialize some variables
readyFile <- read.csv(opt$file)

Nx <- nrow(topo)
Ny <- ncol(topo)
x0 <- Nx+1-which(round(topoLat) == readyFile[["latitude"]][1])[1]
y0 <- which(round(topoLon) == readyFile[["longitude"]][1])[1]
D = 2.5
timeP = round(readyFile[["alleleAge.years"]][1]/290) 
timeSplit = round(2000/290)

# Obtain the spatial distribution for time periods before 2000 years
outOld <- FischerSolver(Nx, Ny, readyFile[["s"]][1], readyFile[["sigmax"]][1], readyFile[["sigmay"]][1], D, x0, y0, readyFile[["xtran"]][1], readyFile[["ytran"]][1], (timeP-timeSplit), oldnew="old")
outOld <- outOld[,2:ncol(outOld)]
initMat <- matrix(unlist(outOld[nrow(outOld),]),nrow=Nx,ncol=Ny)
initMat[which(topo < 0, arr.ind = TRUE)] <- 0

# Obtain the spatial distribution for time periods after 2000 years
outNew <- FischerSolver(Nx, Ny, readyFile[["s"]][4], readyFile[["sigmax"]][4], readyFile[["sigmay"]][4], D, x0, y0, readyFile[["xtran"]][4], readyFile[["ytran"]][4], timeSplit, oldnew="new", initMat=initMat)
outNew <- outNew[,2:ncol(outNew)]

# Create the GIF and save it
savePDF(col_num = 15, image_num = opt$nimage, outName = opt$out, diffusion_mat = rbind(outOld, outNew), timeP = timeP, Nx = Nx, Ny = Ny, lo1, lo2, la1, la2,
        readyFile[["longitude"]][1], readyFile[["latitude"]][1], ccr5_data$longitude, ccr5_data$latitude, ccr5_data$age, datatype[dn])
 
}
