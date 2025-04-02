# Ccr5delta32_analysis

The scripts are based on the method described in [1].
The file data/November_S1.csv is obtained from [2] and ETOPO5.rds from [3].

To reproduce the ccr5delta32 analysis, clone this repository and run the script stepadna.R using R:

# Reproducing the results
Start of by cloning the repository:
```
git clone https://github.com/RasaMukti/ccr5delta32_analysis
```

## Environment
To create a correct environment for the analysis create a conda environment with the correct versions of the packages by running the following:
```
conda create -n ccr5delta32_analysis \
  -c conda-forge -c bioconda -c defaults \
  R-base=4.4.2 \
  graphicsmagick=1.3.26 \
  libgdal=3.10.0 \
  udunits2=2.2.28 \
  r-maptools=1.1_8 \
  r-animation=2.7 \
  r-rgeos=0.6_4 \
  r-raster=3.6_30 \
  r-codetools=0.2_20 \
  r-geomap=2.5_11 \
  r-pak==0.8.0.1
```
Then activate the conda environment
```
conda activate ccr5delta32_analysis
```
Then install the R packages not available on conda using pak, this can be done by running the following in R:
```
pkgs <- c(
  "geomapdata@2.0.2", 
  "ReacTran@1.4.3.1",
  "geosphere@1.5.20",
  "deSolve@1.40",
  "magrittr@2.0.3",
  "cowplot@1.1.3",
  "animation@2.7",
  "RColorBrewer@1.1.3",
  "fields@16.3"
)
# Install packages
pak::pkg_install(pkgs)
```
## Analysis
To rerun the analysis run the stepadna.R script.
The `-f` flag take in a CSV file of the following format  
```
"","Sample","genotype","latitude","longitude","age"
"1","BOT14",0,53.17,"67.67",5263.5
"2","BOT15",0,53.17,"67.67",5134.5
"3","BOT2016",0,53.17,"67.67",5450
```
The `-o` flag determines the outfile  
The `-c` flag determines the number of cores to use  
The `-a` flag the estimated age of the allele  
the `-l` flag determines whether to use pseudohaploid genotypes (ph) or genotype likelihoods (gl)"  

The input files (passed to the program by the `-f` flag) used in our analysis can be found in the `/data` directory  
To rerun our analysis run the following:
```
# HAPI output files without filter
Rscript stepadna.R -f data/HAPI_input.csv -a 16128 -o output/out_HAPI.csv -l gl -i 50 -c 50 
# HAPI output file with permissive filter
Rscript stepadna.R -f data/permissive_input.csv -a 8540 -o output/out_permissive.csv -l gl -i 50 -c 50
# HAPI output file with strict filter
Rscript stepadna.R -f data/strict_input.csv -a 6748 -o output/out_strict.csv -l gl -i 50 -c 50
```
## Results
In order to generate the allele frequency trajectory maps of the results generated in the previous step run the scirpt stepplots.R:

```
Rscript stepplots.R
```

## References
[1] Rasa A Muktupavela, Martin Petr, Laure Ségurel, Thorfinn Korneliussen, John Novembre, Fernando Racimo (2022) Modeling the spatiotemporal spread of beneficial alleles using ancient genomes eLife 11:e73767
https://doi.org/10.7554/eLife.73767

[2] Novembre J Galvani AP Slatkin M (2005) The geographic spread of the CCR5 delta32 HIV-resistance allele PLOS Biology 3:11. https://doi.org/10.1371/journal.pbio.0030339

[3] NOAA, National Geophysical Data Center, B. C (1988) Data Announcement 88-MGG-02,, Digital Relief of the Surface of the Earth Publisher: U.S. Department of Commerce.
