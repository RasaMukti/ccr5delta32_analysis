# ccr5delta32_analysis

The scripts are based on the method described in [1].
The file data/November_S1.csv is obtained from [2] and ETOPO5.rds from [3].

To reproduce the ccr5delta32 analysis, clone this repository and run the script stepadna.R using R:

# Environment
Start of by creating a conda environment with the correct versions of the packages:
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
And lastly install the R packages not available on conda, inside the environemnt.
These are: 
```

```




## Analysis
```
Rscript stepadna.R -f data/HAPI_input.csv -a 16128 -o output/out_HAPI.csv -l gl -i 50 -c 50
Rscript stepadna.R -f data/permissive_input.csv -a 8540 -o output/out_permissive.csv -l gl -i 50 -c 50
Rscript stepadna.R -f data/strict_input.csv -a 6748 -o output/out_strict.csv -l gl -i 50 -c 50
```

The first, second and third line runs the analysis for results using HAPI classification, permissive and strict genotype calls, respectively. The code runs using 50 cores (specified by -c option).

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
