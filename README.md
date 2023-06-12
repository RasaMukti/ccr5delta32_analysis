# ccr5delta32_analysis

The scripts are based on the method described in [1].
The file data/November_S1.csv is obtained from [2] and ETOPO5.rds from [3].

To reproduce the ccr5delta32 analysis, clone this repository and run the script stepadna.R using R:

## Analysis
```
Rscript stepadna.R -f data/HAPI_input.csv -a 17402 -o output/out_HAPI.csv -l gl -i 50 -c 50
Rscript stepadna.R -f data/permissive_input.csv -a 9128 -o output/out_permissive.csv -l gl -i 50 -c 50
Rscript stepadna.R -f data/strict_input.csv -a 7714 -o output/out_strict.csv -l gl -i 50 -c 50
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
