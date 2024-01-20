# MELTRONIC 

_a statistical framework to detect chromatin melting and condensation from GAM data_

<img src="./data/MELTRONIC_schematic.png" width="900">

### Required packages
```r
library(data.table)
library(dplyr)
library(argparser)
library(stringr)q
```

### Available command line applications:
- command_line_apps/matrix_wide_to_long.R: 
    Converts a square matrix into a long matrix for IS calculation. Accepts wildcards for processing of multiple chromosomes.
  - command_line_apps/long_matrix_to_IS.R:
    Calculates insulation scores at multiple distances (default 100kb - 1Mb, steps of 100kb).
- MELTRONIC.R:   
   Compares insulation score (IS) distributions over regions of interest. Was applied to 120 kb sliding windows accross the entire genome and and long genes for the preparation of the manuscript.    
    
type
```bash
Rscript command_line_apps/matrix_wide_to_long.R --help 
```
for explanations
   
   
Developed and tested with R version 3.6.0 Planting of a Tree.  
Developed and maintained by Dominik Szabó [<img src="https://cloud.githubusercontent.com/assets/1810515/4228292/6b03dc88-3958-11e4-9094-d3c1771ccfea.png" width="15">](https://orcid.org/0000-0001-8109-5088).  
Please get in touch for questions and issues: dominik.szabo at mdc-berlin.de  
