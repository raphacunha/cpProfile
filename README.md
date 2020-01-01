# cpProfile

R implementation of core-periphery profile in networks

Calculations detailed at:

Della Rossa, F., Dercole, F., & Piccardi, C. (2013). Profiling core-periphery network structure by random walkers. Scientific Reports, 3, 1467.
DOI: 10.1038/srep01467
URL: http://www.nature.com/articles/srep01467#supplementary-information

Originally implemented in Matlab: http://home.deib.polimi.it/piccardi/PCPNSbyRW.html

## Example

Load the world trade network adjency matrix and node labels from Della Rossa et al.:

```r
DATADIR <- '' # set your data directory
wtn <- read.csv(file.path(DATADIR, "wtn.csv"), header = FALSE)
wtn <- as.matrix(wtn)
labels <- read.table(file.path(DATADIR, "labels.txt"), header = TRUE)
row.names(wtn) <- labels$labels
```r
