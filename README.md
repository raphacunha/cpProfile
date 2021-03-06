# cpProfile

R implementation of core-periphery profile in networks

Calculations detailed at:

Della Rossa, F., Dercole, F., & Piccardi, C. (2013). Profiling core-periphery network structure by random walkers. Scientific Reports, 3, 1467. DOI: 10.1038/srep01467. URL: http://www.nature.com/articles/srep01467#supplementary-information

Originally implemented in Matlab: http://home.deib.polimi.it/piccardi/PCPNSbyRW.html

## Example

Load the world trade network adjency matrix and node labels from Della Rossa et al.:

```r
DATADIR <- '' # set your data directory
wtn <- read.csv(file.path(DATADIR, "wtn_lscc.csv"), header = FALSE)
wtn <- as.matrix(wtn)
```

Load node labels (optional):

```r
wtn_labels <- read.table(file.path(DATADIR, "wtn_labels.txt"), header = TRUE)
row.names(wtn) <- wtn_labels$labels
```

Load required packages:

```r
library(Matrix)
library(network)
library(igraph)
library(pracma)
```

Calculate the core-periphery profile and the core-periphery centralization C of the world trade network:

```r
cp_wtn <- cpProfile(net.object = wtn, directed = TRUE)
```

Obtain the cp-centralizaion C of the world trade network:

```r
> cp_wtn[["C"]]
[1] 0.819075
```

Top nodes based on coreness (alpha_k):

```r
> coreness <- cp_wtn[["alpha"]]
> coreness <- coreness[order(coreness$alpha, decreasing = TRUE), ]
> head(coreness)

        alpha            node_label
181 1.0000000         United States
180 0.8551008               Germany
179 0.7548008 China, P.R.: Mainland
178 0.6940441                France
177 0.6392640        United Kingdom
176 0.5991778                 Japan
```

Finally, plot the core-periphery profile:

```r
library(reshape2)
library(ggplot2)

plot_dat <- cp_wtn[["alpha"]]
plot_dat <- plot_dat[order(plot_dat$alpha),]
N <- dim(wtn)[1] # number of nodes
plot_dat$N <- 1:N
plot_dat$id_line <- (0:(N-1))/(N-1) # identity line
plot_dat <- melt(plot_dat, id.vars = c("N", "node_label"))

p <- ggplot(alpha_plot,
            aes(x = N,
                y = value,
                group = variable,
                colour = variable)) +
  geom_line() +
  scale_colour_manual(values = c("red", "black")) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  theme(legend.position = "none") +
  labs(y = expression("Core-periphery profile"~alpha[k]),
       x = expression("Number of nodes of"~P[k]))
p
```

Which should produce a cp profile plot like this:

<img src="cp_wtn.png" width="40%">

