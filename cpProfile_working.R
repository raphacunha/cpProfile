# set wd
setwd("Z:/sandbox/cpProfile/")
matlab.script <- dir()[grep("\\.m", dir())]

# load function 
library(matconv)
cp_R <-mat2r(matlab.script)
str(cp_R)

# save to text file
file.create("cp_R2.R")
fileConn <- file("cp_R2.R")
writeLines(cp_R$rCode, fileConn)
close(fileConn)
file.show("cp_R2.R")

# create dummy data to test translated snipets
dummy.mat <- matrix(2,length(LETTERS),length(LETTERS),dimnames=list(LETTERS,LETTERS))
dummy.mat2 <- Matrix::Matrix(dummy.mat)
dummy.g <- igraph::graph.adjacency(dummy.mat)
dummy.net <- intergraph::asNetwork(dummy.g)
net.object <- dummy.net


# set directories
ROOTDIR <- "~/Github/cpProfile"

# load world trade network as example
wtn <- read.csv(file.path(ROOTDIR, "wtn_lscc.csv"),
                header = FALSE)
wtn <- as.matrix(wtn)
labels <- read.table(file.path(ROOTDIR, "labels.txt"), header = TRUE)
row.names(wtn) <- labels$labels
directed <- T

###
### START CODE
###

# load required packages
require(Matrix)
require(network)
require(lattice)
require(igraph)
require(intergraph)
require(pracma)

# matlab translated code
cpProfile <- function(net.object, directed = NULL) {   # leave function closed until translation is completed
  
  # create netname from matrix, igraph, or network objects and assign labels
  if(class(net.object) %in% c("matrix", "dsyMatrix", "dscMatrix", "dsparseMatrix", "dsRMatrix", "dtCMatrix", "dtpMatrix", "dtRMatrix", "dtrMatrix")) { 
    print("Matrix object detected. Converting to sparse Matrix format")
    if(is.null(rownames(net.object))) { labels <- 1:nrow(net.object) }
    else { labels <- rownames(net.object) }
    netname <- Matrix::Matrix(net.object, dimnames = list(labels,labels)) }
  
  if(class(net.object) == "igraph") {
    print("igraph object detected. Converting to sparse Matrix format")
    try(label.value <- grep("name", igraph::list.vertex.attributes(net.object), value = T)[1])
    if(length(label.value) == 0 | is.null(label.value)) { labels <- 1:igraph::vcount(net.object) }
    else { labels <- igraph::get.vertex.attribute(net.object, label.value) }
    netname <- Matrix::Matrix(igraph::get.adjacency(net.object), dimnames = list(labels,labels)) }
  
  if(class(net.object) == "network") {
    print("Network object detected. Converting to sparse Matrix format")
    try(label.value <- grep("name", network::list.vertex.attributes(net.object), value = T)[1])
    if(length(label.value) == 0 | is.null(label.value)) { labels <- 1:network::network.size(net.object) }
    else { labels <- network::get.vertex.attribute(net.object, label.value) }
    netname <- Matrix::Matrix(igraph::get.adjacency(intergraph::asIgraph(net.object)), dimnames = list(labels,labels)) }
  
  # check if convertion was sucesscul
  if(!exists("netname")) { stop(print("cpProfile requires a matrix, network, or igraph object")) }
  
  # start routine
  if(any(netname > 1)) { cat("Weighted network detected\n\n") }
  cat("Profiling Core-Periphery\n\n")
  
  # calculate algorithm metrics
  A <- netname
  k_in <- colSums(A)                # row vector of node in-weights (or in-degrees)
  k_out <- rowSums(A)               # column vector of node out-weights (or out-degrees)
  k_tot <- k_out + k_in             # total degree (or twice the degree, if undirected)
  m <- sum(k_in)                    # total weight (or total number of links) in the network
  N <- length(k_in)                 # number of nodes
  Abin <- as.logical(A)             # binary adjacency matrix
  #if(directed == T) {directed <- 1} # use info for directed network
  #if(directed == F) {directed <- 0} # use info for undirected network
  if(is.null(directed)) {
    if(any(rowSums(A) == colSums(A))) {
      directed <- 0 }           # identify undirected network
    else {
      directed <- 1 }           # identify directed network
  }
  if(any(A > 1)) { weighted <- 1 }  # identify weighted networks
  if(!any(A > 1)) { weighted <- 0}  # identify weighted networks
  
  # compute Markov matrix
  cat("Computing the Markov matrix\n\n")
  
  # create the Markov matrix by row-normalizing A
  P <- Matrix::Matrix(0, length(rownames(A)), length(colnames(A)),
                      dimnames = list(rownames(A), colnames(A)))
  P <- (diag(1./k_out)) %*% A
  
  # compute Markov asymptotic distribution (x)
  cat('Computing Markov chain asymptotic distribution...')
  
  if(directed == T) {
    AAA <- diag(N) - t(P)
    # A[N,,drop = FALSE] = 1
    AAA <- as.matrix(AAA)
    AAA[N,] <- 1
    bbb <- matrix(0, N, 1)
    bbb[N] <- 1
    x <- pracma::mldivide(AAA, bbb) # mldivide(A, B) corresponds to A\B in Matlab notation
    rownames(x) <- labels
  } else { x = k_in/sum(k_in) }
  
  # [xP]_ij = x_i * p_ij
  xP <-  diag(as.vector(x), N, N) %*% P
  dimnames(xP) <- list(labels,labels)
  A <- as.matrix(A)
  
  # sorting nodes according to total degree
  if(class(labels) %in% c("integer", "numeric")){
    L <- sort(k_in + k_out)
    nodelist <- as.numeric(names(L))
  }
  if(class(labels) == "character"){
    L <- k_in + k_out
    names(L) <- 1:N
    L <- sort(L)
    nodelist <- as.numeric(names(L))
  }
  
  # current periphery and core
  # starting periphery from the least connected node
  periph <- nodelist[1]
  core <- setdiff(1:N, periph)
  
  # alpha_tmp(i) is the persistence probability of the periphery
  # after the i-th node has been added
  alpha_tmp <- matrix(0, 1, N)
  
  # at each cycle, introducing the node that yields the
  # smallest increase in the pers.prob. of the periphery
  x_sum <- sum(x[periph])
  xP_sum <- sum(sum(xP[periph,periph]))
  #s_sum <- k_out[periph]
  #w_sum <- 0
  
  for (i in 2:(N-1)){
    
    if (pracma::rem(i,100) == 0){
      cat(paste('...adding node', as.character(i), 'of', as.character(N)))
    }
    
    utest <- matrix(0, 1, length(core))
    
    for (j in 1:length(core)){
      # computing the pers. prob. if node j is adedd to periphery
      utest[j] <- (xP_sum +
                     sum(xP[core[j],periph]) +
                     sum(xP[periph,core[j]]) +
                     xP[core[j],core[j]]) / (x_sum + x[core[j]])
    }
    
    colnames(utest) <- 1:ncol(utest)
    uuu <- utest[,order(utest[1,])]
    jjj <- as.numeric(names(uuu))
    
    alpha_tmp[i] <- uuu[1]
    
    # among the core nodes yielding minimal increase in the pers. prob.,
    # select the one with smallest total degree
    listmin <- core[jjj[uuu == min(uuu)]]
    k_listmin <- t(k_tot[listmin])
    
    colnames(k_listmin) <- 1:ncol(k_listmin)
    kkk <- min(k_listmin)
    lll <- as.numeric(colnames(k_listmin)[k_listmin == min(k_listmin)])
    newnode <- listmin[lll]
    
    periph <- c(periph, newnode)
    core <- setdiff(1:N, periph)
    x_sum <- sum(x[periph])
    xP_sum <- sum(sum(xP[periph,periph]))
    
  }
  
  # final step: the current periphery eventually includes the whole network
  alpha_tmp[N] <- 1
  periph <- c(periph, core)
  
  # CP-centralization C
  C <- sum((0:(N-1))/(N-1)-alpha_tmp)/((N-2)/2)
  
  # for export purposes: alpha[i] is the coreness of node i
  alpha <- data.frame(alpha = rep(NA, N),
                      node_label = rep(NA, N))
  
  for (i in 1:N){
    alpha[i,"alpha"] <- alpha_tmp[i]
    alpha[i,"node_label"] <- as.character(labels[periph[i]])
  }
  
  # return list with alpha_i and C 
  cpProfile <- list(alpha = alpha,
                    C = C)
  return(cpProfile)
  
}

cp_object <- cpProfile(net.object = wtn)

alpha[order(alpha$alpha), ]

coreness <- cp_object[["alpha"]]
coreness <- coreness[order(coreness$alpha, decreasing = TRUE), ]
head(coreness)

tail(cp_object[["alpha"]])


# plot

require(reshape2)

alpha_plot <- cp_object[["alpha"]]
alpha_plot <- alpha_plot[order(alpha_plot$alpha), ]
N <- dim(wtn)[1]
alpha_plot$N <- 1:N
alpha_plot$id_line <- (0:(N-1))/(N-1) #identity line
alpha_plot <- melt(alpha_plot, id.vars =  c("N", "node_label"))

require(ggplot2)


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




p <- ggplot(alpha_plot, aes(x = N, y = value,
                            group = variable,
                            colour = variable,
                            shape = variable)) +
  geom_line() +
  geom_point() +
  scale_colour_manual(values = c("red", "black")) +
  scale_shape_manual(values = c(19, NA)) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  theme_bw(base_family = "Myriad Pro") +
  theme(legend.position = "none",
        panel.grid.minor = element_blank()) +
  labs(y = expression("Core-periphery profile"~alpha[k]),
       x = expression("Number of nodes of"~P[k])) +
  NULL
p

<img src="https://imgur.com/yNLr5Nt.png" width="30%">

plot(k_tot,alpha[,1])


disp(['CPU time (main cycle only) [sec] <- ',num2str(ttt)])
disp([' '])
disp(['Press a key to display node coreness...'])
pause
disp(['rank     node label     id     coreness alpha_k'])
for (i in N:-1:1){
    #     disp([int2str(N-i+1),'   ',char(labels(periph(i))),'   ',char(idc(periph(i))),'   ',num2str(alpha_tmp(i))])
    disp([int2str(N-i+1),'   ',char(labels(periph(i))),'   ',num2str(alpha_tmp(i))])
}
disp(['rank     node label     coreness alpha_k'])

#saving the whole workspace
save(strcat('WksCP_',netname,'.mat'))

figure
plot(k_in,alpha,'o')


}




