#'
#' cpProfile implements the core-periphery profile for networks as developed in:
#' 
#' Della Rossa, F., Dercole, F., & Piccardi, C. 2013. "Profiling core-periphery network
#' structure by random walkers." Scientific Reports, 3, 1467. DOI: 10.1038/srep01467
#' 
#' Original code in Matlab provided by the authors at:
#' http://home.deib.polimi.it/piccardi/PCPNSbyRW.html
#'


# load required packages
require(Matrix)
require(network)
require(igraph)
require(pracma)


cpProfile <- function(net.object, directed = NULL) {
  
  # create netname from matrix, igraph, or network objects and assign labels
  if(class(net.object) %in% c("matrix", "dsyMatrix", "dscMatrix", "dsparseMatrix", "dsRMatrix", "dtCMatrix", "dtpMatrix", "dtRMatrix", "dtrMatrix")) { 
    cat("Matrix object detected. Converting to sparse Matrix format\n\n")
    # check if network is connected
    net.object.g <- igraph::graph_from_adjacency_matrix(net.object, mode = directed, weighted = TRUE)
    if(!igraph::is_connected(net.object.g)) { 
      cat("Input is not a connected network. Decomposing to largest connected component\n\n")
      net.object <- igraph::decompose(net.object.g, max.comps = 1)[[1]]
      net.object <- igraph::get.adjacency(net.object)
    }
    if(is.null(rownames(net.object))) { labels <- 1:nrow(net.object) }
    else { labels <- rownames(net.object) }
    netname <- Matrix::Matrix(net.object, dimnames = list(labels,labels)) }
  
  if(class(net.object) == "igraph") {
    cat("igraph object detected. Converting to sparse Matrix format\n\n")
    # check if network is connected
    if(!igraph::is_connected(net.object)) { 
      cat("Input is not a connected network. Decomposing to largest connected component\n\n")
      net.object <- igraph::decompose(net.object, max.comps = 1)[[1]]
    }
    try(label.value <- grep("name", igraph::list.vertex.attributes(net.object), value = TRUE)[1])
    if(length(label.value) == 0 | is.null(label.value)) { labels <- 1:igraph::vcount(net.object) }
    else { labels <- igraph::get.vertex.attribute(net.object, label.value) }
    netname <- Matrix::Matrix(igraph::get.adjacency(net.object), dimnames = list(labels,labels)) }
  
  if(class(net.object) == "network") {
    cat("Network object detected. Converting to sparse Matrix format\n\n")
    try(label.value <- grep("name", network::list.vertex.attributes(net.object), value = TRUE)[1])
    if(length(label.value) == 0 | is.null(label.value)) { labels <- 1:network::network.size(net.object) }
    else { labels <- network::get.vertex.attribute(net.object, label.value) }
    # check if network is connected
    net.object <- intergraph::asIgraph(net.object)
    if(!igraph::is_connected(net.object)) { 
      cat("Input is not a connected network. Decomposing to largest connected component\n\n")
      net.object <- igraph::decompose(net.object, max.comps = 1)[[1]]
    }
    netname <- Matrix::Matrix(igraph::get.adjacency(net.object), dimnames = list(labels,labels)) }
  
  # check if convertion was sucesscul
  if(!exists("netname")) { stop(print("cpProfile requires a matrix, network, or igraph object")) }
  
  # start routine
  if(any(netname > 1)) { cat("Weighted network detected\n\n") }
  cat("Profiling Core-Periphery\n\n")
  
  # calculate algorithm metrics
  A <- netname
  k_in <- colSums(A)                    # row vector of node in-weights (or in-degrees)
  k_out <- rowSums(A)                   # column vector of node out-weights (or out-degrees)
  k_tot <- k_out + k_in                 # total degree (or twice the degree, if undirected)
  m <- sum(k_in)                        # total weight (or total number of links) in the network
  N <- length(k_in)                     # number of nodes
  Abin <- as.logical(A)                 # binary adjacency matrix
  
  if(is.null(directed)) {
    if(any(rowSums(A) == colSums(A))) {
      directed <- 0 }                   # identify undirected network
    else {
      directed <- 1 }                   # identify directed network
  } else {
    if(directed == TRUE) { directed <- 1 } # use info for directed network
    if(directed == FALSE) { directed <- 0 } # use info for undirected network}
  } 
  if(any(A > 1)) { weighted <- 1 }      # identify weighted networks
  if(!any(A > 1)) { weighted <- 0 }     # identify weighted networks
  
  # compute Markov matrix
  cat("Computing the Markov matrix\n\n")
  
  # create the Markov matrix by row-normalizing A
  P <- Matrix::Matrix(0, length(rownames(A)), length(colnames(A)),
                      dimnames = list(rownames(A), colnames(A)))
  P <- (diag(1./k_out)) %*% A
  
  # compute Markov asymptotic distribution (x)
  cat('Computing Markov chain asymptotic distribution...')
  
  if(directed == TRUE) {
    AAA <- diag(N) - t(P)
    # A[N,,drop = FALSE] = 1
    AAA <- as.matrix(AAA)
    AAA[N,] <- 1
    bbb <- matrix(0, N, 1)
    bbb[N] <- 1
    x <- pracma::mldivide(AAA, bbb) # mldivide(A,B) is A\B in Matlab notation
    rownames(x) <- labels
  } else { x = k_in/sum(k_in) }
  
  # [xP]_ij = x_i * p_ij
  xP <-  diag(as.vector(x), N, N) %*% P
  dimnames(xP) <- list(labels,labels)
  A <- as.matrix(A)
  
  # sorting nodes according to total degree
  if(class(labels) %in% c("integer", "numeric")) {
    L <- sort(k_in + k_out)
    nodelist <- as.numeric(names(L))
  }
  if(class(labels) == "character") {
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
  # smallest increase in the pers. prob. of the periphery
  x_sum <- sum(x[periph])
  xP_sum <- sum(sum(xP[periph,periph]))
  
  for (i in 2:(N-1)) {
    
    if (pracma::rem(i,100) == 0) {
      cat(paste('...adding node', as.character(i), 'of', as.character(N)))
    }
    
    utest <- matrix(0, 1, length(core))
    
    for (j in 1:length(core)) {
      # computing the pers. prob. if node j is adedd to periphery
      utest[j] <- (xP_sum +
                     sum(xP[core[j],periph]) +
                     sum(xP[periph,core[j]]) +
                     xP[core[j],core[j]]) / (x_sum + x[core[j]])
    }
    
    try(colnames(utest) <- 1:ncol(utest), silent = TRUE)
    try(uuu <- utest[,order(utest[1,])], silent = TRUE)
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
  
  for (i in 1:N) {
    alpha[i,"alpha"] <- alpha_tmp[i]
    alpha[i,"node_label"] <- as.character(labels[periph[i]])
  }
  
  # return list with alpha_i and C 
  cp_profile <- list(alpha = alpha, C = C)
  return(cp_profile)
}