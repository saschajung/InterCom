set.seed(100)

############################
#  Intercom Functions
############################


#### prepare.data function prepares the data fro downstream analysis by using expression matrix and annotation table. Same as in commap.

#' @description The function prepares the input data for downstream analysis.
#' @param expr A dataframe with TPM normalized gene expression data, gene names as rownames in form Ensemblid_Genename and
#' cell ids as column names.
#' @param ph A dataframe with two columns. First column contains cell ids corresponding to column names in expr file and
#' second column containing cell type information.
#' @param ph.cell.id.col Name of column in ph with cell ids
#' @param ph.cell.type.col Nmae of column in ph with cell types.
#' @return A dataframe with gene names as rownames and cell types as column names.
#' 

prepare.data <- function(expr,ph,placeholder=NULL,ph.cell.id.col="cell.id",ph.cell.type.col="cell.type"){
  ## Prepares the data of a certain format (see below) to be supplied to the commap.
  ## Format:
  ## rownames(expr) have a format ENSMBLID_GENESYMBOL
  ## colnames(expr) are cell id's
  ## phenotype (ph) data have cell id column (ph.cell.id.col)
  ##  and cell type column (ph.cell.type.col)
  ## Remove ENSMBLID part
  ## print(grep("HLA",rownames(expr),value=T))    
  new.rownames <- sapply(strsplit(rownames(expr),"_"),function(x) x[2])
  
  if(!is.null(placeholder)){
    ## Just remove rows with un-identified gene symbols marked with the placeholder
    idx.to.remove <- new.rownames == placeholder
    expr <- expr[!idx.to.remove,]
    rownames(expr) <- new.rownames[!idx.to.remove]
    ## it will have a warning here automatically about the duplicated rownames
  }
  else {
    ## Just remove duplicates in the rownames, leaving potentially a single unknown placeholder
    idx.to.remove <- duplicated(new.rownames)
    cat(" ",length(which(idx.to.remove)),"duplicated entries found in input data and were removed\n")
    expr <- expr[!idx.to.remove,]
    rownames(expr) <- new.rownames[!idx.to.remove]
  }
  ## Substitue cell.id with cell.type information
  if(!(ph.cell.id.col %in% colnames(ph)) || !(ph.cell.type.col %in% colnames(ph)))
    stop("Indicate the proper names of cell.id and cell.type columns")
  new.colnames <- sapply(colnames(expr),function(x) {
    make.names(as.character(ph[[ph.cell.type.col]][x == as.character(ph[[ph.cell.id.col]])]),
               unique=FALSE)
  })
  colnames(expr) <- new.colnames
  return(expr)
}


get.gene.expr <- function(exp.tbl,genes,cell.type=NULL){
  gene.exp.tbl <- exp.tbl[genes,,drop=FALSE]
  
  all.pops <- cell.type
  
  if(length(all.pops) == 1){
    cell.gene.exp <- gene.exp.tbl[,which(grepl(x = colnames(gene.exp.tbl),pattern = paste0("^",cell.type,"[\\.0-9]*$"),ignore.case = F)),drop=FALSE]
    cell.gene.abs.exp <- rowSums(cell.gene.exp)
    cell.gene.abs.exp <- cbind.data.frame(row.names(cell.gene.exp),cell.gene.abs.exp,stringsAsFactors = F)
    colnames(cell.gene.abs.exp) <- c("gene","abs.expr")
    
    cell.gene.bool <- as.data.frame(bool.data(exp.tbl = cell.gene.exp))
    
    cell.gene.cons <- rowSums(cell.gene.bool)
    cell.gene.cons <- cbind.data.frame(row.names(cell.gene.bool),cell.gene.cons,stringsAsFactors = F)
    colnames(cell.gene.cons) <- c("gene","cell.count")
    cell.gene.cons$exp.perc <- cell.gene.cons$cell.count/dim(cell.gene.bool)[2]
    cell.gene.out <- dplyr::inner_join(x = cell.gene.cons,y = cell.gene.abs.exp, by = "gene")
    cell.gene.out$celltype <- cell.type
    out <- cell.gene.out[,c(5,1:4)]
  }else{
    
    out <- do.call("rbind",lapply(X = all.pops,FUN = function(celltype1){
      # print(celltype1)
      cell.gene.exp <- gene.exp.tbl[,which(grepl(x = colnames(gene.exp.tbl),pattern = paste0("^",celltype1,"[\\.0-9]*$"),ignore.case = F))]
      
      cell.gene.abs.exp <- rowSums(cell.gene.exp) 
      cell.gene.abs.exp <- cbind.data.frame(row.names(cell.gene.exp),cell.gene.abs.exp,stringsAsFactors = F)
      colnames(cell.gene.abs.exp) <- c("gene","abs.expr")
      
      cell.gene.bool <- bool.data(exp.tbl = cell.gene.exp)
      
      cell.gene.cons <- rowSums(cell.gene.bool)
      cell.gene.cons <- cbind.data.frame(row.names(cell.gene.bool),cell.gene.cons,stringsAsFactors = F)
      colnames(cell.gene.cons) <- c("gene","cell.count")
      cell.gene.cons$exp.perc <- cell.gene.cons$cell.count/dim(cell.gene.bool)[2]
      
      cell.gene.out <- dplyr::inner_join(x = cell.gene.cons,y = cell.gene.abs.exp, by = "gene")
      cell.gene.out$celltype <- celltype1
      return(cell.gene.out[,c(5,1:4)])
      
    }))
  }
  invisible(gc())
  return(out)
}


#### bool.data function booleanizes the expression matrix based on expression threshold. It also removes genes not expression in any cell.
## Input : 1. exp.tbl : Gene expression count data matrix
##         2. expr.thrs : Expression threshold for booleanization (default 0).
## Output : Booleanized matrix

#' @description The function booleanizes the expression matrix based on expression threshold. 
#' It also removes genes not expression in any cell.
#' @param exp.tbl Gene expression dataframe.
#' @param expr.thrs Expression threshold.
#' @return Gene expression dataframe with booleanized expression.


bool.data <- function(exp.tbl,expr.thrs=0){
  bool.tbl <- matrix(data = as.numeric(exp.tbl > expr.thrs),nrow = nrow(exp.tbl),ncol = ncol(exp.tbl))
  colnames(bool.tbl) <- colnames(exp.tbl)
  row.names(bool.tbl) <- row.names(exp.tbl)
  bool.tbl <- bool.tbl[which(rowSums(bool.tbl) > 0),,drop=FALSE]
  return(bool.tbl)
}


#### maxsub2d function finds the rectangle in the matrix with maximum sum. The function and its fortran implementation has been from the adagio R package (after correction).
## Input : 1. A : Booleanized matrix (the matrix must have booleanization in 1/-1 form and not 1/0 form.)
## Output : Sum, index and the maximum sum sub-matrix.

#' @description The function finds the maximum sum sub-matrix in a matrix.
#' The function and its fortran implementation has been adapted from the adagio R package (after correction).
#' @param A Matrix with booleanization in 1/-1 form and not 1/0 form.
#' @return List containing Sum, index of the sub-matrix and the sub-matrix.
#'


maxsub2d <- function(A) {
  stopifnot(is.numeric(A), is.matrix(A))
  n <- nrow(A); m <- ncol(A)
  
  if (all(A <= 0))
    stop("At least on element of matrix 'A' must be positive.")
  if (all(A >= 0))
    return(list(sum = sum(A), inds = c(1, n, 1, m), submat = A))
  
  mi <- vector("integer", 4)
  S <- matrix(0, nrow = n+1, ncol = m)
  aa <- numeric(m)
  b <- 0.0
  
  fm <- 0.0
  R <- .Fortran("maxsub2f", as.numeric(A), as.numeric(S),
                as.integer(n), as.integer(m),
                fmax = as.numeric(fm), mind = as.integer(mi),
                as.numeric(aa), as.numeric(b))
  
  fm <- R$fmax
  mi <- R$mind
  
  invisible(gc())
  
  return(list(sum = fm, inds = mi,
              submat = A[mi[1]:mi[2], mi[3]:mi[4], drop = FALSE]))
}

get.cons.tfs <- function(exp.tbl,quantile = 0.95){
  tf.tbl <- exp.tbl[which(row.names(exp.tbl) %in% tfs),]
  tf.bool <- bool.data(exp.tbl = tf.tbl)
  freq_df <- apply(tf.bool,1,sum)
  freq_df <- data.frame(Gene = names(freq_df), Freq = freq_df,stringsAsFactors = FALSE)
  cutoff <- unname(quantile(freq_df$Freq[freq_df$Freq > 0],0.95))
  return(list(tf.max.mat.cell = colnames(exp.tbl), tf.count = freq_df[freq_df$Freq >= cutoff,]))
}

#### get.max.cluster function is used to find the sub-matrix with conserved TFs and receptors in the gene expression matrix.
## Input : 1. exp.tbl : Gene expression matrix

#' @description The function computes the sub-matrix with conserved TFs and receptors in the gene expression matrix
#' @param exp.tbl Gene expression matrix.
#' @return list containing TFs and receptors in the sub-matrix and their cell counts.

get.max.cluster <- function(exp.tbl){
  
  tf.tbl <- exp.tbl[which(row.names(exp.tbl) %in% tfs),]
  
  tf.bool <- bool.data(exp.tbl = tf.tbl)
  
  tf.bool[tf.bool == 0] <- -1
  
  tf.clust <- cluster_matrix(x = tf.bool,dim = "both")
  
  maxsum.tf <- maxsub2d(A = tf.clust)
  
  maxsum.tf.idx <- maxsum.tf$inds
  
  maxsum.matrix.tf <- tf.clust[maxsum.tf.idx[1]:maxsum.tf.idx[2],maxsum.tf.idx[3]:maxsum.tf.idx[4]]
  
  maxsum.matrix.tf[maxsum.matrix.tf == -1] <- 0
  
  if(abs(maxsum.tf.idx[3] - maxsum.tf.idx[4]) < 0.05 * as.numeric(dim(exp.tbl)[2])){
    return("No conserved TFs found.")
  }
  
  maxsum.tf.count <- as.data.frame(rowSums(maxsum.matrix.tf),stringsAsFactors = F)
  
  maxsum.tf.count <- cbind.data.frame(row.names(maxsum.tf.count),maxsum.tf.count,stringsAsFactors = F) 
  
  colnames(maxsum.tf.count) <- c("Gene","Cell.count")
  
  maxsum.tf.count$Cell.perc <- maxsum.tf.count$Cell.count/dim(exp.tbl)[2]
  maxsum.tf.count$TF.matrix.frac <- maxsum.tf.count$Cell.count/dim(maxsum.matrix.tf)[2]
  
  invisible(gc())
  
  return(list("tf.max.mat.cell" = colnames(maxsum.matrix.tf),"tf.count" = maxsum.tf.count))
}


gene.coexp <- function(exp.tbl,gene.set.1,gene.set.2,ncores=4){
  out.frame <- do.call("rbind",lapply(X = gene.set.1, FUN = function(x){
    out.row <- do.call("rbind",parallel::mclapply(mc.cores = ncores,X = gene.set.2, FUN = function(y){
      exp.bool <- bool.data(exp.tbl = exp.tbl)
      Combo <- paste0(x,"_",y)
      Set1.exp <- names(which(exp.bool[x,] == 1 ))
      Set2.exp <- names(which(exp.bool[y,] == 1))
      Set1.not.exp <- names(which(exp.bool[x,] == 0))
      Set2.not.exp <- names(which(exp.bool[y,] == 0))
      n11 <- length(intersect(Set1.exp,Set2.exp))
      n01 <- length(intersect(Set1.not.exp,Set2.exp))
      n10 <- length(intersect(Set1.exp,Set2.not.exp))
      n00 <- length(intersect(Set1.not.exp,Set2.not.exp))
      asc.test <- fisheSet1.test(x = matrix(c(n11,n01,n10,n00),nrow=2))
      p.val <- asc.test$p.value
      coexp.count <- n11
      Set1.count <- length(Set1.exp)
      Set2.count <- length(Set2.exp)
      temp.row <- cbind.data.frame(Combo,Set1.count,Set2.count,coexp.count,p.val,stringsAsFactors = F)      
      return(temp.row)
    }))
    return(out.row)
  }))
  invisible(gc())
}

gene.frame.coexp <- function(exp.tbl,gene.set.frame,ncores=4){
  exp.bool <- bool.data (exp.tbl = exp.tbl)
  out.frame <- do.call("rbind",mclapply(X = 1:nrow(gene.set.frame), FUN = function(x){
    
    Set1.exp <- names(which(exp.bool[as.character(gene.set.frame[x,1]),] == 1 ))
    Set2.exp <- names(which(exp.bool[as.character(gene.set.frame[x,2]),] == 1))
    Set1.not.exp <- names(which(exp.bool[as.character(gene.set.frame[x,1]),] == 0))
    Set2.not.exp <- names(which(exp.bool[as.character(gene.set.frame[x,2]),] == 0))
    n11 <- length(intersect(Set1.exp,Set2.exp))
    n01 <- length(intersect(Set1.not.exp,Set2.exp))
    n10 <- length(intersect(Set1.exp,Set2.not.exp))
    n00 <- length(intersect(Set1.not.exp,Set2.not.exp))
    asc.test <- fisher.test(x = matrix(c(n11,n01,n10,n00),nrow=2))
    p.val <- asc.test$p.value
    coexp.count <- n11
    Set1.count <- length(Set1.exp)
    Set2.count <- length(Set2.exp)
    temp.row <- cbind.data.frame(gene.set.frame[x,],Set1.count,Set2.count,coexp.count,p.val,stringsAsFactors = F)
    return(temp.row)
  },mc.cores = ncores))
  invisible(gc())
  return(out.frame)
}

#' @description The function calculates the coexpression of genes in each row of a given dataframe using the gene expression table.
#' @param exp.tbl Gene expression matrix to calculate coexpression of genes.
#' @param gene.frame A dataframe with genes. The coexpression is calculated for genes in each row.
#' @param ncores Number of cores to be used in the function. (Default : 4) 
#' @return gene.frame dataframe with an additional column corresponding to the number of cells coexpressing the genes. 

gene.coexp <- function(exp.tbl,gene.frame,ncores=4){
  gene.frame <- taRifx::remove.factors(gene.frame)
  row.names(gene.frame) <- 1:nrow(gene.frame)
  bool.exp.tbl <- bool.data(exp.tbl = exp.tbl)
  out <- do.call("rbind",lapply(X = 1:nrow(gene.frame), FUN = function(i){
    genes <- as.character(gene.frame[i,])
    len <- length(intersect(genes,row.names(bool.exp.tbl)))
    if (len == length(genes)) {
      bool.gene.tbl <- bool.exp.tbl[genes,]
      bool.sum <- colSums(bool.gene.tbl)
      coexp.count <- length(bool.sum[bool.sum == ncol(gene.frame)])
      out.row <- cbind.data.frame(gene.frame[i,],coexp.count,stringsAsFactors = F)
      return(out.row)
    }else{
      return(NULL)
    }
  }))
  return(out)
}


############################
#  SighotSpotter Functions
############################


#' General pipeline for SigHotSpotter
#'
#' The function computes compatibility scores for signaling intermediates
#'
#' @param species Currently supported species: "HUMAN", "MOUSE"
#' @param input_data File name for input gene expression data
#' @param cutoff Maximum number of zero-value genes, above this cutoff the genes are excluded
#' @param DE_Genes_data Differential expression dataset (1 for up-regulated, -1 for down-regulated genes)
#' @param percentile Predicted intermediates are taken into account above this threshold
#' @param invert_DE If the differential expression should be inverted, default = FALSE
#' @param showprogress shows progress bar in shiny app if set to TRUE, set it to FALSE in batch mode without GUI, default = TRUE
#' @return Compatibility scores
#' @export
SigHotSpotter_pipeline <- function(species, idata, cutoff, DE_Genes, percentile, invert_DE = FALSE, showprogress = TRUE,ncores=4){
  
  subg=Data_preprocessing(input_data = idata,cutoff = cutoff,species = species)
  
  ## Calculate stationary distribution of the MC
  Steady_state_true=Markov_chain_stationary_distribution(subg)
  
  prob.matrix <- summary(Steady_state_true$prob.matrix)
  
  edge.id.name <- cbind.data.frame(as_edgelist(graph = subg,names = T), as_edgelist(graph = subg,names = F),stringsAsFactors = F)
  colnames(edge.id.name) <- c("a","b","c","d")
  
  prob.matrix <- inner_join(x = edge.id.name,y = prob.matrix, by = c("c" = "i" , "d" = "j"))[,c(1,2,5)] 
  
  Steady_state_true <- Steady_state_true$SD
  
  ## Retrieves high probability intermediates
  int=high_probability_intermediates(x = Steady_state_true, intermediates = intermediates,percentile =  percentile)
  gintg=integrate_sig_TF(g = subg,x = Steady_state_true,deg = DE_Genes, non_interface_TFs = non_interface_TFs,TF_TF_interactions = TF_TF_interactions )
  
  # nTF=nonterminal_DE_TFs(g = gintg,deg = DE_Genes,non_interface_TFs = non_interface_TFs)
  
  target.TFs <- inner_join(x = DE_Genes, y = TF_TF_interactions, by = c("Gene" = "Source"))
  target.TFs <- target.TFs[which(target.TFs$Target %in% V(gintg)$name),]
  if(nrow(target.TFs) == 0){
    return(NULL)
  }
  target.coexp <- gene.coexp(exp.tbl = idata,gene.frame = target.TFs[,c(1,3)],ncores = ncores)
  target.coexp <- inner_join(x = target.coexp, y = target.TFs[,-2], by = c("Gene", "Target"))
  target.coexp$perc <- target.coexp$coexp.count/dim(idata)[2]
  target.coexp <- target.coexp[which(target.coexp$perc > 0.05 & target.coexp$Effect == 1),]
  iTF.target.info <- target.coexp[,c(1,2,4)]
  nTF <- iTF.target.info$Target
  
  nTF_scoring <- unique(nTF)
  names(nTF_scoring) <- nTF_scoring
  
  if (length(int) == 0){
    cat("No intermediates found. You may decrease the percentile in order to find intermediates.")
    return(NULL)
  }else if (length(nTF) == 0){
    cat("No non-terminal conserved TFs found at lower conservation also.")
    nTF <- DE_Genes
    return(NULL)
  }
  ## Computing compatibility scores
  score <- lapply(nTF_scoring,function(x){return(comp_score_tf(x,int,gintg))})
  score <- lapply(nTF,function(x){score[[x]]})
  if(is.null(score)){
    cat("No shortest path found. You may decrease the cutoff in order to find shortest path.")
    return(NULL)
  }else{
    #converting the nested list into a matrix whose row sum will give probability of each intermediate
    score_m=(matrix(unlist(score), ncol=length(score), byrow=F))
    score_m_means=as.list(rowMeans(score_m))
    final_score=compatability_score(score_m_means,Steady_state_true,int)
    
    toiintA <- as.character(final_score[which(final_score$Activation_probability > 0.5),]$Gene)
    toiintI <- as.character(final_score[which(final_score$Activation_probability < 0.5),]$Gene)
    
    
    #pruning the integrated networks
    gintg.p=prun.int.g(gintg)
    
    #building networks for all intermediates for active signaling hotspots
    sp_int_A <- lapply(toiintA,to_sp_net_int,gintg.p,nTF,DE_Genes,non_interface_TFs)
    
    #building networks for inactive signaling hotspots
    sp_int_I <- lapply(toiintI,to_sp_net_int,gintg.p,nTF,DE_Genes,non_interface_TFs)
    
    #retrieve receptors for active intermediates & inactive:
    u.gr <- Reduce(graph.union,sp_int_A)
    if(class(u.gr) == "igraph"){
      aa <- incident(u.gr,"NICHE",mode="out")
      bb <- igraph::ends(u.gr,aa)
      active_receptors <- bb[,2]
      active_receptors <- intersect(x = active_receptors,y = LR$Receptor)
    }else{
      active_receptors <- NULL
    }
    
    u.gr <- Reduce(graph.union,sp_int_I)
    if(class(u.gr) == "igraph"){
      aa <- incident(u.gr,"NICHE",mode="out")
      bb <- igraph::ends(u.gr,aa)
      inactive_receptors <- bb[,2]
      inactive_receptors <- intersect(x = inactive_receptors, y = LR$Receptor)
    }else{
      inactive_receptors <- NULL
    }
    
    
    ### Link signaling molecules to receptors
    del <- incident(gintg.p,"NICHE",mode = "out")
    gintg.p.noNiche <- delete.edges(gintg.p,del)
    dists <- distances(gintg.p.noNiche,to = as.character(final_score$Gene), v = as.character(unique(c(active_receptors,inactive_receptors))), mode = "out", weights = NA, algorithm = "johnson")
    dists <- melt(dists)
    dists <- dists[which(!is.infinite(dists$value)),]
    print("Weighted hotspot")
    edge_df <- ends(gintg.p.noNiche,E(gintg.p.noNiche))
    weights <- as.numeric(apply(idata[edge_df[,1],-1],1,mean))*as.numeric(apply(idata[edge_df[,2],-1],1,mean))
    weights[is.na(weights)] <- 0
    g_test <- set_edge_attr(gintg.p.noNiche,"weight", value = weights)
    dists_weighted <- distances(g_test,to = as.character(final_score$Gene), v = as.character(unique(c(active_receptors,inactive_receptors))), mode = "out", weights = NULL, algorithm = "johnson")
    dists_weighted <- melt(dists_weighted)
    dists_weighted <- dists_weighted[which(!is.infinite(dists_weighted$value)),]
    
    recs.to.perturb <- unique(c(active_receptors,inactive_receptors))
    
    cat("   Calculating shortest path weights\n")
    
    path.sums <- path.prob.sum(subg = subg,iTF = unique(iTF.target.info$Gene), Receptors = Receptors,prob.matrix = prob.matrix,ncores = ncores)
    
    # "perturbations" = as.data.frame(perturb.tests,stringsAsFactors = F)
    OUTPUT <- list("active"= active_receptors, "inactive"= inactive_receptors,"iTF.targets" = iTF.target.info,"final.score" = as.data.frame(final_score,stringsAsFactors = F),"path.sums" = path.sums, "rec.hotspot" = dists, "rec.hotspot.weighted" = dists_weighted)
    return(OUTPUT)
  }
}


Data_preprocessing <- function(input_data,cutoff,species){
  
  b=input_data
  ##COMVERT CELLS WITH LESS THAT 1 FPKM TO 0
  b[b < 1] <- 0
  ##Add a new gene Dummy with expression value 1, Only works if Dummy is present in the initial network
  b[nrow(b)+1, ] <- c(dummy.var, rep(1,(ncol(b)-1)))
  #This is to convert chr into numericversion
  b[2:ncol(b)]<-as.data.frame(lapply(b[2:ncol(b)],as.numeric,b[2:ncol(b)]))
  ##Renaming the first column for finding the union
  colnames(b)[1] <- "Source"
  a=Background_signaling_interactome
  
  ## Removing vertices connected to Dummy which are not receptors
  non.recs <- which(a$Source == dummy.var & !(a$Target %in% Receptors))
  recs <- setdiff(c(1:nrow(a)),non.recs)
  a <- a[recs,]
  
  
  ab=join(a,b,by=c("Source"),type="left",match="first")
  colnames(ab)[3:ncol(ab)]="source"
  colnames(b)[1] <- "Target"
  ab1=join(a,b,by=c("Target"),type="left",match="first")
  names(ab1) <- NULL
  names(ab) <- NULL
  ab=ab[,4:ncol(ab)]
  ab1=ab1[,4:ncol(ab1)]
  ab=as.matrix(ab)
  ab1=as.matrix(ab1)
  ########Elementwise product
  g=ab * ab1
  ########Sum of elementwise product
  sum_product_expression=rowSums(g, na.rm = FALSE, dims = 1)
  g3=cbind(a,sum_product_expression)
  ########Calculation of precentage of cells expressed
  h=rowSums(g != 0)
  percent_expressed=(h*100)/ncol(ab)
  g3=cbind(g3,percent_expressed)
  g3[is.na(g3)]<-0

  ######NETWORK preprocessing
  g <- graph.data.frame(as.data.frame(g3))
  del=E(g)[sum_product_expression==0|percent_expressed<as.numeric(cutoff)]
  g <- delete.edges(g,del)
  #SINCE THE TFs AND RECEPTORS ARE ALREADY CONNECTED TO DUMMY, REMOVE ANY NODE THAT HAS ZERO in degree or zero out degree
  #To ensure reachability for the Markov chain
  V(g)$degree=igraph::degree(g, v=V(g), mode = c("in"))
  #Select Nodes to be deleted
  del=V(g)[degree==0]
  #delete vertices from graph
  while(length(del)!=0)
  {
    g <- delete.vertices(g,del)
    V(g)$degree=igraph::degree(g, v=V(g), mode = c("in"))
    del=V(g)[degree==0]
  }
  #Same as above but remove nodes with with zero out degree
  V(g)$degree=igraph::degree(g, v=V(g), mode = c("out"))
  #Select Nodes to be deleted
  del=V(g)[degree==0]
  while(length(del)!=0)
  {
    g <- delete.vertices(g,del)
    V(g)$degree=igraph::degree(g, v=V(g), mode = c("out"))
    del=V(g)[degree==0]
  }
  #####TO EXTRACT THE LARGEST STRONGLY CONNECTED COMPONENT
  members <- membership(clusters(g, mode="strong"))
  SCC <- clusters(g, mode="strong")
  subg <- induced.subgraph(g, which(membership(SCC) == which.max(sizes(SCC))))
  subg=simplify(subg,edge.attr.comb=list("first"))
  subg
}

##STATIONARY DISTRIBUTION IN R ITSELF
##make a sparce-matrix from the edgelist
Markov_chain_stationary_distribution <- function(subg){ #The function takes the edgelist with probabilitys and computes the SS probability
  ####Write this subgraph as edgelist to make normal graph with ids ie.e names=F
  out <- list()
  transition_probability=as.data.frame(as_edgelist(subg,names = F),stringsAsFactors = FALSE)
  transition_probability$probability=paste(E(subg)$sum_product_expression)
  transition_probability[3]=as.numeric(transition_probability[[3]])
  myMatrix = sparseMatrix(i = transition_probability[1:nrow(transition_probability),1], j = transition_probability[1:nrow(transition_probability),2],x = transition_probability[1:nrow(transition_probability),3])
  #Making a stochastic matrix
  myMatrix = (myMatrix)/Matrix::rowSums((myMatrix))
  el=eigs(Matrix::t(myMatrix),1,which="LR")
  SD=(abs(el$vectors))/sum(abs(el$vectors))
  SD=as.data.frame(SD,stringsAsFactors=FALSE)
  SD
  SD=cbind((as.data.frame(V(subg)$name,stringsAsFactors = FALSE)),SD)
  SD=as.data.frame(SD,stringsAsFactors=FALSE)
  colnames(SD)[1] <- "Gene"
  out_SS=paste("Steady_state",sep="")
  colnames(SD)[2] <- out_SS
  out$SD <- SD
  out$prob.matrix <- myMatrix
  return(out)
}

#----------------------------------------------------------------------------------
# TO CALCULATE THE COMPATABILITY SCORE FOR THE HIGH PROBABILITY INTERMEDIATES
#----------------------------------------------------------------------------------

high_probability_intermediates <- function(x, intermediates, percentile)
{
  intermediates=unique(join(intermediates,x,by=c("Gene"),type="inner"))
  if(nrow(intermediates) == 0){
    return(c())
  }
  #Selecting top 90 percentile of intermediates with steady state probabilities
  percentile=as.numeric(percentile)/100
  SS_90percentile=as.vector(quantile(intermediates[,2], as.numeric(percentile)))
  ##Shortlisting high-probability intermediates > 90 percentile of SS probability
  int=as.vector((subset(intermediates, intermediates[2] > SS_90percentile , select=c("Gene")))$Gene)
  int
}


#Function for integrating signaling and TF networks, g=signaling graph, x steady state vector, deg=differentially expressed genes
integrate_sig_TF <- function(g,x,deg, non_interface_TFs, TF_TF_interactions ){
  
  el=as_edgelist(g)
  graph_markov=as.data.frame(cbind(el,E(g)$Effect))
  colnames(graph_markov)=c("Source","Target","Effect")
  colnames(deg)=c("Gene","Bool_DEG")
  non_interface_DE_TFs=join(deg,non_interface_TFs,by=c("Gene"),type="inner")
  #Get the phenotype from the non_interface_DE_TFs and map it to the TF-TF interaction network
  colnames(non_interface_DE_TFs)[1]="Target"
  DE_TF_TF_interactions_target=join(TF_TF_interactions,non_interface_DE_TFs,by=c("Target"),type="left")
  DE_TF_TF_interactions_target=na.omit(DE_TF_TF_interactions_target)
  ab=DE_TF_TF_interactions_target
  names(ab)<-NULL
  ab=as.matrix(ab)
  ab[,3]=as.numeric(ab[,3])*as.numeric(ab[,4])
  ab=as.data.frame(ab)
  names(ab)=c("Source","Target","Effect","DEG")
  graph_markov$Effect=as.numeric(as.character(graph_markov$Effect))
  graph_markov=rbind(graph_markov,ab[1:3]) #merging the nTF interaction with appropriate sign Effect with the original graph
  colnames(x)[1] <- "Source"
  ab=join(graph_markov,x,by=c("Source"),type="left",match="first")
  colnames(ab)[3:ncol(ab)]="source"
  colnames(x)[1] <- "Target"
  ab1=join(graph_markov,x,by=c("Target"),type="left",match="first")
  names(ab1) <- NULL
  names(ab) <- NULL
  ab=ab[,4:ncol(ab)]
  ab1=ab1[,4:ncol(ab1)]
  #creating node SS as the edge property
  weight=as.numeric(as.matrix(ab))
  #edge_P=as.data.frame(weight)
  graph_markov=(cbind(graph_markov,weight))
  graph_markov[is.na(graph_markov)] <- 1  #Making TF-TF interactions dependent only on the expression status
  graph_markov$Effect=as.numeric(as.matrix((graph_markov$Effect)))
  g3 <- graph.data.frame(as.data.frame(graph_markov))
  #updating the graph attribute for the adjacency matrix i.e. product SS (weight) and effect
  E(g3)$weight=E(g3)$weight*E(g3)$Effect
  #deleting TF nodes with no indegree
  V(g3)$degree=igraph::degree(g3, v=V(g3), mode = c("in"))
  #Select Nodes to be deleted
  del=V(g3)[degree==0]
  #delete vertices from graph
  while(length(del)!=0)
  {
    g3 <- delete.vertices(g3,del)
    V(g3)$degree=igraph::degree(g3, v=V(g3), mode = c("in"))
    del=V(g3)[degree==0]
  }
  g3
}



#Function of shortlisting non-terminal differentially expressed genes
nonterminal_DE_TFs <- function(g,deg,non_interface_TFs){
  colnames(deg)=c("Gene","Bool_DEG")
  non_interface_DE_TFs=join(deg,non_interface_TFs,by=c("Gene"),type="inner") #IF THIS IS ZERO NEED TO ABORT
  #load non-terminal TFs
  nTF=non_interface_DE_TFs[1]
  names(nTF)=NULL
  nTF=as.vector(t(nTF))
  nTF<-intersect(nTF,V(g)$name) #Some TFs must still be missing in the final g3
  # nTF
  if (length(nTF) == 0){
    cat("No downstream TF found for the cutoff employed. You may decrease the cutoff in order to find TFs.")
    return(NULL)
  }else{
    return(nTF)
  } 
}

#Function of shortlisting non-terminal differentially expressed genes But without stop command for building networks
#Function for classifying nTF as up or down regulated
up_down_tfs <- function(g,deg,non_interface_TFs)
{
  colnames(deg)=c("Gene","Bool_DEG")
  non_interface_DE_TFs=join(deg,non_interface_TFs,by=c("Gene"),type="inner") #IF THIS IS ZERO NEED TO ABORT
  #load non-terminal TFs
  nTF=non_interface_DE_TFs[1]
  names(nTF)=NULL
  nTF=as.vector(t(nTF))
  nTF<-intersect(nTF,V(g)$name) #Some TFs must still be missing in the final g3
  nTF
}

##Path_weight
#product_path_weight <- function(path, graph) sum(E(graph, path=path)$weight)/length(path)
product_path_weight <- function(path, graph){
  edge_weights=E(graph, path=path)$weight
  v_Edge_weight=sum(edge_weights)
  if (v_Edge_weight == 0) {
    return(0)
  } else{
    return(prod(edge_weights[edge_weights!=0]))
  }
}

path_weights <- function(path, graph) (E(graph, path=path)$weight)

##--Function for spliting and taking product of shortest path res file---
##This function takes ONE shortest path (x) and the adjacency matrix (l) as as input, and returns the product of SS probability of the intermediates in that shortest path.


##Function for Compatability score
##This function takes s="source" (one source gene), t="target" (a vector of target gene),g="graph", l="adjacency matrix" as input and finds the shortest paths and passes the argument to spsplit.
##Then it gets the product of intermediates of each path for a source and all its targets and returs its product as final output.

spcal_path_weight <- function(s,t,g){
  paths=(get.all.shortest.paths(g, s, t, mode = c("out"), weights=NA)$res)
  if (length(paths) == 0){
    return(0)
  } 
  
  weight=lapply(paths,product_path_weight,g)
  #s=skewness(unlist(weight))
  s=weight_probability(unlist(weight))
  #s=sum(unlist(weight))
  return(s)
}

#Function to parallelize the Compatability Score calculations
parallel.function.comp <- function(i){
  mclapply(i,spcal,nTF,g3,l)
}

#FUNCTION FOR SKEWNESS # FROM MOMENTS PACKAGE
skewness <- function (x, na.rm = FALSE)
{
  if (is.matrix(x))
    apply(x, 2, skewness, na.rm = na.rm)
  else if (is.vector(x)) {
    if (na.rm)
      x <- x[!is.na(x)]
    n <- length(x)
    (sum((x - mean(x))^3)/n)/(sum((x - mean(x))^2)/n)^(3/2)
  }
  else if (is.data.frame(x))
    sapply(x, skewness, na.rm = na.rm)
  else skewness(as.vector(x), na.rm = na.rm)
}

weight_probability <- function(x)
{
  x=unlist(x)
  x_pos=x[x>0]
  x_neg=x[x<0]
  x_tot=sum(abs(x_pos),abs(x_neg))
  #probability of the intermediate to be compatible: closer to 1 more compatible it is and closer t zero more incompatible it is
  p_pos=sum(x_pos)/x_tot
}

comp_score_tf <- function(t,s,g) #x is a list of comp score
{
  comp_score=lapply(s,spcal_path_weight,t,g)
}

compatability_score <- function(x,y,int) #where x is compatability score as list and y is steady state
{
  x=unlist(x)
  x=cbind(as.data.frame(int),as.data.frame(unlist(x)))
  colnames(x)=c("Gene", "Activation_probability")
  x=join(x,y,by=c("Gene"),type="inner")
  x=x[order(x$Activation_probability,decreasing = TRUE),]
}

path.prob.sum <- function(subg,iTF,Receptors,prob.matrix,ncores=4){

  Receptors <- intersect(Receptors,V(subg)$name)
  iTF <- intersect(iTF,V(subg)$name)
  
  if(length(iTF) == 0){
    cat("  No interface TFs in graph\n")
    return(NULL)
  }else{
    
    path.sum.frame <- do.call("rbind",mclapply(X = Receptors,FUN = function(rec){
      #print(rec)
      paths <- all_shortest_paths(graph = subg,from = rec, to = iTF,mode = "out",weights = NA)$res
      path.weights <- lapply(X = paths,FUN = function(path){
        path <- names(path)
        prod = 1
        if(length(path) > 1){
          for (j in 1:(length(path)-1)) {
            weight <- as.numeric(prob.matrix[which(prob.matrix$a == path[j] & prob.matrix$b == path[j+1]),][,3])
            prod = prod * weight
          }
        }
        names(prod) <- path[length(path)]
        return(prod)
      })
      names(path.weights) <- sapply(path.weights,function(x){names(x)})
      path.weights.sum <- aggregate(unlist(path.weights),by = list(Gene = names(path.weights)), FUN = function(x){sum(x,na.rm = TRUE)}, simplify = FALSE, drop = FALSE)
      path.weights.sum$Receptor <- rep(rec,nrow(path.weights.sum))
      path.weights.sum <- path.weights.sum[,c("Receptor","Gene","x")]
      colnames(path.weights.sum)[3] <- "path.weight.sum"
      path.weights.sum$path.weight.sum <- as.numeric(path.weights.sum$path.weight.sum)
      return(path.weights.sum)
    },mc.cores = ncores))
    vals  <- path.sum.frame$path.weight.sum
    path.sum.frame$z.score <- (vals-mean(vals))/sd(vals)
    # path.sum.frame <- path.sum.frame[which(path.sum.frame$z.score > 0),]
    return(path.sum.frame)
  }
}

bool.data.sc2i <- function(exp.tbl,expr.thrs=0){
  bool.tbl <- matrix(data = as.numeric(exp.tbl > expr.thrs),nrow = nrow(exp.tbl),ncol = ncol(exp.tbl))
  colnames(bool.tbl) <- colnames(exp.tbl)
  row.names(bool.tbl) <- row.names(exp.tbl)
  # bool.tbl <- bool.tbl[which(rowSums(bool.tbl) > 0),]
  return(bool.tbl)
}


SelectExpGene <- function(data,anno.tbl, exp=0.1,ncores=4) {
  
  all.pops <- unique(anno.tbl[,2])
  
  out.frame <- do.call("cbind",lapply(all.pops,function(celltype){
    #print(celltype)
    temp.data <- data[,which(colnames(data) == celltype)]
    #temp.mean <- rowMeans(x = temp.data,na.rm = F)
    temp.mean <- apply(temp.data,1,function(x){mean(x[x > 0])})
    temp.bool <- bool.data.sc2i(exp.tbl = temp.data)
    temp.expr <- rowMeans(temp.bool)
    out <- rep(0,nrow(temp.data))
    out[temp.expr >= exp] <- temp.mean[temp.expr >= exp]
    out <- as.data.frame(out,stringsAsFactors = F)
    colnames(out) <- celltype
    return(out)
  }))
  
  row.names(out.frame) <-  row.names(data)
  #invisible(gc())
  return(out.frame)
}

InteractStrength.minimal <- function(expressed_data,tissue.LR,ncores = 4){
  out <- apply(tissue.LR,1,function(x){
    expressed_data[x[2],x[1]] * expressed_data[x[4],x[5]]
  })
  return(out)
}

InteractStrength <- function(expressed_data,LR,ncores=4) {
  v <- as.data.frame(gtools::permutations(n=ncol(expressed_data),r=2,v=colnames(expressed_data),repeats.allowed=T))
  v <- tidyr::unite(v, "Combination", c("V1","V2"), remove = FALSE)
  
  LR <- LR[which(LR$Ligand %in% row.names(expressed_data) & LR$Receptor %in% row.names(expressed_data)),]
  
  out.frame <- do.call("cbind",lapply(X = 1:nrow(v),FUN = function(i){
    # print(i)
    out.col <- as.data.frame(expressed_data[LR$Ligand,as.character(v[i,2])]*expressed_data[LR$Receptor,as.character(v[i,3])])
    colnames(out.col) <- v[i,1]
    return(out.col)
  }))
  row.names(out.frame) <- LR[,1]
  colnames(out.frame) <- v[,1]
  invisible(gc())
  return(out.frame)
}

scoringFun <- function(data,tissue.R.TF,tissue.LR,LR,sig.cutoff = 0.9,z.score.cutoff = 0){
  pops <- union(tissue.LR$Lig.pop,tissue.LR$Rec.pop)
  out.tissue.lr.scored <- list()
  for(lig.pop in pops){
    for(rec.pop in pops){
      print(paste0(lig.pop,"_",rec.pop))
      lrsub <- tissue.LR[tissue.LR$Lig.pop == lig.pop & tissue.LR$Rec.pop == rec.pop,]
      if(nrow(lrsub) == 0){
        next
      }
      dat_lig <- data[unique(lrsub$Ligand),which(colnames(data) == lig.pop)]
      dat_rec <- data[unique(lrsub$Receptor),which(colnames(data) == rec.pop)]
      lig_means <- apply(dat_lig,1,function(x){mean(x[x>0])})
      rec_means <- apply(dat_rec,1,function(x){mean(x[x>0])})
      lrsub$score <- apply(lrsub,1,function(x){as.numeric(lig_means[x[2]])*as.numeric(rec_means[x[4]])})
      
      test_lig <- data[LR$Ligand,which(colnames(data) == lig.pop)]
      test_lig_means <- apply(test_lig,1,function(x){mean(as.numeric(x[x>0]))})
      test_rec <- data[LR$Receptor,which(colnames(data) == rec.pop)]
      test_rec_means <- apply(test_rec,1,function(x){mean(as.numeric(x[x>0]))})
      scores <- unname(test_lig_means)*unname(test_rec_means)
      scores[is.na(scores)] <- 0
      e <- ecdf(scores)
      lrsub$significance <- e(lrsub$score)
      out.tissue.lr.scored[[paste0(lig.pop,"_",rec.pop)]] <- lrsub
    }
  }
  out.tissue.lr.scored.joined <- do.call("rbind",out.tissue.lr.scored)
  
  m <- tissue.R.TF[tissue.R.TF$z.score >= z.score.cutoff,c(1,2)]
  m <- m[!duplicated(m),]
  mm <- merge(out.tissue.lr.scored.joined,m,by.x = c("Rec.pop","Receptor"), by.y = c("Celltype","Receptor"))
  mm <- mm[mm$significance >= sig.cutoff,]
  mm <- mm[,c(3,4,5,2,1,6,7,8)]
  
  return(mm)
}






