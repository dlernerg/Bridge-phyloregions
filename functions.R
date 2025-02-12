#functions
is.integer0 <- function(x)
{
  is.integer(x) && length(x) == 0L
}

tip.drop.fun <- function(GBIF_overlap,i){
  keep <- GBIF_overlap[i] %>% unlist()
  if (is.integer0(keep.com)){
    keep.com = 0
  }
  drop <- which(!(all %in% keep))
  tip.drop <- unique(GBIFpolygon_groups$specie[drop] %>% stringr::word(.,1))
  return(tip.drop)
}

make.phylo.rich <- function(globe.sample, GBIF_overlap,tree){
  
  dis.mat <- matrix(nrow = length(globe.sample), ncol = length(globe.sample))
  for (i in 1:length(globe.sample)){
    
    keep.i <- GBIF_overlap[globe.sample[i]] %>% unlist()
    
    if (length(keep.i)<2){
      next
    }else{
      for (j in 1:length(globe.sample)){
        
        keep.j <- GBIF_overlap[globe.sample[j]] %>% unlist()
        keep.com1 <- keep.i[which((keep.i %in% keep.j))]
        keep.com2 <- keep.j[which((keep.j %in% keep.i))]
        keep.com <- unique(c(keep.com1,keep.com2))
        
        if (length(keep.j)<2){
          next
        }else{
          if (is.integer0(keep.com)){
            keep.com = 0
          }
          
          
          tip.i <- unique(GBIFpolygon_groups$specie[keep.i] %>% stringr::word(.,1))
          tip.drop.i <- genus[which(!(genus %in% tip.i))]
          tip.j <- unique(GBIFpolygon_groups$specie[keep.j] %>% stringr::word(.,1))
          tip.drop.j <- genus[which(!(genus %in% tip.j))]
          
          tip.com <- unique(GBIFpolygon_groups$specie[keep.com] %>% stringr::word(.,1))
          tip.drop.com <- genus[which(!(genus %in% tip.com))]
          
          phy.new.i <- drop.tip(tree, tip.drop.i)
          phy.new.j <- drop.tip(tree, tip.drop.j)
          phy.new.com <- drop.tip(tree, tip.drop.com)
          
          B.i <- sum(phy.new.i$edge.length)
          B.j <- sum(phy.new.j$edge.length)
          B.ij <- sum(phy.new.com$edge.length)
          
          dis.mat[i,j] <- 1- B.ij/(0.5*(B.i+B.j))
        }
        
      }
    }
    
  }
  return(dis.mat)
}


delete.na <- function(DF, n=0) {
  DF[which(rowSums(DF, na.rm = T) != n),which(rowSums(DF, na.rm = T) != n)]
}

Test.Kmult<-function(x,phy,iter,plot.hist){
  library(ape)
  Kmult<-function(x,phy){
    x<-as.matrix(x)
    N<-length(phy$tip.label)
    ones<-array(1,N)
    C<-vcv.phylo(phy)
    C<-C[row.names(x),row.names(x)]
    a.obs<-colSums(solve(C))%*%x/sum(solve(C))
    #evol.vcv code
    distmat<-as.matrix(dist(rbind(as.matrix(x),a.obs)))
    MSEobs.d <- sum(distmat[(1:N),(N+1)]^2)
    #sum distances root vs. tips
    eigC <- eigen(C)
    D.mat<-solve(eigC$vectors
                 %*% diag(sqrt(eigC$values))
                 %*% t(eigC$vectors))
    dist.adj<-as.matrix(dist(rbind((D.mat
                                    %*%(x-(ones%*%a.obs))),0)))
    MSE.d<-sum(dist.adj[(1:N),(N+1)]^2)
    #sum distances for transformed data)
    K.denom<-(sum(diag(C))-
                N*solve(t(ones)%*%solve(C)%*%ones)) / (N-1)
    K.stat<-(MSEobs.d/MSE.d)/K.denom
    return(K.stat)
  }
  K.obs<-Kmult(x,phy)
  P.val <- 1
  K.val <- rep(0, iter)
  for (i in 1:iter){
    x.r<-as.matrix(x[sample(nrow(x)),])
    rownames(x.r)<-rownames(x)
    K.rand<-Kmult(x.r,phy)
    P.val<-ifelse(K.rand>=K.obs, P.val+1,P.val)
    K.val[i] <- K.rand
  }
  P.val <- P.val/(iter + 1)
  K.val[iter + 1] = K.obs
  if (plot.hist == T){
    hist(K.val, 30, freq = TRUE, col = "gray",
         xlab = "Phylogenetic Signal")
    arrows(K.obs, 50, K.obs, 5, length = 0.1, lwd = 2)
    
  }
  return(list(phy.signal = K.obs, pvalue = P.val))
}


find.num <- function(cells.num){
  num <- c()
  for (i in 1:length(cells.num)){
    
    GBIFo <- GBIF_overlap[cells.num[i]] %>% unlist()
    keep.i <- GBIFo
    spec.df1 <- unique(GBIFpolygon_groups$specie[keep.i] %>% stringr::word(.,1)) %>% as.data.frame()
    
    num[i] <- nrow(spec.df1)
    
  }
  return(num)
}

find.spec <- function(clus.globe){
  spec.df <- data.frame()
  for (i in 1:nrow(clus.globe)){
    
    GBIFo <- GBIF_overlap[clus.globe$cells[i]] %>% unlist()
    keep.i <- GBIFo
    spec.df1 <- unique(GBIFpolygon_groups$specie[keep.i] %>% stringr::word(.,1)) %>% as.data.frame()
    spec.df1$cluster <- clus.globe$cluster[i]
    spec.df <- bind_rows(spec.df,spec.df1)
  }
  return(spec.df)
}

#sort the different clusters given the expected as a function of the frequency of each cluster in the dataset
#has been changed
residues.fun <- function(df){
  
  ID <- which(!(is.na(df)))
  ID.na <- which((is.na(df)))
  clus.expected <- (cluster.nm$n / sum(cluster.nm$n[ID]))*sum(df, na.rm = T)
  clus.expected[ID.na] <- NA
  z <- (df/clus.expected)
  return(z)
}

residues.fun2 <- function(df){
  df.n <- sum(df, na.rm = T)
  z <- (df/df.n)
  return(z)
}

residues.fun_expected.corrected <- function(df) {
  
  total_counts_per_species <- rowSums(df, na.rm = TRUE)
  
  # Total count of species in each region
  total_counts_per_region <- colSums(df, na.rm = TRUE)
  
  # Expected count for each species in each region
  expected_counts <- sweep(df, 1, total_counts_per_species, "/")
  expected_counts <- sweep(expected_counts, 2, total_counts_per_region, "*")
  
  # Scale expected counts so their sum matches the observed total in each region
  scale_factors <- total_counts_per_region / colSums(expected_counts, na.rm = TRUE)
  expected_counts <- sweep(expected_counts, 2, scale_factors, "*")
  
  # Observed/Expected fractions
  obs_exp_fractions <- (df / expected_counts)
  
}


test.kmult.change.k <- function(spec.df.norm.t,phy){
  spec.df.norm.t[which(is.na(spec.df.norm.t))] = 0
  
  match.df <- spec.df.norm.t[match(phy$tip.label, rownames(spec.df.norm.t)),]
  isna <- which(is.na(match(phy$tip.label, rownames(spec.df.norm.t))))
  if(is.integer0(isna)){
    match.df <- match.df
  }else{
    phy <- drop.tip(phy, tip = phy$tip.label[isna])
    match.df <- match.df[-isna,]
  }
  
  test<- Test.Kmult(match.df,phy, iter = 500,plot.hist = F)
  
  return(test)
}

avg.df.forplot <- function(df, column_name){
  
  mpd.avg <- data.frame()
  mpd.avg[1:5,1] <- c("Trop","Bridge","Desert","Temp","Bor")
  
  mpd.avg$V2 <- NA
  mpd.avg$V2[2] <- mean(c(df[[column_name]][which(df$cluster == 7)], df[[column_name]][which(df$cluster == 8)]), na.rm = TRUE)
  mpd.avg$V2[4] <- mean(c(df[[column_name]][which(df$cluster == 5)], df[[column_name]][which(df$cluster == 4)], df[[column_name]][which(df$cluster == 3)]), na.rm = TRUE)
  mpd.avg$V2[5] <- mean(c(df[[column_name]][which(df$cluster == 1)], df[[column_name]][which(df$cluster == 2)]), na.rm = TRUE)
  mpd.avg$V2[3] <- df[[column_name]][which(df$cluster == 6)]
  mpd.avg$V2[1] <- df[[column_name]][which(df$cluster == 9)]
  
  colnames(mpd.avg) <- c("cluster","value")
  
  if (column_name == "mpd.obs.z"){
    
    mpd.avg$significance <- ifelse(mpd.avg$value < -1.64, "clustered", ifelse(mpd.avg$value> 1.64,"dispersed", "NS"))
    mpd.avg$cluster <- factor(mpd.avg$cluster, levels = rev(mpd.avg$cluster))
    
    mpd.plot2 <- ggplot(mpd.avg, aes(cluster,value))+
      geom_point(aes(col = significance)) + 
      geom_hline(yintercept = -1.64)+
      geom_hline(yintercept = 1.64)+
      scale_color_manual(values = c("red","black","grey"))+
      #scale_x_continuous(labels = c(1:k), breaks = c(1:k)) +
      theme_classic()+
      theme(legend.position = "none")  
  } else {
    
    mpd.avg$cluster <- factor(mpd.avg$cluster, levels = rev(mpd.avg$cluster))
    mpd.plot2 <- ggplot(mpd.avg, aes(cluster, value)) +
      geom_point() + 
      scale_y_continuous(trans = "log2") +  # Apply log2 transformation to the y-axis
      theme_classic() +
      theme(legend.position = "none")
  }
  
  
}
