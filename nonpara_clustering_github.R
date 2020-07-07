Nonpara_clustering<-function(CURVES,q,nknot,range.grid, semimetric="deriv", threshold=0.1, nb.bw=100,
                             nss=0, mspg=0, centrality="mean",kind.of.kernel="quadratic"){
  ############get the number of curves in sample
  #CURVES<-subcurve
  #
  #CURVES<-subcurve2
  #CURVES<-sc1
  #CURVES<-ff1
  #CURVES<-group11
  #CURVES<-group12
  #CURVES<-group2
  #CURVES<-pre_sug
  ###################
  #"CURVES" contains the curves (row by row)
  # q: the order of derivative
  # nknot: the number of knots for b-spline
  # "range.grid" vector of length 2 containing the range of the grid at 
  #                 which the curve are evaluated (i.e. range of the 
  #                 discretization)
  # "kind.of.kernel" gives the kernel used for computing the modes
  # "semimetric" contains the semimetric type (default="deriv")
  # "threshold"  contains the value which gives the minimum gain
  #    accepted in terms of splitting score
  # "nb.bw"      gives the number of bandwidths for building the
  #    corresponding sequence (default="100")
  # "nss"        gives the number of subsamples required for
  #    computing the subsampling heterogeneity index (default="50")
  #    (if "nss=0", the method uses HI instead of SHI which reduces
  #    the computational cost)
  # "mspg"       gives the minimal size per group authorized 
  # "centrality": character string giving the centrality feature 
  #     compared with the mode (centrality="mean" or "median") 
  #     in order to compute the heterogeneity index (default="mean")      
  # Returns a list:
  #  MEANS     matrix containing the mean curve of each group
  #  MEDIANS   matrix containing the median curve of each group
  #  MODES     matrix containing the modal curve of each group
  #  Bw.opt    vector containing the optimal bandwidth of each group
  n <- nrow(CURVES)
  
  ############compute the distance between curves according to selected semimetric
  #semimetric<-"deriv"
  sm <- get(paste("semimetric.", semimetric, sep = ""))
  ####For semimetric.derive function, q is the order of derivative, range.grid is
  #c<-list(q,nknot,range.grid)
  #return(sm)
  SEMIMETRIC <- sm(CURVES,CURVES,nknot, q,range.grid)
  #return(SEMIMETRIC)
  #####Build the sequence of bandwidht
  #nb.bw<-100
  Hrange<- range(SEMIMETRIC)
  Bw.seq <- seq(Hrange[1], Hrange[2]*0.5 , length = nb.bw + 1)[-1]
  #####Compute the probability curves
  PROBCURVES <- matrix(0, n, nb.bw)
  for(i in 1:n){
    PROBCURVES[i,] <- prob.curve(i, SEMIMETRIC, Bw.seq)
  }
  #############Select the optimal bandwidth based on minimizing entropy
  res.classif.bw <- classif.bw(PROBCURVES, Bw.seq, 1:n)
  Bw.opt <- res.classif.bw$bw
  index <- res.classif.bw$index
  ##########################Compute the heterogeneity index of sample
  #kind.of.kernel<-'quadratic'
  if(nss==0){
    shi <- classif.hi(CURVES, SEMIMETRIC, kind.of.kernel='quadratic', Bw.opt, 1:n, semimetric, centrality="mean", q=2,nknot=15, range.grid=c(0,8))
  } else {
    shi <- classif.shi(CURVES, SEMIMETRIC, kind.of.kernel, Bw.opt, 1:n, nss, semimetric, centrality, q,nknot, range.grid)
  }
  
  #####################
  #Partitioning the full sample
  #Split the sample
  Group<-1:n
  Prob.points <- PROBCURVES[Group, index]
  pp.start <- min(Prob.points)
  pp.end <- max(Prob.points)
  subgroup <- list()
  if(pp.end > pp.start){
    est <- density(Prob.points, bw="bcv", from=pp.start, to=pp.end)
    Prob.lim <- est$x[rank.minima(est$y)]
    if(length(Prob.lim)==0){ 
      warning("Zero local minimum ==> this group is a terminal leaf of the classification tree")
    } else {
      kmax <- length(Prob.lim)
      subgroup[[1]] <- Group[Prob.points>Prob.lim[kmax]]
      subgroup[[kmax+1]] <- Group[Prob.points<= Prob.lim[1]]
      if(kmax>1){
        for(k in 1:(kmax-1))
          subgroup[[k+1]] <- Group[(Prob.points>Prob.lim[k]) & (Prob.points<=Prob.lim[k+1])]
      }
    }
  } else {
    remain<-length(Group)
    semi<-SEMIMETRIC[Group,Group]
    m<-1
    Group1<-Group
    while(remain != 0){
      sub<-Group[!is.na(Group[semi[Group,Group] < Bw.opt])]
      Group<-Group[-sub]
      #semi<-semi[Group,Group]
      remain<-remain-length(sub)
      subgroup[[m]]<-sub
      m<-m+1
      
    }
  }
  res.classif.split<-subgroup
  #return(res.classif.split)
  #res.classif.split <- classif.split(PROBCURVES, index, Group=1:n)
  
  ######The number of subgroups being splitted from sample
  nbgroups <- length(res.classif.split)
  res <- list(Groups=list())
  if(nbgroups!=0){
    #get the minimum size of subgroups
    minsize <- min(sapply(res.classif.split, length))
    mspg<-0
    if(minsize>mspg){
      ###
      Bw.opt <- 0
      Index <- 0
      shisub <- 0
      for(k in 1:nbgroups){
        k<-2
        Subgroup <- res.classif.split[[k]]
        subcurve<-CURVES[Subgroup,]
        Subsemi<-SEMIMETRIC[Subgroup,Subgroup]
        subBw.seq<-seq(range(Subsemi)[1], range(Subsemi)[2]*0.5, length=nb.bw+1)[-1]
        subPROBCURVES <- matrix(0, length(Subgroup), nb.bw)
        for(i in 1:length(Subgroup)){
          i<-1
          subPROBCURVES[i,] <- prob.curve(i, Subsemi, subBw.seq)
        }
        res.classif.bw <- classif.bw(subPROBCURVES, subBw.seq, 1:length(Subgroup))
        Bw.opt[k] <- res.classif.bw$bw
        Index[k] <- res.classif.bw$index
        nss<-0
        if(nss==0){
          shisub[k] <- classif.hi(subcurve, Subsemi, kind.of.kernel="quadratic", Bw.opt[k], 1:length(Subgroup), semimetric, centrality="mean", q=2,nknot=15,range.grid=c(0,8))
        } else {
          shisub[k] <- classif.shi(CURVES, SEMIMETRIC, kind.of.kernel, Bw.opt[k], Subgroup, nss, semimetric, centrality, ...)  
        }
      }
      split.score <- 1- sum(sapply(res.classif.split, length)*shisub)/(length(1:n)*shi)
      threshold<-0.1
      if(split.score>threshold){
        res <- list(Groups=res.classif.split, Bw.opt=Bw.opt, Index=Index, SHI=shisub, Ssc=split.score)
      }else {
        res <- list(Groups=res.classif.split, Bw.opt=Bw.opt, Index=Index, SHI=shisub, Ssc=split.score)
        return(res)
        cat("Do not split this sample")
        #return(split.score)
      }
    }
  }
  return(res)
}
