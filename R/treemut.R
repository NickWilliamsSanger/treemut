#' Constructs a genotype summary.
#'
#' The summary is represented as a data.frame where each row corresponds to a branch.
#' @param phylo  An APE tree
#'
#' @return a list with two elements: df - a data.frame where each row represents a branch; samples - a character vector of tip labels.
#' @export
reconstruct_genotype_summary=function(phylo
){
  dat=phylo$edge
  samples=phylo$tip.label
  N=length(samples)
  zeros=rep(0,N)
  profile=sapply(1:length(samples),function(i){tmp=zeros;tmp[i]=1;paste(tmp,collapse="")})
  df=data.frame(profile=profile,edge_length=phylo$edge.length[1:N])
  #Create empty list where each element will correspond to a node
  muts=lapply(1:dim(dat)[1],function(i) c())
  map_node=match(1:max(dat[,2]),dat[,2])
  #Here we loop through the tips (samples) and get all the ancestral nodes for each tip.
  #Then for each node in the list we add the tip to that nodes list of tips, ultimately giving us a list of samples that share the node. 
  for(i in 1:N){
    parents=get_ancestral_nodes(i,edge = dat,exclude_root=TRUE)
    for(j in parents){
      muts[[map_node[j]]]=append(muts[[map_node[j]]],i)
    }
  }
  profiles=sapply(muts,function(x){tmp=zeros;tmp[x]=1;paste(tmp,collapse="")})
  df=data.frame(profile=profiles,edge_length=phylo$edge.length,stringsAsFactors=FALSE)
  df=add_derived_profile_info(df,phylo$tip.label)
  list(df=df,samples=phylo$tip.label)
}

get_ancestral_nodes=function(node,edge,exclude_root=TRUE){
  idx=which(edge[,2]==node)
  parents=node ##Include the node
  while(length(idx)>0){
    if(length(idx)>1){
      stop("multiple parents!")
    }
    parent=edge[idx,1]
    parents=c(parents,parent)
    #This finds the parent of the current parent - thus navigating up to the root.
    idx=which(edge[,2]==parent)
  }
  if(exclude_root){
    parents[-length(parents)] ##The last node is the root.
  }else{
    parents
  }
}

add_derived_profile_info=function(profile_df,samples=sprintf("s%s",0:(nchar(profile_df$profile[1])-1))){
  base=rep(0,nchar(profile_df$profile[1]))
  samples_private=sapply(1:length(base),function(i){this=base;this[i]=1;paste(this,collapse="")})
  missing_samples=setdiff(samples_private,profile_df$profile)
  if(length(missing_samples)>0){
    profile_df=rbind(profile_df,data.frame(profile=missing_samples,edge_length=0))
  }
  
  profile_df$mut_count=nchar(profile_df$profile)-nchar(gsub("1","",profile_df$profile))
  profile_df$profile_int=lapply(profile_df$profile,split_profile)
  profile_df$id=1:dim(profile_df)[1]
  profile_df$label=profile_df$profile
  idx=which(profile_df$mut_count==1)
  profile_df$label[idx]=sapply(idx,function(i) samples[which(profile_df$profile_int[[i]]==1)])
  profile_df
}
#' Assigns mutations to a tree topology and sets the branch lengths based on the hard assignment of mutations to branches.
#'
#' @param tree  phylo. Should  already include the zeros outgroup..
#' @param mtr  Mutant read count matrix. Each row is a mutation and each colomn is a clonal sample BLAH
#' @param depth Depth count matrix. Each row is a mutation and each colomn is a clonal sample
#' @param error_rate  Error rate. This is intended as a coarse error measure and should not take more than 4 distinct values. Useful if a sample is mildly contaminated error->0.1.
#' @param maxits Number of EM iterations. 
#' @param VAF Expected VAF of mutant 
#' @param bverbose Whether to report progress
#' @return Returns a list containing:
#'  tree: the adjusted tree
#'  summary: a data.framethat is aligned with the input matrices and maps the mutations to branches:
#'   summary$edge_ml : Specifies the index of the branch in the tree edge matrix.  
#'   summary$pval :  A heuristic pvalue assessing the hypothesis that the mutation is consistent with the provided tree topology.   
#' @export
#' @useDynLib treemut 
assign_to_tree=function
(tree, 
 mtr,
 depth,##<<
 error_rate=rep(0.01,dim(mtr)[2]),##<<
 maxits=5,
 vaf=0.5,
 bverbose=TRUE){
  df=reconstruct_genotype_summary(tree)
  mtr=cbind(mtr,zeros=0)
  dep=cbind(depth,zeros=10)
  info=assign_to_df(mtr,dep,df,error_rate=rep(0.01,dim(mtr)[2]),maxits=maxits,vaf=vaf,bverbose=bverbose)
  tree$edge.length=info$df$df$edge_length  ##Could be expected edge_length
  info$summary$profile=info$df$df$profile[info$summary$edge_ml]
  info$tree=tree
  info
  ### Returns a list containing:
  ### tree: the adjusted tree
  ### and a data.frame "summary" that is aligned with the input matrices and maps the mutations to branches:
  ### summary$edge_ml : Specifies the index of the branch in the tree edge matrix.  
  ### summary$pval :  A heuristic pvalue assessing the hypothesis that the mutation is consistent with the provided tree topology.   
}

assign_to_df=function(mtr,dep,df,error_rate=rep(0.01,dim(mtr)[2]),maxits=5,vaf=0.5,bverbose=TRUE){
  if(!is.matrix(mtr)){
    nm=names(mtr)
    mtr=matrix(mtr,nrow=1)
    colnames(mtr)=nm
    dep=matrix(dep,nrow=1)
    colnames(dep)=nm
  }
  if(length(vaf)==1){
    vaf=rep(vaf,dim(mtr)[2])
  }
  p.err=error_rate#c(error_rate,1e-10)
  p.err=p.err[match(df$samples,colnames(mtr))]
  vaf=vaf[match(df$samples,colnames(mtr))]
  mtr=mtr[,df$samples]
  dep=dep[,df$samples]
  if(!is.matrix(mtr)){
    nm=names(mtr)
    mtr=matrix(mtr,nrow=1)
    colnames(mtr)=nm
    dep=matrix(dep,nrow=1)
    colnames(dep)=nm
  }
  tree_genotypes=do.call("rbind",df$df$profile_int)
  if(maxits==1){
    if(bverbose)
      cat("maxits=1 - therefore using provided edge lengths as initial edge lengths\n")
    el=df$df$edge_length
  }else{
    if(bverbose)
      cat("Initialising edge lengths to 1\n")
    el=rep(1,length(df$df$edge_length))
  }
  el=ifelse(el<1e-100,1e-100,el)
  loglik=rep(NA,maxits)
  n=dim(mtr)[1]
  for(i in 1:maxits){
    ol=el
    lik=get_likelihood_mtr_C(mtr,dep,tree_genotypes,el,p.error = p.err,vaf=vaf,bverbose=bverbose)
    edge_ml=apply(lik,1,which.max)
    #n=length(edge_ml)
    p=exp(lik-lik[(edge_ml-1)*n+1:n])##Substract max to just control range of lik... comes out in the wash later.
    p=p/rowSums(p)
    loglik[i]=sum(p*lik)
    if(is.na(loglik[i])){
      browser()
    }
    #browser()
    el=colSums(p,na.rm=T)
    ## Don't allow zero length el
    el=ifelse(el<1e-100,1e-100,el)
    epsilon=sum(abs(el-ol))/dim(mtr)[1]
    if(bverbose){
    cat("delta edge length=",epsilon,"\n")
    cat("Loglik=",loglik[i],"\n")
    }
    if(epsilon<1e-6){
      break
    }
  }
  
  df$df$expected_edge_length=el
  df$df$edge_length=sapply(1:length(df$df$edge_length),function(i){length(which(edge_ml==i))})
  
  p_else_where=1-p[(edge_ml-1)*n+1:n]
  if(bverbose){
    cat("Finished assigning mutations\ncalculating pvalues\n")
  }
  pval=rep(NA,length(edge_ml))
  for(i in 1:length(pval)){
    pval[i]= get_mutation_assignment_pval(df$df,edge_ml[i],mtr[i,],dep[i,],p.err,vaf)
    if(i %% 1000 == 0 && bverbose){
      cat("On",i," of ",length(pval),"\n")
    }
  }
  list(df=df,lik=lik,summary=data.frame(edge_ml=edge_ml,pval=pval,p_else_where=p_else_where),p=p)
}



get_likelihood_mtr_C=function(mtr,depth,geno,el,p.error=rep(0.01,dim(mtr)[2]),vaf=rep(0.5,dim(mtr)[2]),bverbose=FALSE){
  if(dim(mtr)[2]!=dim(geno)[2]){
    stop("error: dimension mismatch betwee tree genotypes and mtr")
  }
  nmuts=dim(mtr)[1]
  nsamp=dim(mtr)[2]
  nbranch=dim(geno)[1]
  res = .C("likelihood_mtr", 
           mtr=as.integer(mtr),
           depth=as.integer(depth),
           geno=as.integer(geno),
           el=as.double(el),
           p.error=as.double(p.error),
           nmuts=as.integer(nmuts),
           nsamp=as.integer(nsamp),
           nbranch=as.integer(nbranch),
           vaf=as.double(vaf),
           lik=double(nmuts*nbranch),
           bverbose=as.integer(bverbose)
           
  )
  matrix(res[["lik"]],ncol=nbranch)
}

get_mutation_assignment_pval=function(df,i,mtr,dep,p.error,vaf=0.5){
  profile_int=df$profile_int[[i]]
  idx.inside=which(profile_int==1)
  idx.outside=which(profile_int==0)
  V=vaf*(1-2*p.error)+p.error
  
  ##Calculate probability for observing <= total MTR at variant sites give VAF=0.5.  This is complicated by allowing different error rates for each sample (otherwise just pbinom)
  if(length(idx.inside)>0){
    ##We breakdown errors into unique categories to make calculation more tractable - could further discretise.
    v=unique(V[idx.inside])
    if(length(v)>4){
      stop("Too many error categories in get_pval: not tractable!")
    }
    mtri=mtr[idx.inside]
    depi=dep[idx.inside]
    vi=V[idx.inside]
    mtrii=sapply(v,function(v) sum(mtri[vi==v]))
    depii=sapply(v,function(v) sum(depi[vi==v]))
    p1=cpv(sum(mtrii),depii,v,lower.tail = TRUE)
    d1=1
  }else{
    p1=1
    d1=0
  }
  ##Calculate probability for observing >= total MTR at non-variant sites give VAF=0.  This is complicated by allowing different error rates for each sample (otherwise just pbinom)
  if(length(idx.outside)>0){
    v=unique(p.error[idx.outside])
    if(length(v)>4){
      stop("Too many error categories in get_pval: not tractable!")
    }
    mtri=mtr[idx.outside]
    depi=dep[idx.outside]
    vi=p.error[idx.outside]
    mtrii=sapply(v,function(v) sum(mtri[vi==v]))
    depii=sapply(v,function(v) sum(depi[vi==v]))
    p2=cpv(sum(mtrii),depii,v,lower.tail = FALSE)
    d2=1
  }else{
    p2=1
    d2=0
  }
  ##We are interested in whether either of the tests fails so we combine using conservative bonferroni.  
  ##We are performing 2 tests.  same as p.adjust(c(p1,p2),method = "bonferroni")
  pv=min((d1+d2)*min(c(p1,p2)),1)
  pv
}


cpv=function(mtrtot,depth,probs,lower.tail){
  ##Calculates the probability that there are a total of MTR mutant reads across N bins each with depth d_i and mutant prob p_i
  ##Does this by explicitly calculating  P(M | {d_i},{p_i)=sum_m P(mtr_i=m | d_i, p_i)*P(M-m | )
  n=length(depth)
  idx=order(depth)
  probs=probs[idx]
  depth=depth[idx]
  flag=ifelse(lower.tail,1,0)
  resk=list(mtrtot,probs,depth,n,lower.tail)
  res=.C("cumulate_binomial_pval",as.integer(round(mtrtot)),as.integer(round(depth)),as.double(probs),as.integer(n),as.integer(lower.tail),p=double(1))
  res[["p"]]
}

split_profile=function(profile){
  as.integer(unlist(strsplit(profile,split="")))
}

#' Generate a random tree 
#'
#' Generates a random tree with the specified number of tips. For a more realistic tree use e.g. rsimpop
#' @param nsamples  Number of tips
#'
#' @return phylo object
#' @export
generate_random_tree=function(nsamples){
  tree=rtree(nsamples)
  tree=bind.tree(tree,read.tree(text="(zeros:0);"))
  tree$edge.length=ceiling(200*tree$edge.length) ##Scales the #Muts per edge to the scale 0-200
  tree
}

#' Simulates read counts from the specfified genotype summary using a binomial model
#'
#' @param df  list. genotype summary
#' @param avgdeph  Average depth
#' @param depth Depth count matrix. Each row is a mutation and each column is a clonal sample
#' @param n_artifacts Number of artifact to include
#' @param p.error Error rate
#' @param vaf VAF for determining mutant read probability.
#' @return Returns a list containing:(mtr=mtr,depth=depth,edge=ml,df=df,p.error=p.error
#' @export
simulate_reads_from_tree=function(df,avgdepth,n_artifacts=0,p.error=0.01,vaf=0.5){
  if(length(p.error)>1){
    stop("p.error must be scalar here")
  }
  nsamples=length(df$samples)-1
  N=length(df$df$edge_length)
  dat=lapply(1:N,function(i) get_simulated_reads(df$df$profile_int[[i]],df$df$edge_length[i],avgdepth,p.error=p.error,vaf=vaf))
  
  ml=do.call("c",lapply(1:N,function(i) rep(i,df$df$edge_length[i])))
  if(n_artifacts>0){
    ##Add in some artefacts that aren't simulated according to the tree topology.
    NART=n_artifacts
    dat2=lapply(1:NART,function(i) get_simulated_reads(ifelse(runif(nsamples+1)<0.3,1,0),1,avgdepth,p.error=p.error,vaf = vaf))
    dat=c(dat,dat2)
    ml=c(ml,rep(NA,NART))
  }
  mtr=do.call("rbind",lapply(dat,function(x) x$mtr))
  depth=do.call("rbind",lapply(dat,function(x) x$depth))
  colnames(mtr)=df$samples
  colnames(depth)=df$samples
  p.error=rep(p.error,length(df$samples))
  idx.zeros=match("zeros",df$samples)
  mtr[,idx.zeros]=0
  depth[,idx.zeros]=100
  p.error[idx.zeros]=1e-6 ##This will force the inference code to keep the zeros branch as a zero length outgroup. 
  list(mtr=mtr,depth=depth,edge=ml,df=df,p.error=p.error)
}

get_simulated_reads=function(geno,n,depth,p.error=0.01,vaf=0.5){
  if(n==0){
    emptymat=matrix(1,nrow=0,ncol=length(geno))
    return(list(mtr=emptymat,depth=emptymat))
  }
  geno=rep(geno,n)
  depth=rpois(length(geno),depth)
  mtr=rbinom(length(depth),depth,prob=ifelse(geno==0,p.error,vaf*(1-2*p.error)+p.error))
  if(any(is.na(mtr))){
    browser()
  }
  list(mtr=matrix(mtr,nrow=n,byrow = TRUE),
       depth=matrix(depth,nrow=n,byrow = TRUE)
  )
}
#' Wrapper to run simulation 
#'
#' @export
run_sim=function(nsamp,depth,n_artifacts=0){
  tree=generate_random_tree(nsamp)
  df=reconstruct_genotype_summary(tree)
  dat=simulate_reads_from_tree(df,depth,n_artifacts = n_artifacts )
  res=assign_to_tree(tree,dat$mtr,dat$depth,error_rate=dat$p.error)##
  list(edge_length_orig=dat$df$df$edge_length,
       edge_length_inferred=res$df$df$edge_length,
       expected_edge_length_inferred=res$df$df$expected_edge_length,
       edge_idx_orig=dat$edge,
       edge_idx_ml=res$summary$edge_ml,dat=dat,res=res,tree=tree)
}
#' Plots simulation results
#' @export
plot_sim_result=function(sim){
  par(mfrow=c(2,2))
  mismatch_prop=length(which(sim$edge_idx_orig!=sim$edge_idx_ml))/length(sim$edge_idx_orig)
  plot(sim$edge_length_orig,sim$edge_length_inferred,xlab="Original Edge Length",ylab="Inferred Edge Length",
       main=sprintf("Hard Assigned edge length vs Orig (sd=%3.2f)\n Mismatch Proportion=%5.4f",
                    sd(sim$edge_length_orig-sim$edge_length_inferred),mismatch_prop)
  )
  abline(a=0,b=1)
  plot(sim$edge_length_orig,sim$expected_edge_length_inferred,xlab="Original Edge Length",ylab="Expected Edge Length",
       main=sprintf("Expected edge length (Soft Assigned) \nvs Orig (sd=%3.2f)",sd(sim$edge_length_orig-sim$expected_edge_length_inferred))
  )
  abline(a=0,b=1)
  simdf=data.frame(deviation=sim$expected_edge_length_inferred-sim$edge_length_orig,
                   edge_length_orig=sim$edge_length_orig)
  simdf=simdf[order(simdf$edge_length_orig),]
  
  with(simdf,plot(edge_length_orig,deviation,
       xlab="Original Edge Length",ylab="Deviation",
       main="Deviation vs Original Edge Length"))
  loess_fit <- with(simdf,loess(deviation ~ edge_length_orig))
  pl=predict(loess_fit,se = T)
 
  plu=pl$fit+1.96*pl$se
  pll=pl$fit-1.96*pl$se
  
  
  shade_between(simdf$edge_length_orig,plu,pll,adjustcolor("blue",0.3))
  with(simdf,lines(edge_length_orig, pl$fit, col = "blue",lwd=2))
  
  plot(simdf$edge_length_orig,abs(simdf$deviation),
       xlab="Original Edge Length",ylab="Absolute Deviation",
       main="Absolute Deviation vs Original Edge Length")
  
  
  
}

shade_between=function(x,y1,y2,color){
  polygon(c(x, rev(x), x[1]), c(y1, rev(y2), y1[1]), 
          col = color,border = color) 
}


compare_sim=function(df,res){
  res=data.frame(edge_length_orig=df$df$edge_length,
  edge_length_inferred=res$df$df$edge_length,
  expected_edge_length_inferred=res$df$df$expected_edge_length)
  res$hard_diff=res$edge_length_inferred-res$edge_length_orig
  res$soft_diff=res$expected_edge_length_inferred-res$edge_length_orig
  c(sdhard=sd(res$hard_diff),
    sdsoft=sd(res$soft_diff),
    adhard=mean(abs(res$hard_diff)),
    adsoft=mean(abs(res$soft_diff)))
  
}


