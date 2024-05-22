rm(list=ls())
library('ape')
library('stats4')
cfoutprefix<-'out' ####THIS IS THE PREFIX OF THE CLONALFRAMEML OUTPUT FILES
tree<-read.tree(sprinf('%s.labelled_tree.newick',cfoutprefix))
seqs<-read.dna(sprinf('%s.ML_sequence.fasta',cfoutprefix),format='fasta')
mapping<-scan(sprinf('%s.position_cross_reference.txt',cfoutprefix),sep=',')
labs<-labels(seqs)
ma=max(mapping)
map2=vector("list", ma) 
for (i in 1:ma) {
  map2[[i]]=which(mapping==i)
}

#Transform edge matrix
edges<-tree$edge
treelabs<-c(tree$tip.label,tree$node.label)
for (i in 1:nrow(edges)) for (j in 1:ncol(edges)) edges[i,j]=which(labs==treelabs[edges[i,j]])

#Find branches with mutations
l<-length(seqs[1,])
muts<-matrix(F,l,nrow(edges))
for (i in 1:l) for (b in 1:nrow(edges)) 
  if (seqs[edges[b,1],i]!=seqs[edges[b,2],i]) muts[i,b]<-T

#Direction of mutations
mutsDir<-matrix(NA,l,nrow(edges))
for (i in 1:l) for (b in 1:nrow(edges)) {
  if (seqs[edges[b,1],i]>seqs[edges[b,2],i]) mutsDir[i,b]<-T
  if (seqs[edges[b,1],i]<seqs[edges[b,2],i]) mutsDir[i,b]<-F}

#Types of mutations
mutsTyp<-matrix(NA,l,nrow(edges))
for (i in 1:l) {
  if (length(unique(seqs[,i]))!=2) next # Site is not biallelic
  for (b in 1:nrow(edges)) {
    if (seqs[edges[b,1],i]==seqs[1,i] && seqs[edges[b,2],i]==seqs[1,i]) mutsTyp[i,b]<-1
    if (seqs[edges[b,1],i]==seqs[1,i] && seqs[edges[b,2],i]!=seqs[1,i]) mutsTyp[i,b]<-2
    if (seqs[edges[b,1],i]!=seqs[1,i] && seqs[edges[b,2],i]==seqs[1,i]) mutsTyp[i,b]<-3
    if (seqs[edges[b,1],i]!=seqs[1,i] && seqs[edges[b,2],i]!=seqs[1,i]) mutsTyp[i,b]<-4
  }
}
edgeLengths=tree$edge.length
edgeLengths[which(edgeLengths==0)]=1e-7

#Calculate scores
fout=file('score.csv','w')
s=0
tab=matrix(0,2,2)
npat=dim(muts)[1]
store=rep(1,200*200*200)
for (i in 1:npat) {
  w=which(rowSums(muts[i:npat,which(muts[i,]),drop=F])>=4)+i-1;
  for (j in w) {
    if (is.na(mutsTyp[i,1])||is.na(mutsTyp[j,1])) next#ignore if one of the two sites is not biallelic
    m=muts[i,]*2+muts[j,]
    tab[1,1]=sum(m==3)
    tab[1,2]=sum(m==1)
    tab[2,1]=sum(m==2)
    tab[2,2]=sum(m==0)
    if (tab[1,1]>199||tab[1,2]>199||tab[2,1]>199) {score=(round(log10(fisher.test(tab)$p.value)))} else {
      ind=tab[1,1]*40000+tab[1,2]*200+tab[2,1]
      if (store[ind]==1) store[ind]=round(log10(fisher.test(tab)$p.value))
      score=store[ind]}
    
    site1=mutsTyp[i,]
    site2=mutsTyp[j,]
    root=setdiff(tree$edge[,1],tree$edge[,2])
    bra=which(tree$edge[,1]==root)[1]
    if (site1[bra]==1||site1[bra]==2) root1=1 else root1=2
    if (site2[bra]==1||site2[bra]==2) root2=1 else root2=2

    likelihood <- function(a1,a2,b1,b2,eps) {
      Q0=t(matrix(c(0,a1,b1,0,a2,0,0,b1*eps,b2,0,0,a1*eps,0,b2/eps,a2/eps,0),4,4))
      Q=Q0-diag(rowSums(Q0),4,4)
      Q2=Q%*%Q
      Q3=Q2%*%Q
      Q=t(Q);Q2=t(Q2);Q3=t(Q3)
      map=matrix(c(1,2,5,6,3,4,7,8,9,10,13,14,11,12,15,16),4,4)
      combined=rep(NA,length(site1));for (i in 1:length(site1)) combined[i]=map[site1[i],site2[i]]
      Id=diag(c(1,1,1,1),4,4)
      den=a1*b1*eps*eps+a1*b2+a2*b1+a2*b2
      if (root1==1&&root2==1) num=a2*b2
      if (root1==1&&root2==2) num=a2*b1
      if (root1==2&&root2==1) num=a1*b2
      if (root1==2&&root2==2) num=a1*b1*eps*eps
      ret=log(num/den)
      ret=ret+sum(log(Id[combined]+Q[combined]*edgeLengths+Q2[combined]*edgeLengths^2/2+Q3[combined]*edgeLengths^3/6))
      if (is.na(ret)||is.infinite(ret)) ret=-1e10
      return(-ret)
    }
    
    likelihood0 <- function(a1,a2,b1,b2) {return(likelihood(a1,a2,b1,b2,1))}
    score2=0
    try({
    fit0<-suppressWarnings(mle(likelihood0,start = list(a1=1,a2=1,b1=1,b2=1),method = "L-BFGS-B", lower = rep(0,4)))
    fit <-suppressWarnings(mle(likelihood,start=list(a1=fit0@coef[[1]],a2=fit0@coef[[2]],b1=fit0@coef[[3]],b2=fit0@coef[[4]],eps=1),method = "L-BFGS-B", lower = rep(0,5)))
    lrt=2*(-fit@min+fit0@min)
    score2=pchisq(lrt,df=1,lower.tail = F)
    score2=round(log10(score2))},silent=T)
      if (score<= -10||score2<=-10) {
      s=s+1;writeLines(sprintf('%d,%d,%d,%d,%d,%d,%d,%d',i,j,tab[1,1],tab[1,2],tab[2,1],tab[2,2],score,score2),fout)
      }
  }}
close(fout)
save.image(file='epistasis.RData')
