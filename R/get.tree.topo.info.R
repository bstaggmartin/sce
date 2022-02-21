#edge nedge+1 will always be the root age, even when pruning is FALSE
.get.tree.topo.info<-function(tree,pruning=FALSE){
  #reorder tree
  attr(tree,'order')<-NULL
  if(pruning){
    tree<-reorder(tree,order='pruningwise')
  }else{
    tree<-reorder(tree,order='cladewise')
  }
  #get edge lengths
  hgt<-max(node.depth.edgelength(tree))
  elen<-tree[['edge.length']]
  #get ancestors
  nodes<-tree[['edge']]
  anc<-match(nodes[,1,drop=FALSE],nodes[,2,drop=FALSE])
  nedge<-length(anc)+1
  anc[is.na(anc)]<-nedge #root edge index is number of edges + 1
  #get descendants, note where tips are, make pruning sequence
  sis<-des<-rep(list(integer(0)),nedge)
  topo.seq<-seq_len(nedge)
  tmp<-split(topo.seq[-nedge],anc)
  des[as.numeric(names(tmp))]<-tmp
  tips<-which(lengths(des)==0)
  #reorder tips
  tips.ord<-tree[['tip.label']][nodes[tips,2,drop=FALSE]]
  #finalize nodes matrix
  nodes<-rbind(nodes,c(NA,length(tips)+1))
  #get edge sequence
  if(pruning){
    topo.seq<-topo.seq[-tips]
  }else{
    topo.seq<-topo.seq[-nedge]
  }
  #sister edges...more complicated, but needed for ancestral state recon
  ndes<-lengths(des)
  di<-ndes==2
  di.des<-unlist(des[di])
  if(length(di.des)){
    odds<-seq.int(1,length(di.des),2)
    evens<-odds+1
    odds<-di.des[odds]
    evens<-di.des[evens]
    sis[odds]<-evens
    sis[evens]<-odds
  }
  poly.des<-des[!di]
  if(length(poly.des)){
    ndes<-ndes[!di]
    unlist.poly.des<-unlist(poly.des,use.names=FALSE)
    sis[unlist.poly.des]<-rep(poly.des,ndes)
    foo<-function(x){
      tmp<-sis[[x]]
      tmp[tmp!=x]
    }
    sis[unlist.poly.des]<-lapply(unlist.poly.des,foo)
  }
  out<-list('seq'=topo.seq,
            'nedge'=nedge,
            'hgt'=hgt,
            'elen'=elen,
            'anc'=anc,
            'des'=des,
            'sis'=sis,
            'tips'=tips,
            'tips.ord'=tips.ord,
            'nodes'=nodes)
  if(pruning){
    names(out)[1]<-'prune.seq'
  }else{
    names(out)[1]<-'clade.seq'
  }
  out
}
