#could use better error handling in the future...
#note that the correspondences between parameters and the concept of freq, jump, and skew...
#are rougher in the case of "infinitely active" Levy process like NIG and VG, but I tried my...
#best to make it intuitive
#streamlined into single function that outputs list of DFT functions for specified distributions
#can also output list of characteristic exponent functions
#' @export

get.DFTs<-function(nx,dx,
                  dist=c("BM","JN","VG","NIG","DI","NO"),
                  char.exp=FALSE,
                  x0=0){
  if(is.null(dist)) dist<-1
  if(char.exp){
    dist.choices<-c("BM","JN","VG","NIG")
  }else{
    dist.choices<-c("BM","JN","VG","NIG","DI","NO")
  }
  if(is.character(dist)){
    dist<-pmatch(dist,dist.choices)
  }else if(is.numeric(dist)){
    dist[dist>length(dist.choices)]<-NA
  }else{
    stop("Didn't recognize distribution choices")
  }
  #should probs have a warning about this, but just defaults to BM
  dist[is.na(dist)]<-1

  tmp<-seq(0,-pi/dx,length.out=nx+1)
  bases<-1i*c(tmp,-tmp[nx:2])
  tmp<-tmp^2
  bases<-list(bases,c(tmp,tmp[nx:2]))

  if(char.exp){
    lapply(dist,
           function(ii) switch(ii,
                               function(rate,drift){
                                 -rate*bases[[2]]/2+drift*bases[[1]]
                               },
                               function(freq,jump,skew){
                                 freq*(exp(-jump*bases[[2]]/2+skew*bases[[1]])-1)
                               },
                               function(freq,jump,skew){
                                 if(freq==0){
                                   rep(0,2*nx)
                                 }else{
                                   -freq*log(1-skew*bases[[1]]/freq+jump*bases[[2]]/(2*freq))
                                 }
                               },
                               function(freq,jump,skew){

                                 # if(freq==0){
                                 #   -t*jump*bases[[2]]
                                 # }else{
                                   mskew2<-freq^2*(1-skew^2)
                                   freq*jump*(sqrt(mskew2)-sqrt(mskew2-2*freq*skew*bases[[1]]+bases[[2]]))
                                 # }

                                 # t*jump*(sqrt(freq^2-skew^2)-sqrt(freq^2-skew^2-2*skew*bases[[1]]+bases[[2]]))


                                 # if(freq==0){
                                 #   rep(0,2*nx)
                                 # }else{
                                 #   inv.freq<-1/freq
                                 #   mskew2<-1-skew^2
                                 #   jump*(inv.freq*sqrt(mskew2)-inv.freq*sqrt(mskew2-2*skew*freq*bases[[1]]+freq^2*bases[[2]]))
                                 # }
                               }))
  }else{
    lapply(dist,
           function(ii) switch(ii,
                               function(t,rate,drift){
                                 exp(t*(-rate*bases[[2]]/2+drift*bases[[1]]))
                               },
                               function(t,freq,jump,skew){
                                 exp(t*freq*(exp(-jump*bases[[2]]/2+skew*bases[[1]])-1))
                               },
                               function(t,freq,jump,skew){
                                 if(freq==0){
                                   rep(1,2*nx)
                                 }else{
                                   (1-skew*bases[[1]]/freq+jump*bases[[2]]/(2*freq))^(-t*freq)
                                 }
                               },
                               function(t,freq,jump,skew){

                                 #multiply jump and skew by freq I think...
                                 #works! final parameterization
                                 #decided to not do the if case for consistency of parameterizations
                                 #now doesn't converge to cauchy as freq goes to 0, but appropriately goes to no change
                                 # if(freq==0){
                                 #   exp(-t*jump*bases[[2]])
                                 # }else{
                                   mskew2<-freq^2*(1-skew^2)
                                   exp(t*freq*jump*(sqrt(mskew2)-sqrt(mskew2-2*freq*skew*bases[[1]]+bases[[2]])))
                                 # }

                                 # exp(t*jump*(sqrt(freq^2-skew^2)-sqrt(freq^2-skew^2-2*skew*bases[[1]]+bases[[2]])))

                                 # exp(t*jump*(sqrt(freq^2-skew^2)-sqrt(freq^2-(skew+bases[[1]])^2)))
                                 # if(freq==0){
                                 #   rep(1,2*nx)
                                 # }else{
                                 #   # mskew2<-1-skew^2
                                 #   # exp(t*jump*freq*(sqrt(mskew2)-sqrt(mskew2-2*skew/freq*bases[[1]]+bases[[2]]/freq^2)))
                                 # }
                               },
                               function(loc){
                                 exp((loc-x0)*bases[[1]])
                               },
                               function(loc,scale){
                                 exp((loc-x0)*bases[[1]]-scale*bases[[2]]/2)
                               }))
  }
}

#' @export
unfixed.NO.DFT<-function(nx,dx,x0){
  tmp<-seq(0,-pi/dx,length.out=nx+1)
  bases<-1i*c(tmp,-tmp[nx:2])
  tmp<-tmp^2
  bases<-list(bases,-c(tmp,tmp[nx:2])/2)
  list(base=function(loc) (loc-x0)*bases[[1]],
       modder=function(scale) scale*bases[[2]])
}
