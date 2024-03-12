#hmmm...precomputing integers and integer vectors didn't save any time

#got fftw3 link working, but required makevar file edits...see current one in documents
Rcpp::sourceCpp(code=
  "
  #define ARMA_USE_FFTW3
  #include <RcppArmadillo.h>
  #include <fftw3.h>
  // [[Rcpp::depends(RcppArmadillo)]]
  // [[Rcpp::export]]
  double lik_fun(arma::cx_cube R, arma::cx_cube X,
                 arma::uvec des, arma::uvec des_pos, arma::uvec des_n,
                 arma::uvec prune, arma::vec T,
                 double S,
                 arma::uword l, arma::uword m, arma::uword hm, arma::uword n,
                 arma::uvec ll, arma::uvec rr, arma::uvec np) {

    // arma::uword l = R.n_rows;
    // arma::uword m = R.n_slices;
    // arma::uword hm = m / 2;
    // arma::uword n = prune.n_elem;
    // arma::uvec ll = join_cols(arma::uvec { m - 1 }, arma::regspace<arma::uvec>(0, m - 2));
    // arma::uvec rr = join_cols(arma::regspace<arma::uvec>(1, m - 1), arma::uvec { 0 });
    // nonpad.inds = c((3*res/2+2):(2*res),1:(res/2+1))
    // arma::uvec np = join_cols(arma::regspace<arma::uvec>(3 * hm / 2 + 1, m - 1), arma::regspace<arma::uvec>(0, hm / 2));

    double C;
    double L;

    arma::cx_mat eigvals(l, m);
    arma::cx_vec eigval(l);

    arma::cx_cube eigvecs(l, l, m);
    arma::cx_cube inveigvecs(l, l, m);
    arma::cx_mat eigvec(l, l);

    arma::cx_cube XX = X;

    arma::cx_mat tmp_X(l, m);
    arma::cx_mat mholder(l, m);
    arma::cx_vec vholder(m);
    arma::vec rvholder(m);
    arma::vec vholder0(m);

    // eigendecompose R array
    for(arma::uword i=0; i<m; i++) {
      eig_gen(eigval, eigvec, R.slice(i));
      eigvals.col(i) = eigval;
      eigvecs.slice(i) = eigvec;
      inveigvecs.slice(i) = arma::inv(eigvec);
    }

    // loop over edges
    // something going wrong with indexing...unsure what
    for(arma::uword i=0; i<n; i++) {
      arma::uword e = prune(i);
      arma::mat tmp_x(l, hm, arma::fill::ones);

      // loop over descendants
      for(arma::uword j=0; j<des_n(e); j++){
        arma::uword d = des(des_pos(e) + j);

        tmp_X = XX.row_as_mat(d).st();
        for(arma::uword k=0; k<l; k++) {
          mholder.row(k) = arma::sum(inveigvecs.row_as_mat(k).st() % tmp_X);
        }
        mholder %= arma::exp(T(d) * eigvals);
        // might get messed up because vector/row vector confusion
        for(arma::uword k=0; k<l; k++) {
          vholder = arma::sum(eigvecs.row_as_mat(k).st() % mholder).st();
          C = real(vholder(0));
          rvholder = real(ifft(vholder, m));
          // I think this should work for cleaning...but hard to say
          rvholder.elem(find((rvholder.elem(ll) < 0 && rvholder.elem(rr) < 0) || rvholder < 1e-16)).fill(1e-16);
          rvholder = rvholder * C / arma::sum(rvholder);
          // multiply k-th row of tmp_x by depadded vholder
          tmp_x.row(k) %= rvholder.elem(np).t();
        }
      }

      // here, you have to rescale tmp_x to have max 1 to prevent overflow
      // and keep track of scalar!
      C = tmp_x.max();
      tmp_x = tmp_x / C;
      S += log(C);

      // check if i is n - 1 and do root calculation if so
      if(i == n - 1) {
        L = log(accu(square(tmp_x))) - log(accu(tmp_x)) + S;
      } else {
        // otherwise do this...
        //repad and fft each row/k of tmp_x and insert into XX(e, k, )
        // will probably have to add error-checking capabilities too...
        // but this will be good for first pass
        for(arma::uword k=0; k<l; k++) {
          // not sure if this works...
          // doesn't seem too
          // XX.tube(e, k).elem(np) = fft(tmp_x.row(k));
          // does this work?
          vholder0.elem(np) = tmp_x.row(k);
          XX.tube(e, k) = fft(vholder0, m);
        }
      }
    }

    // return root calculation
    return L;
  }
  ")

#still LOTS of errors, here's what you got for now:
# > you can't just combine regspace with integers using {}
# >> use join_cols instead!
# Rcpp::cppFunction(
#   "
#   List test_fun(arma::uword m){
#   arma::uvec ll = join_cols(arma::uvec { m - 1 }, arma::regspace<arma::uvec>(0, m - 2));
#   arma::uvec rr = join_cols(arma::regspace<arma::uvec>(1, m - 1), arma::uvec { 0 });
#   arma::uword hm = m / 2;
#
#   return List::create(ll, rr, np);
#   }
# ",
#   depends="RcppArmadillo")
# test_fun(100)
# > you don't understand how find and boolean functions interact at the moment...
# >> turned out to just be because you used complex values in find by accident!
# > will need to break down this function into parts I think to fully understand what's happening
# > but this a good draft/template!

#remaining issues:
# > need to declare fill (how?)
# > need to figure out how to insert into specific elements of XX tube...

#all errors fixed and it now compiles!!!
#now I'll just have to debug the actual function as it runs...ooh boy :S
#haven't tested yet, but expecting many, many errors

Rcpp::cppFunction(
  "
  List test_fun(arma::uvec des, arma::uvec des_pos, arma::uvec des_n, arma::uvec prune){


    for(arma::uword i=0; i<n; i++) {
      arma::uword e = prune(i);
      for(arma::uword j=0; j<des_n(e); j++){
        arma::uword d = des(des_pos(e) + j);
      }
    }
  }
",
  depends="RcppArmadillo")

Rcpp::cppFunction(
  "
  double lik_fun(arma::cx_cube R, arma::cx_cube X,
                 arma::uvec des, arma::uvec des_pos, arma::uvec des_n,
                 arma::uvec prune, arma::vec T,
                 double S,
                 arma::uword l, arma::uword m, arma::uword hm, arma::uword n,
                 arma::uvec ll, arma::uvec rr, arma::uvec np) {

    // arma::uword l = R.n_rows;
    // arma::uword m = R.n_slices;
    // arma::uword hm = m / 2;
    // arma::uword n = prune.n_elem;
    // arma::uvec ll = join_cols(arma::uvec { m - 1 }, arma::regspace<arma::uvec>(0, m - 2));
    // arma::uvec rr = join_cols(arma::regspace<arma::uvec>(1, m - 1), arma::uvec { 0 });
    // nonpad.inds = c((3*res/2+2):(2*res),1:(res/2+1))
    // arma::uvec np = join_cols(arma::regspace<arma::uvec>(3 * hm / 2 + 1, m - 1), arma::regspace<arma::uvec>(0, hm / 2));

    double C;
    double L;

    arma::cx_mat eigvals(l, m);
    arma::cx_vec eigval(l);

    arma::cx_cube eigvecs(l, l, m);
    arma::cx_cube inveigvecs(l, l, m);
    arma::cx_mat eigvec(l, l);

    arma::cx_cube XX = X;

    arma::cx_mat tmp_X(l, m);
    arma::cx_mat mholder(l, m);
    arma::cx_vec vholder(m);
    arma::vec rvholder(m);
    arma::vec vholder0(m);

    // eigendecompose R array
    for(arma::uword i=0; i<m; i++) {
      eig_gen(eigval, eigvec, R.slice(i));
      eigvals.col(i) = eigval;
      eigvecs.slice(i) = eigvec;
      inveigvecs.slice(i) = arma::inv(eigvec);
    }

    // loop over edges
    // something going wrong with indexing...unsure what
    for(arma::uword i=0; i<n; i++) {
      arma::uword e = prune(i);
      arma::mat tmp_x(l, hm, arma::fill::ones);

      // loop over descendants
      for(arma::uword j=0; j<des_n(e); j++){
        arma::uword d = des(des_pos(e) + j);

        tmp_X = XX.row_as_mat(d).st();
        for(arma::uword k=0; k<l; k++) {
          mholder.row(k) = arma::sum(inveigvecs.row_as_mat(k).st() % tmp_X);
        }
        mholder %= arma::exp(T(d) * eigvals);
        // might get messed up because vector/row vector confusion
        for(arma::uword k=0; k<l; k++) {
          vholder = arma::sum(eigvecs.row_as_mat(k).st() % mholder).st();
          C = real(vholder(0));
          rvholder = real(ifft(vholder, m));
          // I think this should work for cleaning...but hard to say
          rvholder.elem(find((rvholder.elem(ll) < 0 && rvholder.elem(rr) < 0) || rvholder < 1e-16)).fill(1e-16);
          rvholder = rvholder * C / arma::sum(rvholder);
          // multiply k-th row of tmp_x by depadded vholder
          tmp_x.row(k) %= rvholder.elem(np).t();
        }
      }

      // here, you have to rescale tmp_x to have max 1 to prevent overflow
      // and keep track of scalar!
      C = tmp_x.max();
      tmp_x = tmp_x / C;
      S += log(C);

      // check if i is n - 1 and do root calculation if so
      if(i == n - 1) {
        L = log(accu(square(tmp_x))) - log(accu(tmp_x)) + S;
      } else {
        // otherwise do this...
        //repad and fft each row/k of tmp_x and insert into XX(e, k, )
        // will probably have to add error-checking capabilities too...
        // but this will be good for first pass
        for(arma::uword k=0; k<l; k++) {
          // not sure if this works...
          // doesn't seem too
          // XX.tube(e, k).elem(np) = fft(tmp_x.row(k));
          // does this work?
          vholder0.elem(np) = tmp_x.row(k);
          XX.tube(e, k) = fft(vholder0, m);
        }
      }
    }

    // return root calculation
    return L;
  }
  ",
  depends="RcppArmadillo")

rm(old.des)

conv.dists<-get.DFTs(res,dx,c("JN","VG","NIG","BM"),
                     char.exp=TRUE)
testy<-array(1,c(k,k,2*res))
testy[dg]<- -1+do.call(rbind,lapply(seq_len(4),
                                    function(ii)
                                      conv.dists[[c(1,1,1,1)[ii]]](c(0,3,0,3)[ii],c(0,3,0,9)[ii],0)+
                                      conv.dists[[4]](c(1,0,4,0)[ii],0)))
R<-testy
if(exists("old.des")) des<-old.des
old.des<-des
des_n<-lengths(des)
des_pos<-c(0,cumsum(des_n)[-length(des)])
des<-unlist(des,use.names=FALSE)
prune<-prune.seq
des<-des-1
prune<-prune-1
lik_fun(R, X, des, des_pos, des_n, prune, elen, scalar.init) #seems to be taking a while...
microbenchmark::microbenchmark(lik_fun(R, X, des, des_pos, des_n, prune, elen, scalar.init))
#20ish millisecs for 30 tips with 2 states (plugged in)
#50ish millisecs for 50 tips with 2 states
#130ish millisecs for 100 tips with 3 states
#260ish millisecs for 200 tips with 3 states
#360-390ish millisec for 200 tips with 4 states
#pretty good I think!
#all for res of 1024

l<-k
m<-2*res
hm<-res
n<-length(prune)
ll<-c(m-1,0:(m-2))
rr<-c(1:(m-1),0)
np<-c((3*hm/2+1):(m-1),0:(hm/2))
lik_fun(R, X, des, des_pos, des_n, prune, elen, scalar.init,
        l,m,hm,n,ll,rr,np) #seems to be taking a while...
microbenchmark::microbenchmark(lik_fun(R, X, des, des_pos, des_n, prune, elen, scalar.init,
                                       l,m,hm,n,ll,rr,np))

plot(test.test[1,])
sum(test.test[2,]*2048)
#oops! Made some loops go forever by accident
#yep, size errors...gotta dig in some more
#wish it gave me more informative error messages
#hmmm...getting there; now getting some out of bound issues
#god zero-indexing is conufsing...
#But I think that's the only thing frustrating you at the moment
#Might be best to prep uvecs for 0 indexing?

#Yeah, just need to figure out indexing code
#Holy shit I might have figured it out...
#Maybe? Seems ridiculously high
#so code runs, but there's clearly something wrong with likelihood calcs...

#Just something weird going on with overflow...
#seems to be 1 contrast where things go to hell...
which.max(Re(test.test[,2,1]))
old.des[[56]]
old.des[[54]]
plot(tree)
tmp.tree<-reorder(tree,"pruningwise")
plot(tmp.tree)
edgelabels()
#oh, might just be messing up the indexing...
old.des[[56]]
#no, that's not it...something's going to hell with the fft part I think...
matplot(apply(test.test[47,,],1,IDFT),type="l")
matplot(apply(test.test[13,,],1,IDFT),type="l")
matplot(apply(test.test[14,,],1,IDFT),type="l")

#figured it out! Just zero-padded the vector after ffting rather than before


for(i in seq_along(prune)-1){
  e<-prune[i+1]
  for(j in seq_len(des_n[e+1])-1){
    cat(e,":",des[des_pos[e+1]+j+1],"\n")
  }
}

Rcpp::cppFunction("
List vec_eig(arma::cx_cube x) {

  arma::cx_mat eigvals(x.n_rows, x.n_slices);
  arma::cx_vec eigval(x.n_rows);
  arma::cx_cube eigvecs(x.n_rows, x.n_cols, x.n_slices);
  arma::cx_cube inveigvecs = eigvecs;
  arma::cx_mat eigvec(x.n_rows, x.n_cols);
  // arma::uvec indices;

  for(int i=0; i<x.n_slices; i++) {
    eig_gen(eigval, eigvec, x.slice(i));
    // indices = arma::sort_index(eigval);
    eigvals.col(i) = eigval;
    eigvecs.slice(i) = eigvec;
    inveigvecs.slice(i) = arma::inv(eigvec);
  }

  return List::create(eigvals, eigvecs, inveigvecs);
}
",
depends="RcppArmadillo")


#figured it out!!!
#it seems like, for whatever reason, the complex and imaginary components get "swapped" by Rcpp!
#or something like that...
Rcpp::cppFunction("
arma::cx_mat vec_exp_mult(double t, arma::cx_mat x, arma::cx_mat eigvals, arma::cx_cube eigvecs, arma::cx_cube inveigvecs) {

  arma::cx_mat out1(eigvals.n_rows, eigvals.n_cols);
  arma::cx_mat tmp_exp(eigvals.n_rows, eigvals.n_cols);
  arma::cx_mat out2(eigvals.n_rows, eigvals.n_cols);

  // something seems to be going wrong below...
  // even if you sort the eigenvalues and such beforehand...
  // I figured it out! .t() does the conjugate transpose for complex matrices
  // You want .st() instead!
  for(int i=0; i<eigvals.n_rows; i++) {
    out1.row(i) = arma::sum(inveigvecs.row_as_mat(i).st() % x);
  }
  out1 %= arma::exp(t * eigvals);
  for(int i=0; i<eigvals.n_rows; i++) {
    out2.row(i) = arma::sum(eigvecs.row_as_mat(i).st() % out1);
  }

  return out2;
}
",
depends="RcppArmadillo")
