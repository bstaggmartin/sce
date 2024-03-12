#define ARMA_USE_FFTW3
#include "RcppArmadillo.h"
#include "fftw3.h"

// [[Rcpp::depends(RcppArmadillo)]]

// likelihood function
// [[Rcpp::export]]
double sce_lik(arma::cx_cube R, arma::cx_cube X,
               arma::uvec des, arma::uvec des_pos, arma::uvec des_n,
               arma::uvec prune, arma::vec T,
               double S) {

  arma::uword l = R.n_rows;
  arma::uword m = R.n_slices;
  arma::uword hm = m / 2;
  arma::uword n = prune.n_elem;
  arma::uvec ll = join_cols(arma::uvec { m - 1 },
                            arma::regspace<arma::uvec>(0, m - 2));
  arma::uvec rr = join_cols(arma::regspace<arma::uvec>(1, m - 1),
                            arma::uvec { 0 });
  arma::uvec np = join_cols(arma::regspace<arma::uvec>(3 * hm / 2 + 1, m - 1),
                            arma::regspace<arma::uvec>(0, hm / 2));

  double C;

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
      for(arma::uword k=0; k<l; k++) {
        vholder = arma::sum(eigvecs.row_as_mat(k).st() % mholder).st();
        rvholder = real(ifft(vholder, m));
        // cleaning steps add barely any time --> probs best to keep in
        rvholder.elem(find((rvholder.elem(ll) < 0 && rvholder.elem(rr) < 0) || rvholder < 1e-16)).fill(1e-16);
        rvholder = rvholder * real(vholder(0)) / arma::sum(rvholder);
        // multiply k-th row of tmp_x by depadded vholder
        tmp_x.row(k) %= rvholder.elem(np).t();
      }
    }

    // here, you have to rescale tmp_x to have max 1 to prevent overflow
    // and keep track of scalar!
    C = tmp_x.max();
    tmp_x /= C;
    S += log(C);

    // check if i is n - 1 and do root calculation if so
    if(i == n - 1) {
      S += log(accu(square(tmp_x))) - log(accu(tmp_x));
    } else {
      //outherwise repad and fft each row/k of tmp_x and insert into XX(e, k, )
      for(arma::uword k=0; k<l; k++) {
        vholder0.elem(np) = tmp_x.row(k);
        XX.tube(e, k) = fft(vholder0, m);
      }
    }
  }

  // return root calculation
  return S;
}

// ancestral state reconstruction function
// [[Rcpp::export]]
arma::cube sce_rec(arma::cx_cube R, arma::cx_cube bw_R, arma::cx_cube X,
                   arma::uvec des, arma::uvec des_pos, arma::uvec des_n,
                   arma::uvec prune, arma::vec T,
                   arma::uvec tip) {

  arma::uword l = R.n_rows;
  arma::uword m = R.n_slices;
  arma::uword hm = m / 2;
  arma::uword n = prune.n_elem;
  arma::uvec ll = join_cols(arma::uvec { m - 1 },
                            arma::regspace<arma::uvec>(0, m - 2));
  arma::uvec rr = join_cols(arma::regspace<arma::uvec>(1, m - 1),
                            arma::uvec { 0 });
  arma::uvec np = join_cols(arma::regspace<arma::uvec>(3 * hm / 2 + 1, m - 1),
                            arma::regspace<arma::uvec>(0, hm / 2));

  arma::cx_mat eigvals(l, m);
  arma::cx_mat bw_eigvals(l, m);
  arma::cx_vec eigval(l);

  arma::cx_cube eigvecs(l, l, m);
  arma::cx_cube bw_eigvecs(l, l, m);
  arma::cx_cube inveigvecs(l, l, m);
  arma::cx_cube bw_inveigvecs(l, l, m);
  arma::cx_mat eigvec(l, l);

  arma::uword ne = X.n_rows;
  arma::uword nt = tip.n_elem;
  // store ffts at terminus of each branch based on its descendants ONLY
  arma::cx_cube XX = X;
  // stores iffts at terminus of each branch based on its descendants ONLY
  // eventually gets multiplied by zz and becomes output
  arma::cube d_xx(ne, l, hm, arma::fill::ones);
  // stores iffts at beginning of each branch based on its descendants ONLY
  // used in intermediate calculations of preorder traversal
  arma::cube a_xx(ne, l, hm);
  // stores iffts at terminus of each branch based on its non-descendants ONLY
  arma::cube zz(ne, l, hm);

  arma::cx_mat tmp_X(l, m);
  arma::cx_mat mholder(l, m);
  arma::cx_vec vholder(m);
  arma::vec rvholder(m);
  arma::vec vholder0(m);
  arma::mat tmp_xx(l, hm);

  // eigendecompose R array
  for(arma::uword i=0; i<m; i++) {
    eig_gen(eigval, eigvec, R.slice(i));
    eigvals.col(i) = eigval;
    eigvecs.slice(i) = eigvec;
    inveigvecs.slice(i) = arma::inv(eigvec);

    eig_gen(eigval, eigvec, bw_R.slice(i));
    bw_eigvals.col(i) = eigval;
    bw_eigvecs.slice(i) = eigvec;
    bw_inveigvecs.slice(i) = arma::inv(eigvec);
  }

  // need to ifft tip dists...
  for(arma::uword i=0; i<nt; i++) {
    arma::uword e = tip(i);

    for(arma::uword k=0; k<l; k++) {
      rvholder = real(ifft(vectorise(XX.tube(e, k)), m));
      rvholder.elem(find((rvholder.elem(ll) < 0 && rvholder.elem(rr) < 0) || rvholder < 1e-16)).fill(1e-16);
      rvholder = rvholder * real(XX(e, k, 0)) / arma::sum(rvholder);
      d_xx.tube(e, k) = rvholder.elem(np);
    }
  }

  // loop over edges
  for(arma::uword i=0; i<n; i++) {
    arma::uword e = prune(i);

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
        rvholder = real(ifft(vholder, m));
        // cleaning steps add barely any time --> probs best to keep in
        rvholder.elem(find((rvholder.elem(ll) < 0 && rvholder.elem(rr) < 0) || rvholder < 1e-16)).fill(1e-16);
        rvholder = rvholder * real(vholder(0)) / arma::sum(rvholder);
        // store in a_xx for d
        a_xx.tube(d, k) = rvholder.elem(np);
        // multiply d_xx by d's a_xx for e
        d_xx.tube(e, k) %= a_xx.tube(d, k);
      }
    }

    // here, you have to rescale d_xx to have max 1 to prevent overflow
    d_xx.row(e) /= d_xx.row(e).max();

    // not necessary to treat root specially anymore, I think...
    // repad and fft each row/k of tmp_x and insert into XX(e, k, )
    for(arma::uword k=0; k<l; k++) {
      vholder0.elem(np) = vectorise(d_xx.tube(e, k));
      XX.tube(e, k) = fft(vholder0, m);
    }
  }

  // need to initialize root...
  zz.row(prune(n - 1)).ones();
  // and standardize d_xx for root...
  d_xx.row(prune(n - 1)) /= accu(d_xx.row(prune(n - 1)));

  // loop over edges in reverse order
  for(arma::uword i=0; i<n; i++) {
    arma::uword e = prune(n - 1 - i);

    // loop over descendants
    for(arma::uword j=0; j<des_n(e); j++){
      arma::uword d = des(des_pos(e) + j);

      tmp_xx = zz.row_as_mat(e).t();
      // loop over sisters by skipping current j...
      for(arma::uword jj=0; jj<des_n(e); jj++){
        if(j != jj){
          arma::uword s = des(des_pos(e) + jj);

          tmp_xx %= a_xx.row_as_mat(s).t();
        }
      }
      // fft tmp_xx, store as tmp_X
      for(arma::uword k=0; k<l; k++) {
        vholder0.elem(np) = tmp_xx.row(k);
        tmp_X.row(k) = fft(vholder0, m).st();
      }
      // now use bw_R to solve forward equations!
      // technically, I think d_xx could become the output if you just multiply at each step...?
      // yeah, because you never use d_xx here...
      // that's probably the way to go, but I need to think on it some more
      for(arma::uword k=0; k<l; k++) {
        mholder.row(k) = arma::sum(bw_inveigvecs.row_as_mat(k).st() % tmp_X);
      }
      mholder %= arma::exp(T(d) * bw_eigvals);
      for(arma::uword k=0; k<l; k++) {
        vholder = arma::sum(bw_eigvecs.row_as_mat(k).st() % mholder).st();
        rvholder = real(ifft(vholder, m));
        // cleaning steps add barely any time --> probs best to keep in
        rvholder.elem(find((rvholder.elem(ll) < 0 && rvholder.elem(rr) < 0) || rvholder < 1e-16)).fill(1e-16);
        rvholder = rvholder * real(vholder(0)) / arma::sum(rvholder);
        // store in zz
        zz.tube(d, k) = rvholder.elem(np);
        // multiply d_xx by zz for d
        d_xx.tube(d, k) %= zz.tube(d, k);
      }

      // now you just have to rescale zz to prevent underflow
      zz.row(d) /= zz.row(d).max();
      // and d_xx to sum to 1!
      d_xx.row(d) /= accu(d_xx.row(d));
      // I think that's somehow it...

    }
  }

  // return cube of dists
  return d_xx;
}
