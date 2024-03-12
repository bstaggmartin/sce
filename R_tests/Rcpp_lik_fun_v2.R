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
          // cleaning steps add barely any time --> probs best to keep in
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
        //outherwise repad and fft each row/k of tmp_x and insert into XX(e, k, )
        for(arma::uword k=0; k<l; k++) {
          vholder0.elem(np) = tmp_x.row(k);
          XX.tube(e, k) = fft(vholder0, m);
        }
      }
    }

    // return root calculation
    return L;
  }
  ")
