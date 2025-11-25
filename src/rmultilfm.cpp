// This script implement multi-study robust linear factor model using variational inference
// Date: 2023-09-17
// Here, we exrt the orthogonal conditions of A'A, Bs'Bs after the algorithm converges!
// Revised log:


#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]
#include<ctime>
//#include<boost/math/tools/minima.hpp>

#define INT_MIN1 (-INT_MAX - 1)

using namespace Rcpp;
using namespace arma;
using namespace std;
// using boost::math::tools::brent_find_minima;


// Define global variables
//double X_count_ij, a_i, invLambda_j,Mu_x_ij;
/*
 * Auxiliary
 */
//' @keywords internal
//' @noRd
//'

// diag(W0* Cki * W0.t())
vec decomp(const mat& Cki, const mat& W0){
  vec s, tmp1;
  mat U, V, WC12;
  svd(U, s, V, Cki);
  WC12 = W0 * (U * diagmat(sqrt(s))); // p * q
  tmp1 = sum(WC12 % WC12, 1);
  return tmp1;
}

//
// List irlbaCpp1(const mat& X, const int& q){
//   Rcpp::Environment irlba("package:irlba");
//   Rcpp::Function f = irlba["irlba"];
//
//   return f(Named("A") = X, Named("nv") = q);
//
// }





void update_mu(const field<mat>& Xf, const mat& A, const field<mat>& Bf,
  const field<mat>& Muf_f, const field<mat>& Muf_h, const field<vec>& Phif, mat& bmu){

  int s, S = Xf.n_elem, p= Xf(0).n_cols;
  rowvec tmp_vec;
  for(s=0; s<S; ++s){
    mat dX = (Xf(s) - Muf_f(s)*A.t() - Muf_h(s)*Bf(s).t()) / repmat(Phif(s), 1, p);
    tmp_vec = sum(dX); // 1*p
    bmu.col(s) = tmp_vec.t() /accu(1.0/Phif(s));
  }
}

mat Sphi_fun(const cube& Sfs, const vec& Phifs){ // Sfs \in n_s * q
  int i, n_s = Phifs.n_elem, q = Sfs.n_cols;
  mat Ss(q, q, fill::zeros);
  for(i=0; i<n_s; ++i){
    Ss +=  Sfs.slice(i) / Phifs(i);
  }
  return Ss;
}
// check this function carefully
// void update_A(const field<mat>& Xf, const mat& bmu, const mat& LambdaMat,
//              const field<mat>& Bf, const field<mat>& Muf_f,
//               const field<mat>& Muf_h,const field<cube>& Sf_f, const field<vec>& Phif,
//               mat& A){
//   int s, S = Xf.n_elem, n_s, p = Xf(0).n_cols, q=Muf_f(0).n_cols;
//   mat A1(p*q, p*q, fill::zeros), A2(p, q, fill::zeros), A_tmp;
//
//   for(s=0; s<S; ++s){
//     n_s = Xf(s).n_rows;
//
//     A_tmp = Muf_f(s).t() * (Muf_f(s) / repmat(Phif(s), 1, q)) + Sphi_fun(Sf_f(s), Phif(s)); // q*q
//     // Rprintf("check point 1 \n");
//     A1 += kron(A_tmp, diagmat(1.0/LambdaMat.col(s)));
//     // Rprintf("check point 2 \n");
//     A2 += trans((Xf(s) - repmat(bmu.col(s).t(), n_s, 1) - Muf_h(s) * Bf(s).t())/repmat(LambdaMat.col(s).t(),n_s,1)) *
//       (Muf_f(s) / repmat(Phif(s), 1, q)); // p*q
//     // Rprintf("check point 3 \n");
//   }
//   mat vecA = A1.i()*A2.as_col();
//   // Rprintf("check point 4 \n");
//   vecA.reshape(p, q);
//   // Rprintf("check point 5 \n");
//   A =  vecA;
// }
// Fast version of updating A
void update_A(const field<mat>& Xf, const mat& bmu, const mat& LambdaMat,
              const field<mat>& Bf, const field<mat>& Muf_f,
              const field<mat>& Muf_h,const field<cube>& Sf_f, const field<vec>& Phif,
              mat& A){
  int j, s, S = Xf.n_elem, n_s, p = Xf(0).n_cols, q=Muf_f(0).n_cols;
  mat  A2(p, q, fill::zeros), A_tmp, A1; // A1(p*q, p*q, fill::zeros),
  for(s=0; s<S; ++s){
    n_s = Xf(s).n_rows;
    // Rprintf("check point 2 \n");
    A2 += trans((Xf(s) - repmat(bmu.col(s).t(), n_s, 1) - Muf_h(s) * Bf(s).t())/repmat(LambdaMat.col(s).t(),n_s,1)) *
      (Muf_f(s) / repmat(Phif(s), 1, q)); // p*q
  }
  for(j=0; j<p; ++j){ // linear computational complexity
    A1 = zeros(q, q);
    for(s=0; s<S; ++s){
      n_s = Xf(s).n_rows;
      A_tmp = Muf_f(s).t() * (Muf_f(s) / repmat(Phif(s), 1, q)) + Sphi_fun(Sf_f(s), Phif(s)); // q*q
      // Rprintf("check point 1 \n");
      A1 += A_tmp /LambdaMat(j,s);
    }
    A.row(j) = A2.row(j) * A1.i();
  }
  // Rprintf("check point 4 \n");

}

// update B_s
void update_B(const field<mat>& Xf, const mat& bmu, const mat& LambdaMat,
                          const mat& A, const field<mat>& Muf_f,
                          const field<mat>& Muf_h,const field<cube>& Sf_h, const field<vec>& Phif,
                          field<mat>& Bf){
  int s, n_s, q_s, S = Xf.n_elem, p = A.n_rows;

  for(s=0; s<S; ++s){
    q_s = Muf_h(s).n_cols;
    mat B1(q_s, q_s);
    mat B2(p, q_s);
    n_s = Xf(s).n_rows;

    B1 = Muf_h(s).t() * (Muf_h(s) / repmat(Phif(s), 1, q_s)) + Sphi_fun(Sf_h(s), Phif(s));
    B2 = trans(Xf(s) - repmat(bmu.col(s).t(), n_s, 1) - Muf_f(s) * A.t()) * (Muf_h(s) / repmat(Phif(s), 1, q_s));
    Bf(s) = B2 * B1.i();
  }
}


// update lambda_s
void update_Lambda(const field<mat>& Xf, const mat& bmu,
                         const mat& A, const field<mat>& Bf, const field<mat>& Muf_f, const field<cube>& Sf_f,
                         const field<mat>& Muf_h,const field<cube>& Sf_h, const field<vec>& Phif,
                         const double& nu,mat& LambdaMat){

  int s, n_s, S = Xf.n_elem,  p = A.n_rows;
  double c1 = (nu+p)/nu;
  vec Lambdavec, a_vec, b_vec;
  mat tmpSf, tmpSh, dX;
  for(s = 0; s< S; ++s){
    n_s = Xf(s).n_rows;
    tmpSf = Sphi_fun(Sf_f(s), Phif(s));
    tmpSh = Sphi_fun(Sf_h(s), Phif(s));
    a_vec = decomp(tmpSf, A);
    b_vec = decomp(tmpSh, Bf(s));
    mat dX= Xf(s) - repmat(bmu.col(s).t(), n_s, 1) - Muf_f(s) * A.t() - Muf_h(s)* Bf(s).t();
    Lambdavec  = trans(mean(dX % dX % repmat(1.0/Phif(s), 1, p)))  + a_vec/n_s + b_vec/n_s;
    LambdaMat.col(s) = Lambdavec *c1;
  }

}


void add_IC_Orth(mat& A, field<mat>& Bf){
  // Try the orthogonal matrix method
  int r, q = A.n_cols, qs1 = Bf(0).n_cols, r_max=Bf.n_elem;
  mat U, V;
  vec s;
  mat AB1 = join_rows(A,Bf(0));
  svd(U, s, V, AB1);
  vec signU = sign(U.row(0).t()); // can be speed up by only caculate the first columns
  A = U.cols(0, q-1) * diagmat(s.subvec(0,q-1) % signU.subvec(0, q-1));
  Bf(0) = U.cols(q, q+qs1-1) * diagmat(s.subvec(q,q+qs1-1) % signU.subvec(q,q+qs1-1));

  for(r=1; r<r_max; ++r){
    qs1 = Bf(r).n_cols;
    mat U1, V1;
    vec s1;
    svd(U1, s1, V1, Bf(r));
    vec signU1 = sign(U1.row(0).t());
    Bf(r) = U1.cols(0,qs1-1) * diagmat(s1 % signU1.subvec(0,qs1-1));
  }


}

double logCp(const double& nu, const int& p){

  double p1 = (p+nu)/2.0;
  // Rprintf("Gamma(p+nu/2)= %4f \n", lgamma(p1));
  // Rprintf("Gamma(nu/2)= %4f \n", lgamma(nu/2));
  double y = -p/2*log(nu) + lgamma(p1) - lgamma(nu/2); // log(Gamma(x)) in std library.
  return y;
}

double calELBO(const field<mat>& Xf, const mat& bmu,  const mat& LambdaMat,
                const mat& A, const field<mat>& Bf, const field<mat>& Muf_f,const field<cube>& Sf_f,
                const field<mat>& Muf_h,const field<cube>& Sf_h, const double& nu){

  int s,i, n_s, q_s, S = Xf.n_elem, p = Xf(0).n_cols, q= A.n_cols;
  double ELBO=0.0;
  double Xterm1 = 0.0, dimreduction_term2=0.0;
  mat AtLA, Ssf;
  // Rprintf("Step into loop in calELBO \n");
  for(s=0; s< S; ++s){
    n_s = Xf(s).n_rows;
    q_s = Bf(s).n_cols;
    // Rprintf("Step into loop in calELBO1 \n");
    mat dX = (Xf(s) - repmat(bmu.col(s).t(), n_s, 1) - Muf_f(s) * A.t() -
      Muf_h(s) * Bf(s).t()) / repmat(sqrt(LambdaMat.col(s).t()),n_s,1); // n_s * p
    // Rprintf("Step into loop in calELBO2 \n");
    vec tmpX = sum(dX% dX, 1); // n_s*1
    AtLA = A.t()*(A / repmat(LambdaMat.col(s),1, q));
    // Rprintf("Step into loop in calELBO3 \n");
    mat BtLB = Bf(s).t()*(Bf(s) / repmat(LambdaMat.col(s),1, q_s)); // elementwise division
    rowvec Trvec(n_s, fill::zeros);
    for(i=0; i<n_s; ++i){
      Trvec(i) = trace(AtLA* Sf_f(s).slice(i)) + trace(BtLB* Sf_h(s).slice(i));
    }
    // Rprintf("Step into loop in calELBO4 \n");
    Xterm1 += -0.5*(nu+p)* accu(log(tmpX.t()/nu+ Trvec/nu + 1.0));
    Xterm1 += n_s* logCp(nu, p) - 0.5*n_s *accu(log(LambdaMat.col(s)));
    // Rprintf("Step into loop in calELBO5 \n");
    dimreduction_term2 += -0.5*(accu(Muf_f(s)% Muf_f(s))+accu(Muf_h(s)% Muf_h(s)));
    for(i=0; i<n_s; ++i){
      dimreduction_term2 += -0.5* (trace(Sf_f(s).slice(i))+trace(Sf_h(s).slice(i)));
    }

  }


  double entropy=0.0, val1, val2, sign;
  for(s=0; s< S; ++s){
    n_s = Xf(s).n_rows;
    // Rprintf("Step into loop in calELBO6 \n");
    for(i=0; i<n_s; ++i){
      // entropy += 0.5 * (log_det_sympd(Sf_f(s).slice(i)) + log_det_sympd(Sf_h(s).slice(i)));
      log_det(val1, sign, Sf_f(s).slice(i));
      log_det(val2, sign, Sf_h(s).slice(i));
      entropy += 0.5 * (val1+val2);
    }
    // entropy += 0.5 * (accu(log(Sf_f(s)+(Sf_f(s)==0))) + accu(log(Sf_h(s)+(Sf_h(s)==0))));
  }


  ELBO = Xterm1 + dimreduction_term2 + entropy;
  return ELBO;
}



void VB_EstepVI(const field<mat>& Xf, const mat& bmu, const mat& LambdaMat, const mat& A,
                const field<mat>& Bf,
                field<mat>& Muf_f,field<cube>& Sf_f,field<mat>& Muf_h,field<cube>& Sf_h,
                const double& nu, const field<vec>& Phif){

  int s,i, n_s, q_s, S = Xf.n_elem, q= A.n_cols, p = Xf(0).n_cols;
  double c1 = (nu+p)/nu;
  //  ## VB E-step
  // update posterior variance of y: Sf_y

  double elbo1 = calELBO(Xf, bmu, LambdaMat, A, Bf, Muf_f, Sf_f, Muf_h, Sf_h, nu);
  // update posterior variance of f_i, Sf_f
  // update posterior mean of f_i, Muf_f

  mat Si,Si_inv, Xsi;

  // Rprintf("Update Mu_f and S_f \n");
  for(s=0; s<S; ++s){
    n_s = Xf(s).n_rows;
    Si = A.t() * (A / repmat(LambdaMat.col(s), 1, q));
    Xsi = c1* ((Xf(s) - repmat(bmu.col(s).t(), n_s, 1) - Muf_h(s) * Bf(s).t()) *
      (A / repmat(LambdaMat.col(s), 1, q))); // n_s * q
    // (Xsi.row(0)/p).print();
    for(i=0; i<n_s; ++i){
      Si_inv =  Si * c1 / Phif(s)(i) + eye(q,q); // diagmat(invLambda) * B
      Sf_f(s).slice(i) = Si_inv.i();
      Muf_f(s).row(i) = Xsi.row(i) * (Sf_f(s).slice(i) / Phif(s)(i)); // 1 *q
    }
    // Muf_f(s).rows(0,4).print();
    // Rprintf("The rows of Sf_f(s) is: %d\n", Sf_f(s).n_rows);
    // Rprintf("The rows of Phif(s) is: %d\n", Phif(s).n_elem);


    // vec s2, S2_inv;
    // n_s = Xf(s).n_rows;
    // q_s = Bf(s).n_cols;
    // s2 = diagvec(Bf(s).t() * (Bf(s) / repmat(LambdaMat.col(s), 1, q_s)));
    // for(i=0; i<n_s; ++i){
    //   S2_inv =  s2 * c1 / Phif(s)(i) + 1.0; // diagmat(invLambda) * B
    //   Sf_h(s).row(i) = 1.0/S2_inv.t();
    // }
    // // Rprintf("The rows of Sf_h(s) is: %d\n", Sf_h(s).n_rows);
    // // Rprintf("The rows of Phif(s) is: %d\n", Phif(s).n_elem);
    // Muf_h(s) = ((Xf(s) - repmat(bmu.col(s).t(), n_s, 1) - Muf_f(s) * A.t()) *
    //   (Bf(s) / repmat(LambdaMat.col(s), 1, q_s))) % (Sf_h(s) / repmat(Phif(s),1, q_s))  * c1; // n_s*q_s

  }
  // Rprintf("Finish update Mu_f and S_f \n");
  // double elbo2 = calELBO(Xf, bmu, LambdaMat, A, Bf, Muf_f, Sf_f, Muf_h, Sf_h, nu);
  // Rprintf("VB dF= %4f\n", elbo2 - elbo1);

  // update posterior variance of h_i, S_i
  // update posterior mean of h_i, Mu_h
  mat S2,S2_inv, Xs2;
  // Rprintf("Update Mu_h and S_h \n");
  for(s=0; s<S; ++s){ // can speed by combing two for loops
    n_s = Xf(s).n_rows;
    q_s = Bf(s).n_cols;
    S2 = Bf(s).t() * (Bf(s) / repmat(LambdaMat.col(s), 1, q_s));
    Xs2 = c1 * ((Xf(s) - repmat(bmu.col(s).t(), n_s, 1) - Muf_f(s) * A.t()) *
      (Bf(s) / repmat(LambdaMat.col(s), 1, q_s)));
    for(i=0; i<n_s; ++i){
      S2_inv =  S2 * c1 / Phif(s)(i) + eye(q_s, q_s); // diagmat(invLambda) * B
      Sf_h(s).slice(i) = S2_inv.i();
      Muf_h(s).row(i) = Xs2.row(i) * (Sf_h(s).slice(i) / Phif(s)(i))  ; // 1*q_s
    }
    // Rprintf("The rows of Sf_h(s) is: %d\n", Sf_h(s).n_rows);
    // Rprintf("The rows of Phif(s) is: %d\n", Phif(s).n_elem);
    // Muf_h(s).rows(0,4).print();

  }
  // Rprintf("Finish update Mu_h and S_h \n");
  // double elbo3 = calELBO(Xf, bmu, LambdaMat, A, Bf, Muf_f, Sf_f, Muf_h, Sf_h, nu);
  // Rprintf("VB dH= %4f\n", elbo3 - elbo2);

}

// double afun(const rowvec& x_si, const rowvec& mu_s, const mat& A, const rowvec& m_fsi,
//             const mat& B_s, const rowvec& m_hsi, const vec& lambda_vec, const mat& S_fsi,
//             const mat& S_hsi, const double& nu){
//   int q = A.n_cols, qs = B_s.n_cols;
//   rowvec dX = (x_si - mu_s- m_fsi*A.t() - m_hsi*B_s.t()) / sqrt(lambda_vec.t());
//   mat AtLA = A.t() * (A % repmat(1.0/lambda_vec, 1, q));
//   mat BsLBs = B_s.t() *(B_s % repmat(1.0/lambda_vec, 1, qs));
//   double  y = 1+ 1.0/nu * accu(dX % dX)  + 1.0/nu * (trace(AtLA*S_fsi) + trace(BsLBs*S_hsi));
//   return y;
// }

field<vec> phi_fun(const field<mat>& Xf, const mat& bmu, const mat& A, const mat& LambdaMat,
                  const double& nu, const field<mat>& Muf_f, const field<mat>& Bf,
                  const field<mat>& Muf_h,
                   const field<cube>& Sf_f, const field<cube>& Sf_h){
  int s, i, n_s, qs, S = Xf.n_elem;
  int q = A.n_cols;
  field<vec> Phif(S);

  for(s=0; s<S; ++s){
    // Rprintf("good phi_fun: s= %d!\n", s);
    n_s = Xf(s).n_rows;
    qs = Bf(s).n_cols;
    vec tmpvec(n_s, fill::zeros);
    Phif(s) = tmpvec;
    mat dX= Xf(s) - repmat(bmu.col(s).t(), n_s, 1) - Muf_f(s) * A.t() - Muf_h(s)* Bf(s).t();
    mat AtLA = A.t() * (A / repmat(LambdaMat.col(s), 1, q));
    mat BsLBs = Bf(s).t() *(Bf(s) / repmat(LambdaMat.col(s), 1, qs));
    // mat tmpMat1, tmpMat2;
    for(i = 0; i<n_s; ++i){ // Here can be speeded up!!
      //Rprintf("good afun: s=%d, i=%d!\n", s, i);
      // tmpMat1 = Sf_f(s).slice(i);
      // tmpMat2 = Sf_h(s).slice(i);
      // Phif(s)(i) = afun(Xf(s).row(i), bmu.col(s).t(), A, Muf_f(s).row(i), Bf(s),
      //      Muf_h(s).row(i), LambdaMat.col(s), tmpMat1, tmpMat2, nu);
      tmpvec(i) = 1+(trace(AtLA*Sf_f(s).slice(i)) + trace(BsLBs*Sf_h(s).slice(i))) /nu;
    }
    Phif(s) = tmpvec+ 1.0/nu*sum(dX % (dX / repmat(LambdaMat.col(s).t(), n_s, 1)), 1);

  }
  return Phif;
}
double objfun_nu(const field<mat>& Xf, const mat& bmu, const mat& A, const mat& LambdaMat,
                 const double& nu, const field<mat>& Muf_f, const field<mat>& Bf,
                 const field<mat>& Muf_h,
                 const field<cube>& Sf_f, const field<cube>& Sf_h){

  int s, n_s, S = Xf.n_elem, p = Xf(0).n_cols;
  double obj = 0.0;
  field<vec>  Phif = phi_fun( Xf, bmu, A, LambdaMat, nu, Muf_f, Bf, Muf_h, Sf_f, Sf_h);
  for(s=0; s<S; ++s){
    n_s = Xf(s).n_rows;
    // Rprintf("nu=%1f, term1 = %4f\n", nu, -accu(log(Phif(s))) * (nu+p)*0.5);
    // Rprintf("term2 = %4f\n", n_s* logCp(nu, p));
    obj += -accu(log(Phif(s))) * (nu+p)*0.5 + n_s* logCp(nu, p);

  }
  // Rprintf("obj = %4f\n", obj);
  return(obj);
}
void update_nu(const field<mat>& Xf, const mat& bmu, const mat& A, const mat& LambdaMat,
               const field<mat>& Muf_f, const field<mat>& Bf,
               const field<mat>& Muf_h,
               const field<cube>& Sf_f, const field<cube>& Sf_h, const vec& nu_grid, double& nu){
  int i, n_nu=nu_grid.n_elem;
  vec obj_vec(n_nu);
  for(i=0; i< n_nu; ++i){
    obj_vec(i) = objfun_nu( Xf, bmu, A, LambdaMat, nu_grid(i), Muf_f, Bf, Muf_h, Sf_f, Sf_h);
  }
  nu = nu_grid(index_max(obj_vec));
}

// [[Rcpp::export]]
Rcpp::List rlfm_cpp(const Rcpp::List& XList, const arma::mat& bmu_int,
                      const arma::mat& A_int, const Rcpp::List& BList_int,
                      const double& nu_int, const arma::mat& LambdaMat_int,
                      const Rcpp::List& SList_f_int, const Rcpp::List& MuList_f_int,
                      const Rcpp::List& SList_h_int, const Rcpp::List& MuList_h_int,
                      const double& epsELBO, const int& maxIter, const bool& verbose,
                      const bool& loop_ic = false){


  int s, S = XList.length(); // get the number of data source
  vec nu_grid = round(linspace(3, 100, 100)); //  round(linspace(1, 100, 100))
  // Initialize
  // transfer list to field.
  field<mat> Xf(S),  Muf_f(S),  Muf_h(S), Bf(S);
  field<cube> Sf_f(S),  Sf_h(S); // each S_f_si is a diagonal matrix
  field<vec> af(S);
  // Rprintf("good starting!\n");
  for(s=0; s < S; ++s){
    mat Xtmp1 = XList[s];
    Xf(s) = Xtmp1;
    mat Xtmp3 = MuList_f_int[s];
    Muf_f(s) = Xtmp3;
    cube Xtmp4 = SList_f_int[s];
    Sf_f(s) = Xtmp4;
    mat Xtmp5 = MuList_h_int[s];
    Muf_h(s) = Xtmp5;
    cube Xtmp6 = SList_h_int[s];
    Sf_h(s) = Xtmp6;
    mat tmp = BList_int[s];
    Bf(s) = tmp;
  }

  mat bmu(bmu_int), LambdaMat(LambdaMat_int), A(A_int);
  double nu(nu_int);

  add_IC_Orth(A, Bf);
  //Rprintf("good starting2!\n");
  vec ELBO_vec(maxIter);
  ELBO_vec(0) = INT_MIN1;
  int iter;
  field<vec> Phif;
  for(iter = 1; iter < maxIter; ++iter){


    // Rprintf("Compute Phif!\n");
    Phif = phi_fun( Xf, bmu, A, LambdaMat, nu, Muf_f, Bf, Muf_h, Sf_f, Sf_h);
    // Phif(0).print();
    //Rprintf("Start Estep...");
    // VB E-step
    VB_EstepVI( Xf, bmu, LambdaMat, A, Bf, Muf_f,Sf_f,Muf_h,Sf_h, nu, Phif);
    //Rprintf(", E step Finished!\n");

    //VB M-step
    // double elbo1 = calELBO(Xf, bmu, LambdaMat, A, Bf, Muf_f, Sf_f, Muf_h, Sf_h, nu);
    // Rprintf("ELBO1 = %4f!\n", elbo1);

    Phif = phi_fun( Xf, bmu, A, LambdaMat, nu, Muf_f, Bf, Muf_h, Sf_f, Sf_h); // Update Phif.
    //update bmu
    // Rprintf("Update mu!\n");
    update_mu( Xf, A, Bf, Muf_f, Muf_h, Phif, bmu);
    // double elbo2 = calELBO(Xf, bmu, LambdaMat, A, Bf, Muf_f, Sf_f, Muf_h, Sf_h, nu);
    // Rprintf("dmu= %4f \n", elbo2 - elbo1);

    //update A
    // Rprintf("Update A!\n");
    update_A( Xf, bmu, LambdaMat, Bf, Muf_f, Muf_h, Sf_f, Phif, A);
    // A.row(0).print();
    // double elbo3 = calELBO(Xf, bmu, LambdaMat, A, Bf, Muf_f, Sf_f, Muf_h, Sf_h, nu);
    // Rprintf("dA= %4f \n", elbo3 - elbo2);
    // update B

    // Rprintf("Update B!\n");
    update_B( Xf, bmu, LambdaMat, A, Muf_f, Muf_h, Sf_h, Phif, Bf);
    // double elbo4 = calELBO(Xf, bmu, LambdaMat, A, Bf, Muf_f, Sf_f, Muf_h, Sf_h, nu);
    // Rprintf("dB= %4f \n", elbo4 - elbo3);
    // Add identifiability conditions
    if(loop_ic){
      add_IC_Orth(A, Bf);
    }


    // update Lambda
    // Rprintf("Update Lambda!\n");
    update_Lambda( Xf, bmu, A, Bf, Muf_f, Sf_f, Muf_h, Sf_h, Phif, nu,LambdaMat);
    // double elbo5 =  calELBO(Xf, bmu, LambdaMat, A, Bf, Muf_f, Sf_f, Muf_h, Sf_h, nu);
    // Rprintf("dLambda= %4f\n", elbo5 - elbo4);

    // To speed up the computation, we can update nu at final step
    // Rprintf("Update nu!\n");
    update_nu( Xf, bmu, A, LambdaMat, Muf_f, Bf, Muf_h, Sf_f, Sf_h, nu_grid, nu);
    // double elbo6 =  calELBO(Xf, bmu, LambdaMat, A, Bf, Muf_f, Sf_f, Muf_h, Sf_h, nu);
    // Rprintf("dnu= %4f \n", elbo6 - elbo5);

    ELBO_vec(iter) = calELBO( Xf, bmu, LambdaMat, A, Bf, Muf_f, Sf_f, Muf_h, Sf_h, nu);


    if(verbose){
      Rprintf("iter = %d, ELBO= %4f, dELBO=%4f \n",
              iter +1, ELBO_vec(iter), abs(ELBO_vec(iter)  - ELBO_vec(iter-1))/ abs(ELBO_vec(iter-1)));
    }
    if(abs((ELBO_vec(iter)  - ELBO_vec(iter-1))/ ELBO_vec(iter-1)) < epsELBO) break;


  }

  // Add identifiability conditions
  if(!loop_ic){
    add_IC_Orth(A, Bf);
    // re-estimate the parameters.
    VB_EstepVI( Xf, bmu, LambdaMat, A, Bf, Muf_f,Sf_f,Muf_h,Sf_h, nu, Phif);
  }
  // update nu: TODO
  // Rprintf("Update nu!\n");
  // update_nu( Xf, bmu, A, LambdaMat, Muf_f, Bf, Muf_h, Sf_f, Sf_h, nu_grid, nu);
  // double elbo6 =  calELBO(Xf, bmu, LambdaMat, A, Bf, Muf_f, Sf_f, Muf_h, Sf_h, nu);
  // Rprintf("dnu= %4f \n", elbo6 - elbo5);

  // output return value
  List resList = List::create(
    Rcpp::Named("F") = Muf_f,
    Rcpp::Named("H") = Muf_h,
    Rcpp::Named("Sf") =Sf_f,
    Rcpp::Named("Sh") = Sf_h,
    Rcpp::Named("mu") = bmu,
    Rcpp::Named("A") = A,
    Rcpp::Named("B") = Bf,
    Rcpp::Named("nu") = nu,
    Rcpp::Named("Lambda") = LambdaMat,
    Rcpp::Named("ELBO") = ELBO_vec(iter-1),
    Rcpp::Named("dELBO") = ELBO_vec(iter-1)  - ELBO_vec(iter-2),
    Rcpp::Named("ELBO_seq") = ELBO_vec.subvec(0, iter-1)
  );
  return(resList);

}
