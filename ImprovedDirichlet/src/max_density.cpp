#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <cmath>
using namespace Rcpp;
using namespace arma;
using namespace Eigen;

// Digamma and trigamma functions
inline double digamma(double x) {
    return R::digamma(x);
}

inline double trigamma(double x) {
    return R::trigamma(x);
}

// [[Rcpp::export]]
NumericVector max_density_for_beta_cpp(double c, double v = -1, double a_init = -1, double b_init = -1,
                                       double tol = 1e-8, int maxiter = 100, double stepsize = 0.5,
                                       bool verbose = false, double alpha = -1) {

    if (a_init == -1)  a_init = c;
    if (b_init == -1)  b_init = 1 - c;

    for (int attempt = 0; attempt < 10; attempt++) {
        double a = a_init, b = b_init;
        double a_old, b_old, f;

        for (int iter = 0; iter < maxiter; iter++) {
            double df_a = digamma(a) - digamma(a + b) - log(c);
            double df_b = digamma(b) - digamma(a + b) + log(1 - c);

            double h, dh_a, dh_b;
            if (alpha == -1) {
                h = log(a) + log(b) - 2 * log(a + b) - log(a + b + 1) - log(v);
                dh_a = 1 / a - 2 / (a + b) - 1 / (a + b + 1);
                dh_b = 1 / b - 2 / (a + b) - 1 / (a + b + 1);
            } else {
                h = (a + b) / alpha - 1;
                dh_a = dh_b = 1 / alpha;
            }

            // Hessian matrix
            double df_aa = trigamma(a) - trigamma(a + b);
            double df_bb = trigamma(b) - trigamma(a + b);
            double df_ab = - trigamma(a + b);

            mat H = {{df_aa, df_ab}, {df_ab, df_bb}};
            rowvec J = {dh_a, dh_b};

            mat A = join_vert(H, J);
            A.insert_cols(2, 1);
            A.row(2) = join_horiz(J, rowvec(1));

            vec y = {-df_a, -df_b, -h};
            vec x = solve(A, y);

            if (verbose) {
                f = - (a - 1) * log(c) - (b - 1) * log(1 - c) - lgamma(a + b) + lgamma(a) + lgamma(b);
                Rprintf("iter = %d \t a = %.10f \t b = %.10f \t f = %.10f \t h = %.10f\n", iter, a, b, f, h);
            }

            a_old = a;
            b_old = b;

            a += stepsize * x[0];
            b += stepsize * x[1];

            if (a < 0)  a = a_old / 2;
            if (b < 0)  b = b_old / 2;

            if (std::abs(a / a_old - 1) + std::abs(b / b_old - 1) + std::abs(h) < tol)  return NumericVector::create(a, b);
        }
        stepsize /= 5;
        maxiter *= 5;
        if (verbose)  Rprintf("Failed to converge. Retrying with step.size = %.5f and max.iter = %d", stepsize, maxiter);
    }
    stop("Failed to converge after all attempts.");
}

// // [[Rcpp::export]]
// NumericVector max_density_for_dirichlet_cpp(NumericVector c, double theta, std::string version,
//                                             NumericVector a_init = NumericVector(), double tol = 1e-8,
//                                             int maxiter = 100, double stepsize = 0.5, bool verbose = false) {
//
//     int k = c.size();
//
//     if (a_init.size() == 0) {
//         a_init = clone(c);
//     }
//
//     for (int attempt = 0; attempt < 5; attempt++) {
//         NumericVector a = clone(a_init);
//         NumericVector a_old(k);
//
//         for (int iter = 0; iter < maxiter; iter++) {
//             NumericVector g(k);
//             for (int i = 0; i < k; ++i) {
//                 g[i] = digamma(a[i]) - digamma(sum(a)) - log(c[i]);
//             }
//
//             NumericMatrix H(k, k);
//             for (int i = 0; i < k; i++) {
//                 H(i, i) = trigamma(a[i]);
//                 for (int j = 0; j < k; j++) {
//                     H(i, j) -= trigamma(sum(a));
//                 }
//             }
//
//             double h;
//             NumericVector J(k);
//             if (version == "concentration") {
//                 h = sum(a) / theta - 1;
//                 for (int i = 0; i < k; ++i) J[i] = 1 / theta;
//             }
//
//             NumericMatrix A(k + 1, k + 1);
//             for (int i = 0; i < k; ++i) {
//                 for (int j = 0; j < k; ++j) {
//                     A(i, j) = H(i, j);
//                 }
//                 A(i, k) = J[i];
//                 A(k, i) = J[i];
//             }
//             A(k, k) = 0;
//
//             NumericVector y(k + 1);
//             for (int i = 0; i < k; ++i) y[i] = -g[i];
//             y[k] = -h;
//
//             NumericVector x = Rcpp::solve(A, y);
//
//             a_old = clone(a);
//             for (int i = 0; i < k; ++i)  a[i] += stepsize * x[i];
//
//             for (int i = 0; i < k; ++i) {
//                 if (a[i] < 0)  a[i] = a_old[i] / 2;
//             }
//
//             if (sum(abs(a / a_old - 1)) + std::abs(h) < tol)  return a;
//         }
//         stepsize /= 5;
//         maxiter *= 5;
//     }
//     stop("Failed to converge after all attempts.");
// }
