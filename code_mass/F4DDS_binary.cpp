// Rcpp code
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double F4DDS_binary_cpp(NumericMatrix B, NumericMatrix Xn, NumericVector Yn, double hn) {
  int d = B.ncol();
  IntegerVector idx1;
  IntegerVector idx2;
  for (int i = 0; i < Yn.length(); i++) {
    if (Yn[i] == 0) {
      idx1.push_back(i);
    } else if (Yn[i] == 1) {
      idx2.push_back(i);
    }
  }

  NumericMatrix X1(idx1.size(), d);
  NumericMatrix X2(idx2.size(), d);
  for (int i = 0; i < idx1.size(); i++) {
    for (int j = 0; j < d; j++) {
      X1(i, j) = sum(Xn(idx1[i], _) * B(_, j)); // Assuming B is a column vector
    }
  }
  for (int i = 0; i < idx2.size(); i++) {
    for (int j = 0; j < d; j++) {
      X2(i, j) = sum(Xn(idx2[i], _) * B(_, j)); // Assuming B is a column vector
    }
  }

  int n1 = X1.nrow();
  int n2 = X2.nrow();
  int n = n1 + n2;
  double p1 = static_cast<double>(n1) / n;
  double p2 = static_cast<double>(n2) / n;
  NumericMatrix Zn(Xn.nrow(), d);
  
  for (int i = 0; i < Xn.nrow(); i++) {
    for (int j = 0; j < d; j++) {
      Zn(i, j) = 0.0;
      for (int k = 0; k < Xn.ncol(); k++) {
        Zn(i, j) += Xn(i, k) * B(k, j);
      }
    }
  }
  
  NumericMatrix W(n, n);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      if (i != j) {
        double norm_sq = 0.0;
        
        for (int k = 0; k < d; k++) {
          double diff = Zn(i, k) - Zn(j, k);
          norm_sq += diff * diff;
        }
        
        double exponent = -0.5 * pow(hn, -2 * d) * norm_sq;
        W(i, j) = exp(exponent);
      }
    }
  }

  NumericVector f1hat(n);
  NumericVector f2hat(n);
  for (int i = 0; i < n; i++) {
    f1hat[i] = 0;
    f2hat[i] = 0;
    for (int j = 0; j < n1; j++) {
      f1hat[i] += W(idx1[j], i);
    }
    for (int j = 0; j < n2; j++) {
      f2hat[i] += W(idx2[j], i);
    }
    f1hat[i] *= 1.0 / (n1 - 1) * pow(hn, -d) * pow(2.0 * M_PI, -0.5 * d);
    f2hat[i] *= 1.0 / (n2 - 1) * pow(hn, -d) * pow(2.0 * M_PI, -0.5 * d);
  }
  
  double f = 0.0;
  for (int i = 0; i < n; i++) {
    if (p1 * f1hat[i] + p2 * f2hat[i] != 0) {
      f += sqrt(f1hat[i] * f2hat[i]) / (p1 * f1hat[i] + p2 * f2hat[i]);
    }
  }
  
  
  
  return f;
}
