#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix dF4DDS_pairwise_cpp(NumericMatrix B, NumericMatrix Xn, NumericVector Yn, double hn) {

  int K = max(Yn);
  int d = B.ncol();
  int p = Xn.ncol();
  NumericMatrix df(p, d);
  df.fill(0);
  for (int km = 1; km <= K-1; km++) {
    IntegerVector idx10;
    IntegerVector idx20;
    for (int i = 0; i < Yn.length(); i++) {
      if (Yn[i] == km) {
        idx10.push_back(i);
      } else if (Yn[i] == km+1) {
        idx20.push_back(i);
      }
    }
    
    NumericMatrix Xnew(idx10.size() + idx20.size(), p);
    NumericMatrix Ynew(idx10.size() + idx20.size(), 1);
    int counter = 0;
    for (int i = 0; i < idx10.size(); i++) {
      for (int j = 0; j < p; j++) {
        Xnew(counter, j) = Xn(idx10[i] , j);
      }
      Ynew[counter] = Yn[idx10[i] ];
      counter++;
    }
    
    for (int i = 0; i < idx20.size(); i++) {
      for (int j = 0; j < p; j++) {
        Xnew(counter, j) = Xn(idx20[i], j);
      }
      Ynew[counter] = Yn[idx20[i]];
      counter++;
    }
    
    IntegerVector idx1;
    IntegerVector idx2;
    
    for (int i = 0; i < Ynew.size(); i++) {
      if (Ynew[i] == km) {
        idx1.push_back(i);
      }
      if (Ynew[i] == km+1) {
        idx2.push_back(i);
      }
    }
    
    NumericMatrix X1(idx1.size(), d);
    NumericMatrix X2(idx2.size(), d);
    for (int i = 0; i < idx1.size(); i++) {
      for (int j = 0; j < d; j++) {
        X1(i, j) = sum(Xnew(idx1[i], _) * B(_, j)); // Assuming B is a column vector
      }
    }
    for (int i = 0; i < idx2.size(); i++) {
      for (int j = 0; j < d; j++) {
        X2(i, j) = sum(Xnew(idx2[i], _) * B(_, j)); // Assuming B is a column vector
      }
    }
    int n1 = X1.nrow();
    int n2 = X2.nrow();
    int n = n1 + n2;
    double p1 = static_cast<double>(n1) / n;
    double p2 = static_cast<double>(n2) / n;
    NumericMatrix Znew(Xnew.nrow(), d);
    
    for (int i = 0; i < Xnew.nrow(); i++) {
      for (int j = 0; j < d; j++) {
        Znew(i, j) = 0.0;
        for (int k = 0; k < Xn.ncol(); k++) {
          Znew(i, j) += Xnew(i, k) * B(k, j);
        }
      }
    }
    
    NumericMatrix W(n, n);
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        if (i != j) {
          double norm_sq = 0.0;
          for (int k = 0; k < d; k++) {
            double diff = Znew(i, k) - Znew(j, k);
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
    
    NumericVector W1(n);
    NumericVector W2(n);
    NumericMatrix Wt = clone(W);
    
    for (int j = 0; j < n; j++) {
      W1[j] = 0.0;
      for (int i = 0; i < n1; i++) {
        W1[j] += W(idx1[i], j);
      }
    }
    
    
    for (int j = 0; j < n; j++) {
      W2[j] = 0.0;
      for (int i = 0; i < n2; i++) {
        W2[j] += W(idx2[i], j);
      }
    }
    
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        if (Yn[j] == km) {
          Wt(i, j) = W(i, j) / W1[i] * pow(hn, -2 * d);
        } else {
          Wt(i, j) = -W(i, j) / W2[i] * pow(hn, -2 * d);
        }
      }
    }
    
    for (int i = 0; i < n; i++) {
      NumericMatrix Xik(n, p);
      for (int j = 0; j < n; j++) {
        for (int k = 0; k < p; k++) {
          Xik(j, k) = Xnew(j, k) - Xnew(i, k);
        }
      }
      
      NumericMatrix Sik(p, p);
      Sik.fill(0);
      for (int j = 0; j < p; j++) {
        for (int k = 0; k < p; k++) {
          for (int l = 0; l < n; l++) {
            Sik(j, k) += Xik(l, j) * Wt(i, l) * Xik(l, k);
          }
        }
      }
      
      for (int j = 0; j < p; j++) {
        for (int k = 0; k < d; k++) {
          for (int l=0; l<p; l++){
            df(j, k) += (p1 * f1hat[i] - p2 * f2hat[i]) * sqrt(f1hat[i] * f2hat[i]) / (2 * pow(p1 * f1hat[i] + p2 * f2hat[i], 2)) * Sik(j, l) * B(l, k);
          }
        }  
      }
    }
  }

  return df;  
}
  
  