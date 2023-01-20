#ifndef _HEADER_
#include "cgalMesh.h"
#endif

// [[Rcpp::export]]
Rcpp::IntegerMatrix meshTopology(
  int nu, int nv, bool uperiodic, bool vperiodic
) {
  if(uperiodic && !vperiodic) {
    int shift = nu * (nv - 1);
    Rcpp::IntegerMatrix triangles(3, 2 * shift);
    for(int i = 1; i <= nu; i++) {
      int im1nv = (i - 1) * nv;
      int ip1m1nv = i == nu ? 0 : i * nv;
      int k = (i - 1) * (nv - 1);
      for(int j = 1; j <= nv-1; j++) {
        int jp1 = j + 1;
        Rcpp::IntegerVector tri1 = {im1nv + j, im1nv + jp1, ip1m1nv + j};
        Rcpp::IntegerVector tri2 = {im1nv + jp1, ip1m1nv + jp1, ip1m1nv + j};
        triangles(Rcpp::_, k + j - 1) = tri1;
        triangles(Rcpp::_, k + j - 1 + shift) = tri2;
      }
    }
    return triangles;
  }

  if(uperiodic && vperiodic) {
    int shift = nu * nv;
    Rcpp::IntegerMatrix triangles(3, 2 * shift);
    for(int i = 1; i <= nu; i++) {
      int im1nv = (i - 1) * nv;
      int ip1m1nv = i == nu ? 0 : i * nv;
      int k = (i - 1) * nv;
      for(int j = 1; j <= nv; j++) {
        int jp1 = j == nv ? 1 : j + 1;
        Rcpp::IntegerVector tri1 = {im1nv + j, im1nv + jp1, ip1m1nv + j};
        Rcpp::IntegerVector tri2 = {im1nv + jp1, ip1m1nv + jp1, ip1m1nv + j};
        triangles(Rcpp::_, k + j - 1) = tri1;
        triangles(Rcpp::_, k + j - 1 + shift) = tri2;
      }
    }
    return triangles;
  }

  if(!uperiodic && vperiodic) {
    int shift = (nu - 1) * nv;
    Rcpp::IntegerMatrix triangles(3, 2 * shift);
    for(int i = 1; i <= nu-1; i++) {
      int im1nv = (i - 1) * nv;
      int ip1m1nv = i * nv;
      int k = (i - 1) * nv;
      for(int j = 1; j <= nv; j++) {
        int jp1 = j == nv ? 1 : j + 1;
        Rcpp::IntegerVector tri1 = {im1nv + j, im1nv + jp1, ip1m1nv + j};
        Rcpp::IntegerVector tri2 = {im1nv + jp1, ip1m1nv + jp1, ip1m1nv + j};
        triangles(Rcpp::_, k + j - 1) = tri1;
        triangles(Rcpp::_, k + j - 1 + shift) = tri2;
      }
    }
    return triangles;
  }

  int shift = (nu - 1) * (nv - 1);
  Rcpp::IntegerMatrix triangles(3, 2 * shift);
  for(int i = 1; i <= nu-1; i++) {
    int im1nv = (i - 1) * nv;
    int ip1m1nv = i * nv;
    int k = (i - 1) * (nv - 1);
    for(int j = 1; j <= nv-1; j++) {
      int jp1 = j + 1;
      Rcpp::IntegerVector tri1 = {im1nv + j, im1nv + jp1, ip1m1nv + j};
      Rcpp::IntegerVector tri2 = {im1nv + jp1, ip1m1nv + jp1, ip1m1nv + j};
      triangles(Rcpp::_, k + j - 1) = tri1;
      triangles(Rcpp::_, k + j - 1 + shift) = tri2;
    }
  }
  return triangles;

}
