#ifndef _HEADER_
#include "cgalMesh.h"
#endif


bool isEqual(double l, double r) {
  return r == std::nextafter(l, r);
}


// [[Rcpp::export]]
Rcpp::List gatherVertices(
    Rcpp::NumericMatrix Vertices, Rcpp::IntegerMatrix Faces
) {
  int nvertices = Vertices.ncol();
  std::vector<bool> duplicated(nvertices, false);
  std::map<int, int> newindices;
  int newindex = 1;
  
  for(int index = 0; index < nvertices - 1; index++) {
    if(!duplicated[index]) {
      newindices[index + 1] = newindex;
//      Rcpp::NumericVector vertex = Vertices(Rcpp::_, index);
      const Rcpp::NumericMatrix::Column& vertex = Vertices.column(index);
      for(int j = index + 1; j < nvertices; j++) {
        //Rcpp::NumericVector vertex_j = Vertices(Rcpp::_, j);
        const Rcpp::NumericMatrix::Column& vertex_j = Vertices.column(j);
        bool equal = Rcpp::max(Rcpp::abs(vertex - vertex_j)) < 1e-16;
        // bool equal = isEqual(vertex(0), vertex_j(0)) && 
        //              isEqual(vertex(1), vertex_j(1)) && 
        //              isEqual(vertex(2), vertex_j(2));
        if(equal) {
          duplicated[j] = true;
          newindices[j + 1] = newindex;
        }
      }
      newindex++;
    }
  }
  
  if(!duplicated[nvertices - 1]) {
    newindices[nvertices] = newindex++;
  }
  
  Rcpp::NumericMatrix NewVertices(3, newindex-1);
  int j = 0;
  for(auto const& [key, value] : newindices) {
    if(!duplicated[key - 1]) {
      NewVertices(Rcpp::_, j++) = Vertices(Rcpp::_, key - 1);
    }
  }
  
  int nfaces = Faces.ncol();
  Rcpp::IntegerMatrix NewFaces(3, nfaces);
  for(int i = 0; i < nfaces; i++) {
    Rcpp::IntegerVector face = Faces(Rcpp::_, i);
    Rcpp::IntegerVector newface = 
      {newindices[face(0)], newindices[face(1)], newindices[face(2)]};
    NewFaces(Rcpp::_, i) = newface;
  }
  
  return Rcpp::List::create(
    Rcpp::Named("vertices") = NewVertices,
    Rcpp::Named("faces")    = NewFaces
  );
}