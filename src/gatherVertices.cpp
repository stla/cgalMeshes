#ifndef _HEADER_
#include "cgalMesh.h"
#endif

typedef std::tuple<double, double, double> Vertex3;


// -------------------------------------------------------------------------- //
// author: Mikko Marttila --------------------------------------------------- //
std::vector<std::vector<int>> duplicatedIndices(Rcpp::NumericMatrix Vertices) {
  std::map<Vertex3, std::vector<int>> positions;
  std::set<Vertex3> duplicates;
  
  const int nvertices = Vertices.ncol(); 
  
  for(int i = 0; i < nvertices; i++) {
    Rcpp::NumericVector vi = Vertices(Rcpp::_, i);
    Vertex3 vertex = { vi(0), vi(1), vi(2) };
    if(positions.count(vertex) == 0) {
      std::vector<int> position = { i + 1 };
      positions.emplace(vertex, position);
    } else {
      duplicates.insert(vertex);
      positions.at(vertex).push_back(i + 1);
    }
  }
  
  std::vector<std::vector<int>> result;
  int nduplicates = duplicates.size();
  result.reserve(nduplicates);
  for(const Vertex3& vertex : duplicates) {
    result.push_back(positions.at(vertex));
  }
  
  return result;
}


// bool isEqual(double l, double r) {
//   return r == std::nextafter(l, r);
// }


// -------------------------------------------------------------------------- //
// -------------------------------------------------------------------------- //
// [[Rcpp::export]]
Rcpp::List gatherVertices(
    Rcpp::NumericMatrix Vertices, Rcpp::IntegerMatrix Faces
) {
  std::vector<std::vector<int>> dupIndices = duplicatedIndices(Vertices);
  const int ndups = dupIndices.size();
  const int nvertices = Vertices.ncol(); 
  std::vector<int> duplicated(nvertices, 0);
  
  for(int i = 0; i < ndups; i++) {
    std::vector<int> indices = dupIndices[i];
    const int nindices = indices.size();
    for(int j = 1; j < nindices; j++) {
      duplicated[indices[j]-1] = indices[0];
    }
  }

  std::map<int, int> newindices;
  int newindex = 1;
  
  for(int index = 0; index < nvertices; index++) {
    if(duplicated[index] == 0) {
      newindices[index] = newindex++;
    } else {
      newindices[index] = duplicated[index];
    }
  }

  Rcpp::IntegerVector NewIndices(newindex - 1);
  Rcpp::NumericMatrix NewVertices(3, newindex - 1);
  int j = 0;
  for(auto const& [key, value] : newindices) {
    if(duplicated[key] == 0) {
      NewIndices(j) = key + 1;
      NewVertices(Rcpp::_, j++) = Vertices(Rcpp::_, key);
    }
  }

  int nremoved = nvertices - newindex + 2;
  std::string msg;
  if(nremoved == 0) {
    msg = "No duplicated vertex.";
  } else if(nremoved == 1) {
    msg = "One duplicated vertex has been removed.";
  } else {
    msg = std::to_string(nremoved) + " duplicated vertices removed.";
  }
  Message(msg);
  
  int nfaces = Faces.ncol();
  Rcpp::IntegerMatrix NewFaces(3, nfaces);
  for(int i = 0; i < nfaces; i++) {
    Rcpp::IntegerVector face = Faces(Rcpp::_, i);
    Rcpp::IntegerVector newface = 
      {newindices[face(0)-1], newindices[face(1)-1], newindices[face(2)-1]};
    NewFaces(Rcpp::_, i) = newface;
  }
  
  return Rcpp::List::create(
    Rcpp::Named("vertices") = NewVertices,
    Rcpp::Named("faces")    = NewFaces,
    Rcpp::Named("indices")  = NewIndices
  );
}
