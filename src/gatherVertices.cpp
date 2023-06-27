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
Rcpp::List gatherVertices(Rcpp::NumericMatrix Vertices) {
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

  Rcpp::IntegerVector NewIndices(nvertices);
  int newindex = 1;
  
  for(int index = 0; index < nvertices; index++) {
    if(duplicated[index] == 0) {
      NewIndices(index) = newindex++;
    } else {
      NewIndices(index) = duplicated[index];
    }
  }

  Rcpp::IntegerVector Extraction(newindex - 1);
  int j = 0;
  for(int i = 0; i < nvertices; i++) {
    if(duplicated[i] == 0) {
      Extraction(j++) = i + 1;
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
  
  return Rcpp::List::create(
    Rcpp::Named("extraction") = Extraction,
    Rcpp::Named("newindices") = NewIndices
  );
}


// -------------------------------------------------------------------------- //
// -------------------------------------------------------------------------- //
// [[Rcpp::export]]
Rcpp::IntegerVector facesToDelete(Rcpp::IntegerMatrix Faces) {
  int nfaces = Faces.ncol();
  std::vector<int> todelete;
  for(int j = 0; j < nfaces; j++) {
    Rcpp::IntegerVector F = Faces(Rcpp::_, j);
    if(F(0) == F(1) || F(0) == F(2) || F(1) == F(2)) {
      todelete.push_back(j + 1);
    }
  }
  Rcpp::IntegerVector ToDelete(todelete.begin(), todelete.end());
  return ToDelete;
}