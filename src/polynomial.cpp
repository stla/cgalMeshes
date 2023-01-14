#ifndef _HEADER_
#include "cgalMesh.h"
#endif

#include <CGAL/Polynomial.h>
#include <CGAL/Polynomial_traits_d.h>
#include <CGAL/Polynomial_type_generator.h>

typedef CGAL::Polynomial_type_generator<double, 3>::Type Poly3;
typedef CGAL::Polynomial_traits_d<Poly3>                 PT3;
//typedef PT_2::Coefficient_type                       Poly_1;
typedef PT3::Innermost_coefficient_type                  Real;

// [[Rcpp::export]]
void polynomial(Rcpp::IntegerMatrix powers, Rcpp::NumericVector coeffs, double x, double y, double z) {
  CGAL::IO::set_pretty_mode(std::cout);

  PT3::Construct_polynomial construct_polynomial;

  std::list<std::pair<CGAL::Exponent_vector, Real>> innermost_coeffs;
  for(int i = 0; i < coeffs.size(); i++) {
    Rcpp::IntegerVector pows = powers(i, Rcpp::_);
    innermost_coeffs.push_back(
      std::make_pair(CGAL::Exponent_vector(pows(0), pows(1), pows(2)), coeffs(i))
    );
  }
  Poly3 G = 
    construct_polynomial(innermost_coeffs.begin(),innermost_coeffs.end());
  std::cout << "The trivariate polynomial G: " << G << std::endl;

  // evaluation
  PT3::Substitute substitute;
  std::list<Real> replacements;
  replacements.push_back(x); 
  replacements.push_back(y);
  replacements.push_back(z);

  std::cout << "G(x,y,z): "  
            << substitute(G,replacements.begin(),replacements.end())
            << std::endl;

}
