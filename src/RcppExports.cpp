// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "cgalMeshes_types.h"
#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// AFSreconstruction_cpp
Rcpp::XPtr<EMesh3> AFSreconstruction_cpp(const Rcpp::NumericMatrix pts, const unsigned nneighs);
RcppExport SEXP _cgalMeshes_AFSreconstruction_cpp(SEXP ptsSEXP, SEXP nneighsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix >::type pts(ptsSEXP);
    Rcpp::traits::input_parameter< const unsigned >::type nneighs(nneighsSEXP);
    rcpp_result_gen = Rcpp::wrap(AFSreconstruction_cpp(pts, nneighs));
    return rcpp_result_gen;
END_RCPP
}
// cxhull
Rcpp::XPtr<EMesh3> cxhull(Rcpp::NumericMatrix pts);
RcppExport SEXP _cgalMeshes_cxhull(SEXP ptsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type pts(ptsSEXP);
    rcpp_result_gen = Rcpp::wrap(cxhull(pts));
    return rcpp_result_gen;
END_RCPP
}
// cxhullsIntersection
Rcpp::XPtr<EMesh3> cxhullsIntersection(Rcpp::List Pts, Rcpp::Nullable<Rcpp::NumericVector> origin_);
RcppExport SEXP _cgalMeshes_cxhullsIntersection(SEXP PtsSEXP, SEXP origin_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type Pts(PtsSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericVector> >::type origin_(origin_SEXP);
    rcpp_result_gen = Rcpp::wrap(cxhullsIntersection(Pts, origin_));
    return rcpp_result_gen;
END_RCPP
}
// gTriangle
Rcpp::XPtr<EMesh3> gTriangle(Rcpp::NumericVector A, Rcpp::NumericVector B, Rcpp::NumericVector C, double s, int iterations);
RcppExport SEXP _cgalMeshes_gTriangle(SEXP ASEXP, SEXP BSEXP, SEXP CSEXP, SEXP sSEXP, SEXP iterationsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type A(ASEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type B(BSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type C(CSEXP);
    Rcpp::traits::input_parameter< double >::type s(sSEXP);
    Rcpp::traits::input_parameter< int >::type iterations(iterationsSEXP);
    rcpp_result_gen = Rcpp::wrap(gTriangle(A, B, C, s, iterations));
    return rcpp_result_gen;
END_RCPP
}
// AlgebraicMesh
Rcpp::XPtr<EMesh3> AlgebraicMesh(Rcpp::IntegerMatrix powers, Rcpp::NumericVector coeffs, double isolevel, Rcpp::NumericVector sphereCenter, double sphereRadius, double angle_bound, double radius_bound, double distance_bound, double error_bound);
RcppExport SEXP _cgalMeshes_AlgebraicMesh(SEXP powersSEXP, SEXP coeffsSEXP, SEXP isolevelSEXP, SEXP sphereCenterSEXP, SEXP sphereRadiusSEXP, SEXP angle_boundSEXP, SEXP radius_boundSEXP, SEXP distance_boundSEXP, SEXP error_boundSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::IntegerMatrix >::type powers(powersSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type coeffs(coeffsSEXP);
    Rcpp::traits::input_parameter< double >::type isolevel(isolevelSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type sphereCenter(sphereCenterSEXP);
    Rcpp::traits::input_parameter< double >::type sphereRadius(sphereRadiusSEXP);
    Rcpp::traits::input_parameter< double >::type angle_bound(angle_boundSEXP);
    Rcpp::traits::input_parameter< double >::type radius_bound(radius_boundSEXP);
    Rcpp::traits::input_parameter< double >::type distance_bound(distance_boundSEXP);
    Rcpp::traits::input_parameter< double >::type error_bound(error_boundSEXP);
    rcpp_result_gen = Rcpp::wrap(AlgebraicMesh(powers, coeffs, isolevel, sphereCenter, sphereRadius, angle_bound, radius_bound, distance_bound, error_bound));
    return rcpp_result_gen;
END_RCPP
}
// AlgebraicMeshesIntersection
Rcpp::XPtr<EMesh3> AlgebraicMeshesIntersection(Rcpp::List Rpolynomials, Rcpp::NumericVector sphereCenter, double sphereRadius, double angle_bound, double radius_bound, double distance_bound, double error_bound);
RcppExport SEXP _cgalMeshes_AlgebraicMeshesIntersection(SEXP RpolynomialsSEXP, SEXP sphereCenterSEXP, SEXP sphereRadiusSEXP, SEXP angle_boundSEXP, SEXP radius_boundSEXP, SEXP distance_boundSEXP, SEXP error_boundSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type Rpolynomials(RpolynomialsSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type sphereCenter(sphereCenterSEXP);
    Rcpp::traits::input_parameter< double >::type sphereRadius(sphereRadiusSEXP);
    Rcpp::traits::input_parameter< double >::type angle_bound(angle_boundSEXP);
    Rcpp::traits::input_parameter< double >::type radius_bound(radius_boundSEXP);
    Rcpp::traits::input_parameter< double >::type distance_bound(distance_boundSEXP);
    Rcpp::traits::input_parameter< double >::type error_bound(error_boundSEXP);
    rcpp_result_gen = Rcpp::wrap(AlgebraicMeshesIntersection(Rpolynomials, sphereCenter, sphereRadius, angle_bound, radius_bound, distance_bound, error_bound));
    return rcpp_result_gen;
END_RCPP
}
// AlgebraicMeshesUnion
Rcpp::XPtr<EMesh3> AlgebraicMeshesUnion(Rcpp::List Rpolynomials, Rcpp::NumericVector sphereCenter, double sphereRadius, double angle_bound, double radius_bound, double distance_bound, double error_bound);
RcppExport SEXP _cgalMeshes_AlgebraicMeshesUnion(SEXP RpolynomialsSEXP, SEXP sphereCenterSEXP, SEXP sphereRadiusSEXP, SEXP angle_boundSEXP, SEXP radius_boundSEXP, SEXP distance_boundSEXP, SEXP error_boundSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type Rpolynomials(RpolynomialsSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type sphereCenter(sphereCenterSEXP);
    Rcpp::traits::input_parameter< double >::type sphereRadius(sphereRadiusSEXP);
    Rcpp::traits::input_parameter< double >::type angle_bound(angle_boundSEXP);
    Rcpp::traits::input_parameter< double >::type radius_bound(radius_boundSEXP);
    Rcpp::traits::input_parameter< double >::type distance_bound(distance_boundSEXP);
    Rcpp::traits::input_parameter< double >::type error_bound(error_boundSEXP);
    rcpp_result_gen = Rcpp::wrap(AlgebraicMeshesUnion(Rpolynomials, sphereCenter, sphereRadius, angle_bound, radius_bound, distance_bound, error_bound));
    return rcpp_result_gen;
END_RCPP
}
// Mandelbulb
Rcpp::XPtr<EMesh3> Mandelbulb(int maxloop, Rcpp::NumericVector sphereCenter, double sphereRadius, double angle_bound, double radius_bound, double distance_bound, double error_bound);
RcppExport SEXP _cgalMeshes_Mandelbulb(SEXP maxloopSEXP, SEXP sphereCenterSEXP, SEXP sphereRadiusSEXP, SEXP angle_boundSEXP, SEXP radius_boundSEXP, SEXP distance_boundSEXP, SEXP error_boundSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type maxloop(maxloopSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type sphereCenter(sphereCenterSEXP);
    Rcpp::traits::input_parameter< double >::type sphereRadius(sphereRadiusSEXP);
    Rcpp::traits::input_parameter< double >::type angle_bound(angle_boundSEXP);
    Rcpp::traits::input_parameter< double >::type radius_bound(radius_boundSEXP);
    Rcpp::traits::input_parameter< double >::type distance_bound(distance_boundSEXP);
    Rcpp::traits::input_parameter< double >::type error_bound(error_boundSEXP);
    rcpp_result_gen = Rcpp::wrap(Mandelbulb(maxloop, sphereCenter, sphereRadius, angle_bound, radius_bound, distance_bound, error_bound));
    return rcpp_result_gen;
END_RCPP
}
// meshTopology
Rcpp::IntegerMatrix meshTopology(int nu, int nv, bool uperiodic, bool vperiodic);
RcppExport SEXP _cgalMeshes_meshTopology(SEXP nuSEXP, SEXP nvSEXP, SEXP uperiodicSEXP, SEXP vperiodicSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< int >::type nv(nvSEXP);
    Rcpp::traits::input_parameter< bool >::type uperiodic(uperiodicSEXP);
    Rcpp::traits::input_parameter< bool >::type vperiodic(vperiodicSEXP);
    rcpp_result_gen = Rcpp::wrap(meshTopology(nu, nv, uperiodic, vperiodic));
    return rcpp_result_gen;
END_RCPP
}
// sTriangle
Rcpp::XPtr<EMesh3> sTriangle(Rcpp::NumericVector A, Rcpp::NumericVector B, Rcpp::NumericVector C, Rcpp::NumericVector center, double radius, unsigned int iterations);
RcppExport SEXP _cgalMeshes_sTriangle(SEXP ASEXP, SEXP BSEXP, SEXP CSEXP, SEXP centerSEXP, SEXP radiusSEXP, SEXP iterationsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type A(ASEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type B(BSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type C(CSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type center(centerSEXP);
    Rcpp::traits::input_parameter< double >::type radius(radiusSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type iterations(iterationsSEXP);
    rcpp_result_gen = Rcpp::wrap(sTriangle(A, B, C, center, radius, iterations));
    return rcpp_result_gen;
END_RCPP
}

RcppExport SEXP _rcpp_module_boot_class_CGALmesh();

static const R_CallMethodDef CallEntries[] = {
    {"_cgalMeshes_AFSreconstruction_cpp", (DL_FUNC) &_cgalMeshes_AFSreconstruction_cpp, 2},
    {"_cgalMeshes_cxhull", (DL_FUNC) &_cgalMeshes_cxhull, 1},
    {"_cgalMeshes_cxhullsIntersection", (DL_FUNC) &_cgalMeshes_cxhullsIntersection, 2},
    {"_cgalMeshes_gTriangle", (DL_FUNC) &_cgalMeshes_gTriangle, 5},
    {"_cgalMeshes_AlgebraicMesh", (DL_FUNC) &_cgalMeshes_AlgebraicMesh, 9},
    {"_cgalMeshes_AlgebraicMeshesIntersection", (DL_FUNC) &_cgalMeshes_AlgebraicMeshesIntersection, 7},
    {"_cgalMeshes_AlgebraicMeshesUnion", (DL_FUNC) &_cgalMeshes_AlgebraicMeshesUnion, 7},
    {"_cgalMeshes_Mandelbulb", (DL_FUNC) &_cgalMeshes_Mandelbulb, 7},
    {"_cgalMeshes_meshTopology", (DL_FUNC) &_cgalMeshes_meshTopology, 4},
    {"_cgalMeshes_sTriangle", (DL_FUNC) &_cgalMeshes_sTriangle, 6},
    {"_rcpp_module_boot_class_CGALmesh", (DL_FUNC) &_rcpp_module_boot_class_CGALmesh, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_cgalMeshes(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
