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
Rcpp::XPtr<EMesh3> AFSreconstruction_cpp(const Rcpp::NumericMatrix pts);
RcppExport SEXP _cgalMeshes_AFSreconstruction_cpp(SEXP ptsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix >::type pts(ptsSEXP);
    rcpp_result_gen = Rcpp::wrap(AFSreconstruction_cpp(pts));
    return rcpp_result_gen;
END_RCPP
}
// AlgebraicMesh
Rcpp::XPtr<EMesh3> AlgebraicMesh(Rcpp::IntegerMatrix powers, Rcpp::NumericVector coeffs, double isolevel, Rcpp::NumericVector sphereCenter, double sphereRadius, double angle_bound, double radius_bound, double distance_bound);
RcppExport SEXP _cgalMeshes_AlgebraicMesh(SEXP powersSEXP, SEXP coeffsSEXP, SEXP isolevelSEXP, SEXP sphereCenterSEXP, SEXP sphereRadiusSEXP, SEXP angle_boundSEXP, SEXP radius_boundSEXP, SEXP distance_boundSEXP) {
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
    rcpp_result_gen = Rcpp::wrap(AlgebraicMesh(powers, coeffs, isolevel, sphereCenter, sphereRadius, angle_bound, radius_bound, distance_bound));
    return rcpp_result_gen;
END_RCPP
}
// mandelbulb
Rcpp::XPtr<EMesh3> mandelbulb(double angle_bound, double radius_bound, double distance_bound);
RcppExport SEXP _cgalMeshes_mandelbulb(SEXP angle_boundSEXP, SEXP radius_boundSEXP, SEXP distance_boundSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type angle_bound(angle_boundSEXP);
    Rcpp::traits::input_parameter< double >::type radius_bound(radius_boundSEXP);
    Rcpp::traits::input_parameter< double >::type distance_bound(distance_boundSEXP);
    rcpp_result_gen = Rcpp::wrap(mandelbulb(angle_bound, radius_bound, distance_bound));
    return rcpp_result_gen;
END_RCPP
}
// spikes
Rcpp::XPtr<EMesh3> spikes(double angle_bound, double radius_bound, double distance_bound);
RcppExport SEXP _cgalMeshes_spikes(SEXP angle_boundSEXP, SEXP radius_boundSEXP, SEXP distance_boundSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type angle_bound(angle_boundSEXP);
    Rcpp::traits::input_parameter< double >::type radius_bound(radius_boundSEXP);
    Rcpp::traits::input_parameter< double >::type distance_bound(distance_boundSEXP);
    rcpp_result_gen = Rcpp::wrap(spikes(angle_bound, radius_bound, distance_bound));
    return rcpp_result_gen;
END_RCPP
}
// VoxelToMesh
Rcpp::XPtr<EMesh3> VoxelToMesh(std::string filename, double isovalue, Rcpp::NumericVector center, double radius, double angle_bound, double radius_bound, double distance_bound);
RcppExport SEXP _cgalMeshes_VoxelToMesh(SEXP filenameSEXP, SEXP isovalueSEXP, SEXP centerSEXP, SEXP radiusSEXP, SEXP angle_boundSEXP, SEXP radius_boundSEXP, SEXP distance_boundSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type filename(filenameSEXP);
    Rcpp::traits::input_parameter< double >::type isovalue(isovalueSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type center(centerSEXP);
    Rcpp::traits::input_parameter< double >::type radius(radiusSEXP);
    Rcpp::traits::input_parameter< double >::type angle_bound(angle_boundSEXP);
    Rcpp::traits::input_parameter< double >::type radius_bound(radius_boundSEXP);
    Rcpp::traits::input_parameter< double >::type distance_bound(distance_boundSEXP);
    rcpp_result_gen = Rcpp::wrap(VoxelToMesh(filename, isovalue, center, radius, angle_bound, radius_bound, distance_bound));
    return rcpp_result_gen;
END_RCPP
}

RcppExport SEXP _rcpp_module_boot_class_CGALmesh();

static const R_CallMethodDef CallEntries[] = {
    {"_cgalMeshes_AFSreconstruction_cpp", (DL_FUNC) &_cgalMeshes_AFSreconstruction_cpp, 1},
    {"_cgalMeshes_AlgebraicMesh", (DL_FUNC) &_cgalMeshes_AlgebraicMesh, 8},
    {"_cgalMeshes_mandelbulb", (DL_FUNC) &_cgalMeshes_mandelbulb, 3},
    {"_cgalMeshes_spikes", (DL_FUNC) &_cgalMeshes_spikes, 3},
    {"_cgalMeshes_VoxelToMesh", (DL_FUNC) &_cgalMeshes_VoxelToMesh, 7},
    {"_rcpp_module_boot_class_CGALmesh", (DL_FUNC) &_rcpp_module_boot_class_CGALmesh, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_cgalMeshes(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
