// #include <RcppArmadillo.h>
// #include <RcppEigen.h>
// 
// // [[Rcpp::depends(RcppArmadillo)]]
// // [[Rcpp::depends(RcppEigen)]]
// 
// // [[Rcpp::export]]
// SEXP eigenMatMult(Eigen::MatrixXd A, Eigen::MatrixXd B){
//     Eigen::MatrixXd C = A * B;
//     
//     return Rcpp::wrap(C);
// }
// 
// // [[Rcpp::export]]
// SEXP eigenMapMatMult(const Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> B){
//     Eigen::MatrixXd C = A * B;
//     
//     return Rcpp::wrap(C);
// }
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List rcpp_hello_world() {
    
    CharacterVector x = CharacterVector::create( "foo", "bar" )  ;
    NumericVector y   = NumericVector::create( 0.0, 1.0 ) ;
    List z            = List::create( x, y ) ;
    
    return z ;
}
