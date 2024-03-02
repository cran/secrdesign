/*
   External procedures for secrdesign package 2022-10-21
 
   2023-09-15 better protection from divide-by-zero in Lambdacpp

*/

#include "secrdesign.h"

/*==============================================================================*/

// [[Rcpp::export]]
arma::mat hazmatcpp (
        const arma::vec &par, 
        const arma::mat &d, 
        const int &detectfn) {
    
    arma::mat H = d;
    
    int kk = d.n_rows;    /* number of traps */
    int mm = d.n_cols;    /* number of points on mask */
    int k,m;

    //-------------------------------------------------------------------------
    
    if (detectfn == 0 || detectfn == 14) {       // HHN
        H = par(0) * arma::exp(-arma::square(d) / 2 / par(1) / par(1));
    }
    else if (detectfn == 1 || detectfn == 15) {  // HHR
        H = par(0) * (1 - arma::exp(- arma::pow(d /par(1), - par(2))));
    } 
    else if (detectfn == 2 || detectfn == 16) {  // HEX
        H = par(0) * arma::exp(-d / par(1));
    } 
    else if (detectfn == 6 || detectfn == 17) {  // HAN
        H = par(0) * arma::exp(-arma::square(d-par(2)) / 2 / par(1) / par(1));
    } 
    else if (detectfn == 8 || detectfn == 18) {  // HCG
        for (k=0; k<kk; k++) {
            for (m=0; m<mm; m++) {
                boost::math::gamma_distribution<> gam(par(2), par(1)/par(2));
                H(k,m) = par(0) * boost::math::cdf(complement(gam,d(k,m)));
            }
        }
    } 
    else if (detectfn == 19) {                   // HVP
        H = par(0) * arma::exp(- arma::pow(d / par(1), par(2)));
    }
    else Rcpp::stop ("detectfn not implemented");
    
    return (H);
}
/*==============================================================================*/

// [[Rcpp::export]]
Rcpp::List Lambdacpp (
        const int &type, 
        const arma::vec &par, 
        const arma::mat &d, 
        const int &detectfn)
{
    // traps x mask hazard matrix
    arma::mat h = hazmatcpp(par, d, detectfn);
    
    // column sums
    arma::rowvec sumhk = arma::sum(h, 0); 
    // cast as Rcpp NumericVector
    Rcpp::NumericVector outsumhk = Rcpp::NumericVector(sumhk.begin(), sumhk.end());

    Rcpp::NumericVector outsumpk (1);
    Rcpp::NumericVector outsumq2 (1);
    outsumpk[0] = NA_REAL;
    outsumq2[0] = NA_REAL;
    
    // column sums of squared elements
    arma::rowvec sumhk2 = arma::sum(arma::square(h), 0); 
    
    // protect against 1/0 - revised 2023-09-15
    arma::rowvec sqsumhk = arma::square(sumhk);
    sqsumhk.replace(0, arma::datum::eps);
    arma::rowvec sumq2 = sumhk2 / sqsumhk;
    // Rprintf("%8.6g sum(sumq2)\n", arma::accu(sumq2));

    outsumq2 = Rcpp::NumericVector(sumq2.begin(), sumq2.end());
    
    // multi
    if (type == 0) {
        arma::rowvec sumpk = 1 - arma::exp(- sumhk);
        outsumpk = Rcpp::NumericVector(sumpk.begin(), sumpk.end());
    }
    // proximity
    else if (type == 1) {
        arma::rowvec sumpk = sum(1 - arma::exp(- h), 0);
        outsumpk = Rcpp::NumericVector(sumpk.begin(), sumpk.end());
    }
    // count: sumpk not used
    
    return Rcpp::List::create(
        Rcpp::Named("sumhk") = outsumhk,
        Rcpp::Named("sumpk") = outsumpk,
        Rcpp::Named("sumq2") = outsumq2
    );
}
/*==============================================================================*/

// [[Rcpp::export]]
Rcpp::List Qpmcpp (
        const arma::vec &par,
        const arma::rowvec &D,
        const arma::mat &d, 
        const int &detectfn,
        const int &noccasions)
{
    double Qp, Qpm;
    double G = arma::accu(D);
    arma::mat pij;
    arma::rowvec pd, p0, p1, p2;

    // traps x mask matrix hazard per occasion
    arma::mat h = hazmatcpp(par, d, detectfn);
    
    // mask probability detection
    pd = 1 - arma::exp(- arma::sum(h,0) * noccasions);
    p0 = 1-pd;
    pd = pd % D;
    Qp = arma::accu(pd) / G;

    pij = 1 - arma::exp(- h * noccasions);
    pij = pij/(1-pij);
    p1 = p0 % arma::sum(pij, 0);   // sum over traps
    p2 = D % (1 - p0 - p1);
    Qpm = arma::accu(p2) / G;

    return (Rcpp::List::create(
            Rcpp::Named("Qp") = Qp,
            Rcpp::Named("Qpm") = Qpm));
}
/*==============================================================================*/

// [[Rcpp::export]]
Rcpp::List En2cpp (
        const int &type, 
        const arma::vec &par,
        const arma::rowvec &D,
        const arma::mat &d, 
        const int &detectfn,
        const int &noccasions)
{
    // suffix k refers to detectors, m to mask cells
    arma::mat hkm = hazmatcpp(par, d, detectfn);
    
    arma::mat pkm = d;
    arma::rowvec p0, p1, D1, D2;
    arma::uword i;
    
    arma::rowvec Hm = sum(hkm,0);       // hazard summed over traps
    Hm.replace(0, arma::datum::eps);    // protect divide by zero
    
    p0 = arma::exp(-Hm * noccasions);   // Pr not detected
    
    double En, En2;                     // return values
    
    // multi-catch traps, by occasion
    if (type == 0) {
        arma::mat Hkm = hkm;
        for (i = 0; i < hkm.n_rows; i++ ) {
            Hkm.row(i) = Hm;
        }
        arma::mat Pkm = (1 - arma::exp(-Hkm));
        pkm = hkm/Hkm % Pkm;
        pkm = (1 - arma::pow(1 - pkm, noccasions)) % 
            arma::pow(1 - Pkm + pkm, noccasions-1);
        p1  = arma::sum(pkm,0);        // trapped at only one site
    }
    
    // binary and count proximity detectors, all occasions
    else {
        pkm = 1 - arma::exp(-hkm * noccasions);
        pkm = pkm/(1-pkm);              // fails if any pkm==1?
        p1  = p0 % arma::sum(pkm,0);    // detected at only one site       
    }
    
    D1 = D % (1 - p0);
    D2 = D % (1 - p0 - p1);
    En  = arma::accu(D1);
    En2 = arma::accu(D2);
    
    return (Rcpp::List::create(
            Rcpp::Named("En") = En,
            Rcpp::Named("En2") = En2));
    
}
/*==============================================================================*/
