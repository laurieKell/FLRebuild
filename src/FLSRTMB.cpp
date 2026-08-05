// Use TMB_LIB_INIT to register TMB symbols (getParameterOrder, etc.)
#define TMB_LIB_INIT R_init_FLRebuild

#include <TMB.hpp>
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

// Forward declaration for Rcpp registration
extern void R_init_FLRebuild_Rcpp(DllInfo *dll);

// Copyright Henning Winker (JRC) & Iago MOSQUEIRA (WMR), 2021
// Authors:  Henning Winker (JRC) <henning.winker@ec.europa.eu>
// Iago MOSQUEIRA (WMR) <iago.mosqueira@wur.nl>


// Space time
template<class Type>

// objective function
Type objective_function<Type>::operator() () {

  // Data
  DATA_VECTOR( rec );
  DATA_VECTOR( ssb );
  DATA_VECTOR( prior_s ); // Prior vector for s, [logit(mean), stdev in logit, useflag]
  DATA_VECTOR(prior_r0); // Optional prior for log_r0: [mean, sd]
  DATA_VECTOR( spr0 );
  DATA_INTEGER(use_b0); // 0 = v(t)=r0*spr0(t); 1 = v(t)=b0(t)
  DATA_VECTOR( b0yr );
  DATA_INTEGER(nyears);
  DATA_INTEGER(Rmodel); // Recruitment model

  // Parameters
  PARAMETER(log_r0);
  PARAMETER(log_sigR);
  PARAMETER(logit_s);
  PARAMETER(log_inflect);
  // Derived quantities
  Type r0 = exp(log_r0);
  Type sigR = exp(log_sigR);
  Type s = 0.2001 + 0.7998*1/(1+exp(-logit_s));
  vector<Type> log_rec_hat(nyears);
  Type a = 0;
  Type b = 0;

  // Objective function
  Type ans=0;
  vector<Type> v(nyears);

  for(int t=0; t<nyears; t++){
    if(use_b0 == 1)
      v(t) = b0yr(t);
    else
      v(t) = r0 * spr0(t);
  }

  if(Rmodel==0){ // bevholtSV()
   for( int t=0; t< nyears; t++){
     log_rec_hat(t) = log(4.0 * s * r0 *ssb(t) / (v(t)*(1.0-s)+ssb(t)*(5.0*s-1.0)));//-pow(sigR,2)/2.0;
     }}

   if(Rmodel==1){ // rickerSV()
     for( int t=0; t< nyears; t++){
       b = log(5.0*s)/(0.8*v(t));
       if(use_b0 == 1)
         a = r0 * exp(b*v(t))/v(t);
       else
         a = exp(b*v(t))/spr0(t);
       log_rec_hat(t) = log(a*ssb(t)*exp(-b*ssb(t)));
     }}
   
   if(Rmodel==2){ // segreg() aka Hockey Stick
     for( int t=0; t< nyears; t++){
       log_rec_hat(t) = log(r0)+log(2.5*s/v(t)*(ssb(t)+exp(log_inflect)-pow(pow(ssb(t)-exp(log_inflect),2.0),0.5)));//-pow(sigR,2)/2.0;
     }}


   // Depensatory Beverton and Holt
   if(Rmodel==3){ 
     for( int t=0; t< nyears; t++){
       log_rec_hat(t) = log(r0)+log(2.5*s/v(t)*(ssb(t)+0.2*v(t)/s-pow(pow(ssb(t)-0.2*v(t)/s,2.0),0.5)));//-pow(sigR,2)/2.0;
     }}
   
   vector<Type> rec_hat = exp(log_rec_hat);
 
   // LL
   for( int t=0; t<nyears; t++){
     ans -= dnorm( log(rec(t)), log_rec_hat(t), sigR, true );
   }
 
  // loglAR1 function(obs, hat, rho = 0) 
 
   //prior s
   ans -= dnorm(logit_s, prior_s(0), prior_s(1), 1); // Prior for logn
  
  // prior for log_r0 if provided
  if(prior_r0.size() == 2) {
    ans -= dnorm(log_r0, prior_r0(0), prior_r0(1), true);
  }
//   
//   if(Rmodel==0){
//   a = Type(4)*v*s/(spr0*(Type(5)*s-Type(1)));
//   b = v*(Type(1)-s)/(Type(5)*s-Type(1));
//   }
//   
//   if(Rmodel==1){
//     //b = log(5.0*s)/(0.8*v);
//     //a = exp(b*v)/spr0;
//   }
//   
//   if(Rmodel==2){
//     b = 0.2*v/s;
//     a = r0/b;
//   }
//   
//   // prior s
//   
//   
   // Reporting
   REPORT( rec_hat );
   REPORT( nyears );
   REPORT( sigR );
   REPORT( r0 );
   REPORT( v );
   REPORT( log_inflect );
//   REPORT( a );
//   REPORT( b );
   REPORT( s );
   ADREPORT(r0);
   if(use_b0 == 0) ADREPORT(v);
// 
   return ans;}
