// #include <Rcpp.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

double g(const double& llambda, const NumericVector& time, const int& pos1, const int& pos2) {
  double res = 0;

  int max2 = 0;
  if(pos2 == time.size() - 1) max2 = 1;

  res += std::exp(llambda * (time[pos1] - time[pos2]));

  if(max2 == 1 && pos1 == 0) {
    // if(debug == 1) Rcout<<"first pos1, last pos2";
    return res;
  }
  if(max2 == 1 && pos1 != 0) {
    // if(debug == 1) Rcout<<"last pos2";
    res += (-1) * std::exp(llambda * (time[pos1 - 1] - time[pos2]));
    return res;
  }
  if(max2 == 0 && pos1 == 0) {
    // if(debug == 1) Rcout<<"first pos1";
    res += (-1) * std::exp(llambda * (time[pos1] - time[pos2 + 1]));
    return res;
  }
  if(max2 == 0 && pos1 != 0) {
    res += (-1) * std::exp(llambda * (time[pos1 - 1] - time[pos2])) +
      (-1) * std::exp(llambda * (time[pos1] - time[pos2 + 1])) +
      std::exp(llambda * (time[pos1 - 1] - time[pos2 + 1]));
    return res;

  }

  return 0;


}

double logfactorial(const int &k) {
  double res = 0.0;
  if(k==0 || k==1) return res; else
    for(int i = 2; i<= k; i++)
      res += std::log(i);
  return res;
}

// psi functions
double psi_gamma(const double& ggamma, const double& c, const int& nderiv) {
  if(nderiv == 0) {
    return (std::log(ggamma + c) - std::log(ggamma));
  } else
    return pow(ggamma + c, -nderiv) * std::pow(std::exp(1.0), logfactorial(nderiv - 1)) *
      pow(-1.0, nderiv - 1);
}

double psi_stab(const double& ggamma,
                     const double& c,
                     const int& nderiv) {
  if(nderiv == 0) {
    return (std::pow(c , ggamma));
  } else

    return pow(c, ggamma - nderiv) *
      ggamma * std::exp(lgamma(nderiv - ggamma) - lgamma(1-ggamma)) *
      pow(-1.0, nderiv + 1);
}

double psi_pvf(const double& ggamma, const double &pvfm,
                    const double& c, const int& nderiv) {

  double sign = 1.0;
  if(pvfm < 0) sign = -1.0;

  if(nderiv == 0) {
    return (std::pow(ggamma / ggamma, pvfm) - std::pow(ggamma / (ggamma + c), pvfm)) * sign;
    //return alpha * (1 - std::pow( bbeta / (bbeta + c), pvfm)) * sign;
  } else  {
    // Rcout<<pow(bbeta, pvfm)<<" * "<<std::pow(bbeta + c, -pvfm - nderiv)<<" * "<<std::exp(lgamma(pvfm + nderiv) - lgamma(pvfm));
    // Rcout<<" = "<<alpha * std::pow(bbeta, pvfm) * std::pow(bbeta + c, -pvfm - nderiv) * std::exp(lgamma(pvfm + nderiv) - lgamma(pvfm))<<std::endl;
    // Rcout<<"sign = "<<pow(-1.0, nderiv) * sign;
    //return alpha * std::pow(bbeta, pvfm) * std::pow(bbeta + c, -pvfm - nderiv) * std::exp(lgamma(pvfm + nderiv) - lgamma(pvfm)) *  pow(-1.0, nderiv + 1);// * sign;
    return std::pow(ggamma, pvfm) * std::pow(ggamma + c, -pvfm - nderiv) *
      std::exp(lgamma(pvfm + nderiv) - lgamma(pvfm)) *  pow(-1.0, nderiv + 1);
  }

}


// dist: 0 gamma, 1 stable, 2 pvf
double psi(const int &dist, const double &pvfm, const double &ggamma, const int &deriv, const double &ck1k2) {
  if(dist == 0) return psi_gamma(ggamma, ck1k2, deriv); else
    if(dist == 1) return psi_stab(ggamma, ck1k2, deriv); else
      return psi_pvf(ggamma, pvfm, ck1k2, deriv);
}

// logic: S takes all sums of cvec on positions containing {left right}
double S(const int& deriv, const int& left, const int& right, const double& ggamma,
         const int& dist, const double& pvfm,
         const NumericVector& times, const double& llambda, const NumericVector& cvec) {

  double res = 0.0;

  for(unsigned int i = 0; i < left; i++)
    for(unsigned int j = right - 1; j < cvec.size(); j++) {
      // Rcout<<"psi^("<<deriv<<")(sum[c_"<<i<<"^"<<j<<"])";
      double ck1k2 = std::accumulate(cvec.begin() + i, cvec.begin() + j + 1, 0.0);
      res += psi(dist, pvfm, ggamma, deriv, ck1k2) * g(llambda, times, i, j);
    }

    // Rcout<<std::endl;
    return res;
}



double divideSum(const std::vector<int> &events, const NumericVector &cvec, const double &aalpha, const double &ggamma,
               const int &dist, const double &pvfm, const NumericVector &times, const double &llambda) {

  int n = events.size();
  if(n==0) return 1.0;
  double res = 0.0;


  // get the machinery to list all the partitions
  std::vector<int> kappa(n, 0);
  std::vector<int> M(n,0);

  // First partition is here: 0 0 0 (all in the same set)

  res += (-1.0) * aalpha * S(n, events[0], events[n-1],
   ggamma, dist, pvfm, times, llambda, cvec);

  // Rcout<<"res is"<<res<< "  ////";

  int i = n-1;

  while(i!=0) {

    if(kappa[i] <= M[i-1]) {

      ++kappa[i];

      const int new_max = std::max(M[i], kappa[i]);
      M[i] = new_max;

      for(int j=i+1; j<n; ++j) {
        kappa[j] = 0;
        M[j] = new_max;
      }

      // int kappa_size = M[n - 1] + 1;

      // This is a new partition now (kappa).

      // Rcout<<"The partition: ";
      // for(std::vector<int>::iterator it2 = kappa.begin(); it2!=kappa.end(); it2++) {
      //   Rcout<<*it2<<" ";
      // }
      // Rcout<<std::endl;

      // Get arguments for S

      std::vector<std::vector<int> > indexes(M[n - 1] + 1, std::vector<int>(3,0));

      for(int j = 0; j != kappa.size(); ++j) {

        //the length of interval goes here (for the detemrination of the function.)

        int nset = kappa[j];

        indexes[nset][0]++;

        if(indexes[nset][1]==0 || indexes[nset][1] > events[j]) indexes[nset][1] = events[j];

        if(indexes[nset][2] < events[j]) indexes[nset][2] = events[j];

      }

      // for(int p = 0; p < indexes.size(); p++)
      //   Rcout<<"("<<indexes[p][0]<<","<<indexes[p][1]<<","<<indexes[p][2]<<")  ";

      // Rcout<<std::endl;

      double res_part = 1.0;

      for(int j = 0; j < indexes.size(); j++) {
        res_part = res_part * (-1.0) * aalpha * S(indexes[j][0], indexes[j][1], indexes[j][2],
          ggamma, dist, pvfm, times, llambda, cvec);
      }


      res += res_part;

      if(i+1<n) i = n-1;


    } else i--;
  }

return res;
}

// [[Rcpp::export]]
NumericVector Estep_id(NumericVector events, NumericVector cvec, double aalpha, double ggamma,
                int dist, double pvfm, NumericVector times, double llambda) {

  std::vector<int> events_add(events.begin(), events.end());
  double denom = divideSum(events_add, cvec, aalpha, ggamma, dist, pvfm, times, llambda);


  // Rcout<<std::endl<<"--- now the numerator"<<std::endl;
  std::vector<double> num;

  std::vector<int> base = events_add;
  // Rcout<<"base: ";
  // for(std::vector<int>::iterator it = base.begin(); it != base.end(); it++)
  //   Rcout<<*it<<" ";

  //Rcout<<std::endl;
  // For each time point...
  for(int i = 1; i<= cvec.size(); i++) {

    // Create the new vector...
    std::vector<int> newvec = events_add;

    std::vector<int>::iterator pos = std::upper_bound(newvec.begin(), newvec.end(), i);

    newvec.insert(pos, i);
    // for(std::vector<int>::iterator it = newvec.begin(); it != newvec.end(); it++)
    //   Rcout<<*it<<" ";
    //
    // Rcout<<" // "<<divideSum(newvec, cvec, aalpha, ggamma, dist, pvfm, times, llambda);
    num.push_back(divideSum(newvec, cvec, aalpha, ggamma, dist, pvfm, times, llambda));

    // Rcout<<std::endl;
  }

  NumericVector out(num.begin(), num.end());

  out.push_back(denom);


  // Also add the log-Laplace transform
  double logLaplace = 0.0;

  // Rcout<<"cvec size: "<<cvec.size()<<" times size "<<times.size();
  // Rcout<<std::endl;

  for(unsigned int i = 0; i < cvec.size(); i++)
    for(unsigned int j = i; j < cvec.size(); j++) {

      double ck1k2 = std::accumulate(cvec.begin() + i, cvec.begin() + j + 1, 0.0);
      logLaplace += psi(dist, pvfm, ggamma, 0, ck1k2) * g(llambda, times, i, j);
      // Rcout<<"("<<i + 1<<" "<<j + 1<<")="<<g(llambda, times, i, j)<<" / ";

      // Rcout<<p;
    }

  logLaplace *= (-1) * aalpha;

  out.push_back(logLaplace);

  return out;

}

void Vcov_adj_id(const NumericVector &events,
                  const NumericVector &cvec,
                  const double &aalpha,
                  const double &ggamma,
                  const int &dist,
                  const double &pvfm,
                  const NumericVector &times,
                  const double &llambda,
                  const NumericVector& elp,
                  const arma::mat &xelph,
                  const arma::mat &tau,
                  const Rcpp::List &interval_rows,
                  const NumericVector &ez,
                  const int &n_times,
                  arma::mat &betabeta, // this will get modified
                  std::vector<double> &lambdalambda, // this will get modified
                  arma::mat &betalambda // this will get modified
                 ) {

  // make the List into a vector of vectors
  std::vector<std::vector<int> > int_rows(interval_rows.size());

  for(int i = 0; i<interval_rows.size(); i++) {
    Rcpp::NumericVector s1(interval_rows[i]);
    int_rows[i] = std::vector<int>(s1.begin(), s1.end());
  }

  std::vector<int> events_add(events.begin(), events.end());

  double denom = divideSum(events_add, cvec, aalpha, ggamma, dist, pvfm, times, llambda);

  // For each interval
  for(int i = 1; i<= cvec.size(); i++) {


    // the vector of events; this is from 1 to L
    std::vector<int> newvec = events_add;

    std::vector<int>::iterator pos = std::upper_bound(newvec.begin(), newvec.end(), i);

    newvec.insert(pos, i);

    for(int j = 1; j<= cvec.size(); j++) {

      std::vector<int> newvec2 = newvec;

      std::vector<int>::iterator pos2 = std::upper_bound(newvec2.begin(), newvec2.end(), j);

      newvec2.insert(pos2, j);


      // Rcout<<"Adding intervals "<<i<<" and "<<j<<", ";
      // Rcout<<"numerator: ";
      // for(std::vector<int>::iterator it = newvec2.begin(); it != newvec2.end(); it++)
      //   Rcout<<*it<<" ";
      // Rcout<<std::endl;


      //e[z1z2]
      double exy_exey = divideSum(newvec2, cvec, aalpha, ggamma, dist, pvfm, times, llambda)/denom - ez[i-1] * ez[j-1];

      // Rcout<<"exey: "<<exy_exey;
      // Rcout<<"i="<<i<<", j="<<j<<std::endl;

      // std::vector<double> tmp1(xelph.n_cols, 0.0);
      // std::vector<double> tmp2(xelph.n_cols, 0.0);
      arma::rowvec tmp_i(xelph.n_cols, arma::fill::zeros);
      arma::rowvec tmp_j(xelph.n_cols, arma::fill::zeros);

      // l1 is the rows of interval i
      for(std::vector<int>::iterator l1 = int_rows[i-1].begin(); l1 != int_rows[i-1].end(); l1++) {

        tmp_i += xelph.row(*l1 - 1);
        // Rcout<<"add line "<<*l1<<" to tmp1";
        //
        // Rcout<<std::endl;

        // l2 is the rows of interval j
        for(std::vector<int>::iterator l2 = int_rows[j-1].begin(); l2 != int_rows[j-1].end(); l2++) {

          // Rcout<<std::endl;
          // Rcout<<"Looking at lines: "<<*l1<<" "<<*l2<<" /// ";
          // Rcout<<"["<<tau(*l1-1, 0)<<", "<<tau(*l1-1, 1)<<") // ";
          // Rcout<<"["<<tau(*l2-1, 0)<<", "<<tau(*l2-1, 1)<<") //";
          // Rcout<<std::endl;

          // p is the times within the line l1

          for(int p = tau(*l1 - 1, 0); p < tau(*l1 - 1, 1); p++)  {


            betalambda.row(p) += exy_exey * elp[*l1-1] * xelph.row(*l2-1);

            // when i!=j this would normally get added twice
            // if(i == j) {
            //   betalambda.row(p) += exy_exey * elp[*l1-1] * xelph.row(*l2-1);
            // } else
            //   betalambda.row(p) += exy_exey * elp[*l1-1] * xelph.row(*l2-1) +
            //     exy_exey * xelph.row(*l1 -1) * elp[*l2 - 1];

            // q is the times within line l2

            for(int q = tau(*l2 -1, 0); q < tau(*l2 - 1, 1); q++)

              if(p <= q) {
                lambdalambda[n_times * p - p * (p+1)/2 + q] += exy_exey * elp[*l1-1] * elp[*l2-1];

                // if (i==j) {
                //   lambdalambda[n_times * p - p * (p+1)/2 + q] += exy_exey * elp[*l1-1] * elp[*l2-1];
                // } else {
                //   lambdalambda[n_times * p - p * (p+1)/2 + q] += exy_exey * 2.0 * elp[*l1-1] * elp[*l2-1];
                // }
              }


          }

          // this must be added just once; I loop through these values more times unfortunately
          if(l1 == int_rows[i-1].begin()) {
            tmp_j += xelph.row(*l2 - 1);
            // Rcout<<"add line "<<*l2<<" to tmp2";
            // Rcout<<std::endl;
          }



        }


      }

      betabeta += exy_exey * trans(tmp_i) * tmp_j;

        // if(i == j) {
        //   betabeta += exy_exey * trans(tmp1) * tmp2;
        // } else
        //   betabeta += exy_exey * trans(tmp1) * tmp2 +
        //     exy_exey * trans(tmp2) * trans(tmp1);

        // point being: because we don't loop through ALL the combination (
        //  j >= i


    }

  }
}


// [[Rcpp::export]]
List Vcov_adj(List events_l,
              List cvec_l,
              double aalpha,
              double ggamma, int dist, double pvfm,
              List times_l,
              double llambda,
              List elp_l,
              List xelph_l,
              List tau_l,
              List interval_rows_l,
              List ez_l,
              int n_times,
              int n_covs) {

  // outcome

  arma::mat betabeta(n_covs, n_covs, arma::fill::zeros);
  std::vector<double> lambdalambda(n_times * (n_times+1)/2, 0.0);
  arma::mat betalambda(n_times, n_covs, arma::fill::zeros);

    for(int indiv = 0; indiv < events_l.size(); indiv++)
      Vcov_adj_id(events_l[indiv],
                  cvec_l[indiv],
                  aalpha, ggamma, dist, pvfm,
                  times_l[indiv], llambda,
                  elp_l[indiv],
                  xelph_l[indiv],
                  tau_l[indiv],
                  interval_rows_l[indiv],
                  ez_l[indiv],
                  n_times,
                  betabeta,
                  lambdalambda,
                  betalambda);


    return Rcpp::List::create(Rcpp::Named("betabeta") = betabeta,
                              Rcpp::Named("betalambda") = betalambda,
                              Rcpp::Named("lambdalambda") = lambdalambda);

}
/*** R
set.seed(1)
# Last one is the denominator
# First ones are numerators
# Also need the Laplace transform
# Estep_id(events = c(2,3), cvec = rexp(6), aalpha = 2, ggamma = 2,
#            dist = 0, pvfm = -1/2, times = 1:6)
*/
