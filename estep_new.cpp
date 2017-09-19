#include <Rcpp.h>
using namespace Rcpp;

double g(const double& llambda, NumericVector& time, const int& pos1, const int& pos2) {
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
double S(int& deriv, int& left, int& right, double& ggamma,
         int& dist, double& pvfm,
         NumericVector& times, double& llambda, NumericVector& cvec) {

  double res = 0.0;

  for(unsigned int i = 0; i < left; i++)
    for(unsigned int j = right - 1; j < cvec.size(); j++) {
      // Rcout<<"("<<i + 1<<" "<<j + 1<<")";
      double ck1k2 = std::accumulate(cvec.begin() + i, cvec.begin() + j + 1, 0.0);
      res += psi(dist, pvfm, ggamma, deriv, ck1k2) * g(llambda, times, i, j);
    }

    return res;
}



double divideSum(std::vector<int> &events, NumericVector &cvec, double &aalpha, double &ggamma,
               int &dist, double &pvfm, NumericVector &times, double &llambda) {

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

        //Rcout<<" adding to set number... "<<nset<<" pertaining to time point "<<timepoints[i]<<std::endl;
        indexes[nset][0]++;

        if(indexes[nset][1]==0 || indexes[nset][1] > events[j]) indexes[nset][1] = events[j];

        if(indexes[nset][2] < events[j]) indexes[nset][2] = events[j];

      }

      // for(int i = 0; i < indexes.size(); i++)
      //   Rcout<<"("<<indexes[i][0]<<","<<indexes[i][1]<<","<<indexes[i][2]<<")  ";
      //
      // Rcout<<std::endl;

      //


      // double S(int deriv, int left, int right,
      //          double ggamma,
      //          double lambda,
      //          int dist, double pvfm,
      //          NumericVector times, NumericVector cvec)

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

  std::vector<double> num;

  // For each time point...
  for(int i = 1; i<= cvec.size(); i++) {

    // Create the new vector...
    std::vector<int> newvec = events_add;

    std::vector<int>::iterator pos = std::upper_bound(newvec.begin(), newvec.end(), i);

    newvec.insert(pos, i);
    // for(std::vector<int>::iterator it = newvec.begin(); it != newvec.end(); it++)
    //   Rcout<<*it<<"  ";
    num.push_back(divideSum(newvec, cvec, aalpha, ggamma, dist, pvfm, times, llambda));

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

/*** R
set.seed(1)
# Last one is the denominator
# First ones are numerators
# Also need the Laplace transform
# Estep_id(events = c(2,3), cvec = rexp(6), aalpha = 2, ggamma = 2,
#            dist = 0, pvfm = -1/2, times = 1:6)
*/
