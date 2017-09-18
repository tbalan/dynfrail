//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
//#include <Rcpp.h>
using namespace Rcpp; using namespace arma;

#include <iostream>
#include <vector>



/*
repr - generates a partition 'structure' of size n. Each number is a coding for the partition, e.g.
for a set of size 4, 01203 means the partition (1,4)(2)(3), and a total of 3 subsets. The last element of the returned vector is the number of sets.

voidrepr - similar with rep, does not return the vectors, just calculates them.

printvecvec - prints vectors<vectors>

part_to_index : takes a partition code and makes it into arguments for an S function.

*/
void printvecvec(std::vector<std::vector<int> > vec ) {
  std::vector<std::vector<int> >::iterator it;

  for(it = vec.begin(); it < vec.end(); it++) {

    std::vector<int> minivec = *it;
    std::vector<int>::iterator it2;

    for(it2 = minivec.begin(); it2 < minivec.end(); it2++) {
      Rcout<<*it2<<" ";

    }
    Rcout<<std::endl;
  }
}

void printvec(std::vector<int> v) {
  Rcout<<std::endl;
  Rcout<<"( ";
  for(std::vector<int>::iterator it = v.begin(); it != v.end(); it++)
    Rcout<<*it<<" ";
  Rcout<<")";
}

unsigned int factorial(unsigned int n)
{
  unsigned int retval = 1;
  for (int i = n; i > 1; --i)
    retval *= i;
  return retval;
}

double logfactorial(const int &k) {
  double res = 0.0;
  if(k==0 || k==1) return res; else
    for(int i = 2; i<= k; i++)
      res += std::log(i);
  return res;
}


// Generates a representation object.
std::vector<std::vector < int> > repr(int n) {

  std::vector<std::vector<int> > res;

  std::vector<int> kappa(n, 0);
  std::vector<int> M(n,0);

  // add the first partition representation & length
  res.push_back(kappa);
  res[0].push_back(1);


  int i=n-1;

  // int done=0;

  while(i!=0) {



    if(kappa[i] <= M[i-1]) {

      ++kappa[i];

      const int new_max = std::max(M[i], kappa[i]);
      M[i] = new_max;

      for(int j=i+1; j<n; ++j) {
        kappa[j] = 0;
        M[j] = new_max;
      }


      std::vector<int> kappa_size = kappa;
      kappa_size.push_back(M[n - 1] + 1);

      res.push_back(kappa_size);


      if(i+1<n) i = n-1;


    } else i--;
  }

  return res;

}




// timepoints is a vector of timepoints corresponding to the set being partitioned.
std::vector<std::vector<int> > part_to_index(std::vector<int>& partition, const std::vector<int>& timepoints) {

  std::vector<std::vector<int> > indexes(partition.back(), std::vector<int>(3,0));
  partition.pop_back();

  for(int i = 0; i != partition.size(); ++i) {
    //*it tells us the number of the set that we are at.

    //the length of interval goes here (for the detemrination of the function.)


    int nset = partition[i];

    //Rcout<<" adding to set number... "<<nset<<" pertaining to time point "<<timepoints[i]<<std::endl;
    indexes[nset][0]++;

    // Assume that the std::vector timepoints is already ordered ascending, it it contains no 0! (only from 1 on)
    if(indexes[nset][1]==0 || indexes[nset][1] > timepoints[i]) indexes[nset][1]  = timepoints[i];

    // the maximal time point, always updated.
    if(indexes[nset][2] < timepoints[i]) indexes[nset][2] = timepoints[i];


  }

  return indexes;

}


// [[Rcpp::export]]
double S(const int& len, const int& posstart, const int& posstop,  const mat& A,  const mat& H, const double& bbeta, const int& nev) {

  //Rf_PrintValue(A);
  //Rcout<<std::endl<<"======"<<std::endl;

  // NumericMatrix subA = A(Range(0, posstart-1), Range(posstop-1, nev-1));
  // NumericMatrix subH = H(Range(0, posstart-1), Range(posstop-1, nev-1));

  //Rf_PrintValue(subA);
  //Rf_PrintValue(subH);
  //  factorial(m-1)* A / (H+bbeta)^m

  //X.submat( first_row, first_col, last_row, last_col )

  // double res = sum(factorial(len-1) * subA / pow(subH + bbeta, len)) ;

  //uvec XXX(Range(0, posstart - 1));
  //mat Asub = A((Range(0, posstart-1)), as<uvec>(Range(posstop-1, nev-1)));

  // Works:
  // double res = accu(factorial(len-1) * A.submat(0, posstop-1, posstart -1, nev-1) /
  //                   pow(A.submat(0, posstop-1, posstart -1, nev-1) + bbeta, len)) ;

  // Rcout<<"S calculating: "<<len - 1<<" factorial times A[0:,"<<posstart - 1<<"]["<<posstop - 1<<", "<<nev - 1<<"]"<<std::endl;
  // for(int i = posstart - 1; i<=posstop - 1; i++) {
  //   for( int j = posstop - 1; j < nev - 1; j++) {
  //     Rcout<<"e: "<<A[i,j];
  //   }
  // }
  double res = accu(factorial(len-1) * A.submat(arma::span(0, posstart - 1), arma::span(posstop-1, nev-1)) /
                    (pow(H.submat(arma::span(0, posstart - 1), arma::span(posstop-1, nev-1)) + bbeta, len)) ) ;

  // Rcout<<"S returning "<<res<<std::endl;
  return res;

}

// This function takes the vector of time points and gives the result. with derivative at only those points.
double vec_to_repr_to_number(const std::vector<int>&  tp_g, const mat& A, const mat& H, const double& bbeta, const int& nev) {

  Rcout<<"vec2repr2number"<<std::endl;
  double res = 0.0;

  int n = tp_g.size();


  std::vector<int> kappa(n, 0);
  std::vector<int> M(n,0);


  Rcout<<"nevents = "<<n<<"; start at = "<<tp_g[0]<<"; end at = "<<tp_g.back()<<"total events = "<<nev;


  Rcout<<std::endl<<"First contribution:"<<S(n, tp_g[0], tp_g.back(), A, H, bbeta, nev)<<std::endl;
  res += S(n, tp_g[0], tp_g.back(), A, H, bbeta, nev);



  int i=n-1;

  Rcout<<std::endl<<"start looping"<<std::endl;


  // start looping through all the possible partitions

  while(i!=0) {

    if(kappa[i] <= M[i-1]) {

      ++kappa[i];

      const int new_max = std::max(M[i], kappa[i]);
      M[i] = new_max;

      for(int j=i+1; j<n; ++j) {
        kappa[j] = 0;
        M[j] = new_max;
      }


      std::vector<int> kappa_size = kappa;
      kappa_size.push_back(M[n - 1] + 1);

      // This is a new partition now (kappa) with the last number the number of sets.

      //
      //
      // Rcout<<"The partition: ";
      // for(std::vector<int>::iterator it2 = kappa_size.begin(); it2!=kappa_size.end(); it2++) {
      //   Rcout<<*it2<<" ";
      // }
      // Rcout<<std::endl;
      //

      // From this particular partition we puke out an argument for the S function.

      std::vector<std::vector<int> > x = part_to_index(kappa_size, tp_g);


      // Basically for example for (12)(3)(4) I will have S2(12)
      double temp = 1;

      for(std::vector<std::vector<int> >::iterator it = x.begin(); it!=x.end(); it++) {
        std::vector<int> position = *it;
        temp = temp * S(position[0], position[1], position[2], A, H, bbeta, nev);
      }
      // Rcout<<temp<<std::endl;
      res+= temp;

      if(i+1<n) i = n-1;


    } else i--;
  }

  return res;
}




// This function takes an existing representation structure, a vector, and pukes out the number.
double repr_to_number(std::vector<std::vector <int> >  partitions, const std::vector<int>& tps, const mat& A, const mat& H, const double& bbeta, const int& nev) {

  double res = 0;

  //Rcout<<"repr_to_number here. I received as partitions:";
  //printvecvec(partitions);

  //Rcout<<std::endl<<"and as timepoints: ";
  //for(std::vector<int>::iterator it = tps.begin(); it<tps.end(); it++) Rcout<<" "<<*it;

  for(std::vector<std::vector<int> >::iterator itpart = partitions.begin(); itpart!=partitions.end(); itpart++) {

    // every element is a partition
    std::vector<int> part = *itpart;

    //Rcout<<std::endl<<"Selected partition: ";
    //for(std::vector<int>::iterator pl = part.begin(); pl<part.end(); pl++) Rcout<<*pl<<" ";

    // that we make an index for... or a bnch of indexes
    std::vector<std::vector<int> > x = part_to_index(*itpart, tps);

    //Rcout<<" The corresponding indexes are:"<<std::endl;
    //printvecvec(x);

    // Rcout<<std::endl<<"=========";
    double temp = 1;

    // from the indexes we compute the function.
    for(std::vector<std::vector<int> >::iterator it = x.begin(); it!=x.end(); it++) {
      std::vector<int> position = *it;
      temp = temp * S(position[0], position[1], position[2], A, H, bbeta, nev);
    }

    res+= temp;

  }


  return res;
}




// [[Rcpp::export]]
NumericVector SdivideC(const NumericVector& tmp, const NumericVector& all_tmp, const mat& A, const mat& H, const double& bbeta, const int& nev) {

  // Denominator:
  std::vector<int> tp_g(tmp.begin(), tmp.end());

  double denom = vec_to_repr_to_number(tp_g, A,  H, bbeta,  nev);

  Rcout<<endl<<"denominator: "<<denom;
  //Rcout<<"Done with the denominator."<<std::endl;

  // Generate a structure for the size that we need.
  std::vector<std::vector<int> > partitions_ext = repr(tmp.size() + 1);

  // Introduce sequentially all_tmp in the proper place in the tmp vector.
  std::vector<int> tp_all(all_tmp.begin(), all_tmp.end());
  std::vector<double> num;



  // For each time point...
  for(std::vector<int>::iterator it = tp_all.begin(); it < tp_all.end(); ++it) {

    std::vector<int> newvec = tp_g;

    std::vector<int>::iterator pos  = std::upper_bound(newvec.begin(), newvec.end(), *it);

    newvec.insert(pos, *it);

    // Calculate the sum of S from the representation structure already made
    num.push_back(repr_to_number(partitions_ext, newvec, A, H, bbeta, nev));

  }

  NumericVector out(num.begin(), num.end());

  out.push_back(denom);
  return out;
}

// [[Rcpp::export]]
List SdivideC_V(NumericVector tmp, NumericVector all_tmp, mat A, mat H, double bbeta, int nev) {

  // Denominator:
  std::vector<int> tp_g(tmp.begin(), tmp.end());

  double denom = vec_to_repr_to_number(tp_g, A,  H, bbeta,  nev);

  //Rcout<<"Done with the denominator."<<std::endl;

  // Generate a structure for the size that we need.
  std::vector<std::vector<int> > partitions_ext = repr(tmp.size() + 1);
  std::vector<std::vector<int> > partitions_ext_v = repr(tmp.size() + 2);

  // Introduce sequentially all_tmp in the proper place in the tmp vector.
  std::vector<int> tp_all(all_tmp.begin(), all_tmp.end());
  std::vector<double> num;
  std::vector<double> v;


  // For each time point...

  for(std::vector<int>::iterator it = tp_all.begin(); it < tp_all.end(); ++it) {

    std::vector<int> newvec = tp_g;

    // Put in the current time point at the right place....
    std::vector<int>::iterator pos  = std::upper_bound(newvec.begin(), newvec.end(), *it);

    newvec.insert(pos, *it);

    // Calculate the sum of S from the representation structure already made
    num.push_back(repr_to_number(partitions_ext, newvec, A, H, bbeta, nev));

    // Calculation of the V part

    for(std::vector<int>::iterator it2 = it; it2 < tp_all.end(); ++it2) {
      std::vector<int> newvec2 = newvec;
      std::vector<int>::iterator pos2  = std::upper_bound(newvec2.begin(), newvec2.end(), *it2);
      newvec2.insert(pos2, *it);
      v.push_back(repr_to_number(partitions_ext_v, newvec2, A, H, bbeta, nev));
    }


  }


  NumericVector out(num.begin(), num.end());
  out.push_back(denom);

  //return out;


  List ret = List::create(wrap(out), wrap(v));

  return ret;
}



// [[Rcpp::export]]
NumericVector SdivideC_iter(const NumericVector& tmp,
                            const NumericVector& all_tmp,
                            const mat& A, const mat& H,
                            const double& bbeta,
                            const int& nev) {

  // Denominator:
  std::vector<int> tp_g(tmp.begin(), tmp.end());
  //arma::mat matA = as<mat>(A);
  //arma::mat matH = as<mat>(H);

  double denom = vec_to_repr_to_number(tp_g, A,  H, bbeta,  nev);

  std::vector<int> tp_all(all_tmp.begin(), all_tmp.end());

  std::vector<double> num;

  // For each time point...
  for(std::vector<int>::iterator it = tp_all.begin(); it < tp_all.end(); ++it) {

    // Create the new vector...
    std::vector<int> newvec = tp_g;

    std::vector<int>::iterator pos  = std::upper_bound(newvec.begin(), newvec.end(), *it);

    newvec.insert(pos, *it);

    // Calculate the representation, give the sum of S.
    num.push_back(vec_to_repr_to_number(newvec, A,  H, bbeta,  nev));

  }
  NumericVector out(num.begin(), num.end());

  out.push_back(denom);
  return out;
}
