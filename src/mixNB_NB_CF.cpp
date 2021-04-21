#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]
NumericVector reorder(NumericVector vec, NumericVector ord){ 
  int n(vec.size());
  NumericVector reord(n);
  for(int i=0;i<n;i++) reord[i]=vec[ord[i]-1];
  return reord;
}
arma::vec vec2arma(Rcpp::NumericVector x) {
  return arma::vec(x.begin(), x.size(),false);// this code originate from P159 Seamless Rcpp
}
template <typename T>
Rcpp::NumericVector arma2vec(const T& x) {
  return Rcpp::NumericVector(x.begin(), x.end());
}
mat EmisPr(mat X,mat N, mat M){
  mat res(X.n_rows,X.n_cols);
  for(unsigned int i=0;i<X.n_rows;i++){
    for(unsigned int j=0;j<X.n_cols;j++){
      double m=M(i,j),x=X(i,j),n=N(i,j);
      double p=1.0*n/(n+m);
      double tmp=R::dnbinom(x,n,p,false);
      res(i,j)=tmp;
        
    }
  }
  return res;
}
mat digamma(mat X){
  mat res(X.n_rows,X.n_cols);
  for(unsigned int i=0;i<X.n_rows;i++){
    for(unsigned int j=0;j<X.n_cols;j++){
      double tmp=R::digamma(X(i,j));
      res(i,j)=tmp;
      
    }
  }
  return res;
}
/*-------------------------------------------------------*
*              --- function : mixNB_NB_CF ---
// Input:
// X, sample
// nbS, number of subpopulation
// avg, initial value of mean
// var, initial value of variance
// iterMax, maximum number of iterations
// EPSILON, precision/error aloowed
*-------------------------------------------------------*/
// [[Rcpp::export]]
List mixNB_NB_CF(NumericVector X, int nbS,
                 int iterMax=15000, double EPSILON=1e-7){
  NumericVector eWeight  (nbS, 1.0/nbS); // the initial value of weight
  int nbT = X.size(); // amount of samples
  NumericVector esAvg;
  NumericVector esVar;
  
  Rcpp::Function kmeans("kmeans");
  Rcpp::Function order("order");
  List res_kmeans = kmeans(X, nbS);
  
  esAvg={1};
  esVar={0};
  while(is_true(all(esAvg>esVar))){
    
    NumericVector esAvg_init(nbS);
    NumericVector esVar_init(nbS);
    NumericVector centers=res_kmeans[1];
    NumericVector withinss=res_kmeans[3];
    NumericVector size=res_kmeans[6];
    for(int i=0;i<nbS;i++)
    {
      esAvg_init[i] = centers[i];
      esVar_init[i] = withinss[i]/ (size[i]-1) ;
    }
    esAvg = clone(esAvg_init);
    esVar = clone(esVar_init);
    NumericVector ord = order(esVar);
    esVar=reorder(esVar,ord);
    esAvg=reorder(esAvg,ord);
  }
  
  NumericVector eSize = pow(esAvg,2)/(esVar-esAvg); // calculate initial value of size
  NumericVector eProb = esAvg/esVar; // calculate initial value of probability
  
  mat ones_nbS = ones(nbS,1); // a vector with nbS elements equals 1
  mat ones_nbT = ones(nbT,1); // a vector with nbT elements equals 1
  mat ones_nbTS = ones(nbT,nbS); // a matrix with nbT*nbS elements equals 1
  mat ones_nbST = ones(nbS,nbT); // a matrix with nbS*nbT elements equals 1
  mat ones_nbSS = ones(nbS,nbS); // a matrix with nbS*nbS elements equals 1
  mat Pmat,postPr;
  //iteration starts
  for(int iter=0;iter<iterMax;iter++){
  // for(int iter=0;iter<1;iter++){
    // mat prevParam = join_cols(vec2arma(esAvg), vec2arma(eSize),vec2arma(eWeight)); // define iterative parameters in an arrray
    NumericVector prevParam;
    for(int i=0;i<esAvg.size();i++)prevParam.push_back(esAvg[i]);
    for(int i=0;i<eSize.size();i++)prevParam.push_back(eSize[i]);
    for(int i=0;i<eWeight.size();i++)prevParam.push_back(eWeight[i]);
    mat Nmat = vec2arma(eSize)*ones_nbT.t(); // matrix of size
    mat Mmat = vec2arma(esAvg)*ones_nbT.t(); // matrix of mean
    mat Xmat = ones_nbS*vec2arma(X).t();     // matrix of sample X
    // Estep --------
     Pmat = EmisPr(Xmat,Nmat, Mmat); // matrix of emission Probability
    mat rawPost = Pmat%(vec2arma(eWeight)*ones_nbT.t());  postPr = rawPost/(ones_nbSS*rawPost);// matrix of post Probability
    NumericVector beta = 1 - 1/(1-eProb) - 1/log(eProb); // value of beta
    mat Bmat = vec2arma(beta)*ones_nbT.t(); // matrix of beta
    mat Dmat = Nmat%(digamma(Nmat+Xmat)-digamma(Nmat)); //matrix of delta

    //// Mstep ---------
    mat tmp1=(postPr*Dmat.t()),tmp2=(postPr*(Xmat-((ones_nbST-Bmat)%Dmat)).t()),tmp3=postPr*ones_nbTS;
    mat eProb_mat = vec2arma(beta)%(tmp1.diag())/(tmp2.diag()); // new value of probability
    eProb=arma2vec(eProb_mat);
    mat eSize_mat = -tmp1.diag()/tmp3.diag()/vec2arma(log(eProb)); // new value of size
    eSize=arma2vec(eSize_mat);
    vec sumWeight = sum(postPr,1);

    mat eWeight_mat=sumWeight/sum(sumWeight); // new value of weight
    eWeight=arma2vec(eWeight_mat);
    esAvg = eSize * (1-eProb)/eProb;  // new value of mean
    esVar = eSize * (1-eProb)/pow(eProb,2);  // new value of variance

    Rcpp::Function order("order");
    NumericVector ordAvg = order(esAvg); // order of mean
    esAvg = reorder(esAvg,ordAvg); // rearrange mean with increasing order
    eSize = reorder(eSize,ordAvg); // rearrange size in increasing order of mean
    eWeight= reorder(eWeight,ordAvg); // rearrange weight in increasing order of mean
    eProb = reorder(eProb,ordAvg); // rearrange probability in increasing order of mean
    esVar = reorder(esVar,ordAvg); // rearrange variance in increasing order of mean

    // mat newParam = join_cols(vec2arma(esAvg), vec2arma(eSize),vec2arma(eWeight)); // update iterative parameters in an arrray
    NumericVector newParam;
    for(int i=0;i<esAvg.size();i++)newParam.push_back(esAvg[i]);
    for(int i=0;i<eSize.size();i++)newParam.push_back(eSize[i]);
    for(int i=0;i<eWeight.size();i++)newParam.push_back(eWeight[i]);
      // throw warning when any element in the iterative parameters can't be calculated
      if(any(is_na(esAvg)) | any(is_na(esAvg))){
        warning("mixNB_NB_CF - NA is appeared!");
        NumericVector a;
        a=seq_len(nbS);
        esAvg=prevParam[a-1];
        eSize=prevParam[a+nbS-1];
        break;
      }
      double D = sum(pow((newParam - prevParam),2)); // when precision reach the desired level, stop iteration, otherwise reach the maximum limit
      if(D < EPSILON){break;};
    }

    // return with a list
    List parameter; parameter["eSize"]=eSize; parameter["eProb"]=eProb; parameter["esAvg"]=esAvg;
    parameter["esVar"]=esVar; parameter["eWeight"]=eWeight;
    return List::create(_["parameter"]=parameter,_["emisPr"]=Pmat,_["postPr"]=postPr);

  }
  
  /*** R
  res<-mixNB_NB_CF(X,nbS)
    */