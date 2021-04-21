#include <Rcpp.h>
using namespace Rcpp;
NumericVector reorder(NumericVector vec, NumericVector ord){ 
  int n(vec.size());
  NumericVector reord(n);
  for(int i=0;i<n;i++) reord[i]=vec[ord[i]-1];
  return reord;
}
NumericVector rowSums(NumericMatrix a){
  NumericVector ans(a.nrow());
  for(int i=0;i<a.nrow();i++){
    ans[i]=sum(a(i,_));
  }
  return ans;
}
NumericVector colSums(NumericMatrix a){
  NumericVector ans(a.ncol());
  for(int i=0;i<a.ncol();i++){
    ans[i]=sum(a(_,i));
  }
  return ans;
}
NumericMatrix EmisPr_NB(NumericVector X, NumericVector eSize, NumericVector esAvg) {
  NumericVector eProb=eSize/(eSize+esAvg);
  int nbT = X.size(), nbS = eProb.size();
  NumericMatrix res(nbS, nbT);
  
  for(int r=0; r<nbS; r++){
    for(int t=0; t<nbT; t++){
      res(r,t) =  R::dnbinom(X(t),eSize(r), eProb(r),false);
    }
  }
  return res;
}
double objFunc(NumericVector par,NumericVector X,NumericVector postPrv){
  double size=par[0];
  double mu=par[1];
  double prob=size/(size+mu);
  int nbT=X.size();
  NumericVector tmp(nbT);
  for(int i=0;i<nbT;i++){
    tmp[i]=R::dnbinom(X[i],size, prob,false);
  }
  return -sum(tmp * postPrv);
}
/*-------------------------------------------------------*
 *              --- function : mixNB_NB_NS ---
 // Input:
 // X, sample
 // nbS, number of subpopulation
 // avg, initial value of mean
 // var, initial value of variance
 // iterMax, maximum number of iterations
 // EPSILON, precision/error aloowed
 *-------------------------------------------------------*/
// [[Rcpp::export]]
List mixNB_NB_NS(NumericVector X, int nbS,
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
  NumericMatrix emisPr(nbS,nbT),rawPost(nbS,nbT),postPr(nbS,nbT);
  
  for(int iter=0;iter<iterMax;iter++){
    NumericVector prevParam;
    for(int i=0;i<esAvg.size();i++)prevParam.push_back(esAvg[i]);
    for(int i=0;i<eSize.size();i++)prevParam.push_back(eSize[i]);
    for(int i=0;i<eWeight.size();i++)prevParam.push_back(eWeight[i]);
    // cat("iter=",iter,"\n")
    //// Estep ————
    emisPr = EmisPr_NB(X, eSize, esAvg);
    
    // rawPost = emisPr*eWeight; 

    for(int r=0; r<nbS; r++){
      for(int t=0; t<nbT; t++){
        rawPost(r,t) = emisPr(r,t)*eWeight[r] ;
      }
    }
    NumericVector rawPost_colSums=colSums(rawPost);
    for(int r=0; r<nbS; r++){
      for(int t=0; t<nbT; t++){
        postPr(r,t) = rawPost(r,t)*rawPost_colSums[t] ;
      }
    }
    // postPr = t(t(rawPost)/colSums(rawPost);
    
    // return List::create(_["postPr"]=postPr);
    //// Mstep ————
    //// BFGS
    List res;
    Rcpp::Function optim("optim");
    for(int k=0;k<nbS;k++){
      NumericVector par={eSize[k], esAvg[k]};
      NumericVector postPrv=postPr(k,_);
      // objFunc <- function(par){-sum(dnbinom(X, size=par[1], mu=par[2], log=TRUE) * postPr[k,])}
      res=optim(Rcpp::_["par"]    = par,
                Rcpp::_["fn"]     = Rcpp::InternalFunction(&objFunc),
                Rcpp::_["method"] = "BFGS",
                Rcpp::_["X"]=X,
                Rcpp::_["postPrv"]=postPrv);
      // eSize[k] = eParam[1]; esAvg[k] = eParam[2]
      NumericVector newpar=res["par"];
      eSize[k] = newpar[0]; esAvg[k] = newpar[1];

    }

    NumericVector sumWeight = rowSums(postPr); eWeight = sumWeight/sum(sumWeight);

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
  return List::create(_["parameter"]=parameter,_["emisPr"]=emisPr,_["postPr"]=postPr);
  
}
/*** R
res<-mixNB_NB_NS(X,nbS)
  */