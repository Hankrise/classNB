#include <Rcpp.h>
using namespace Rcpp;
/*-------------------------------------------------------*
 *              --- function : LogSumExp ---
 * output: sum of a vector in log scale 
 *-------------------------------------------------------*/
// [[Rcpp::export]]
double LogSumExp(NumericVector logV) {
  bool posInflogV = is_true(any(logV==R_PosInf));
  bool negInflogV = is_true(all(logV==R_NegInf));
  double res = 0.0;
  // NumericVector logVM = logV[logV>R_NegInf];
  // double maxlogVM = max(logVM);
  double maxlogV = max(logV);
  if(posInflogV){
    stop("existing positive inifinite value in logVector \n");
  }else if(negInflogV){
    res = R_NegInf;
  }else{
    // res = maxlogVM+ log(sum(exp(logVM-maxlogVM)));
    res = maxlogV + log(sum(exp(logV-maxlogV)));
  }
  return res;
}


/*-------------------------------------------------------*
 *              --- function : C_sum2Mat ---
 * output: A+B
 *-------------------------------------------------------*/
// [[Rcpp::export]]
NumericMatrix C_sum2Mat(NumericMatrix A, NumericMatrix B){
  NumericMatrix res(A.nrow(),A.ncol());
  for(int i=0; i<A.nrow(); i++){
    for(int j=0; j<A.ncol(); j++){
      res(i,j) = A(i,j) + B(i,j);
    }
  }
  return res;
}

/*-------------------------------------------------------*
 *              --- function : ParamVec ---
 * to gather initPr, transPr, eParam into a vector
 *-------------------------------------------------------*/
// [[Rcpp::export]]
SEXP ParamVec(SEXP A_, SEXP B_, SEXP C_, SEXP D_){
  NumericVector A(A_), B(B_), C=(C_), D=(D_);
  int nbA = A.size(), nbB = B.size(), nbC = C.size(), nbD = D.size();
  NumericVector res = no_init_vector(nbA+nbB+nbC+nbD);
  res[seq_len(nbA)-1] = A;
  for(int i=0; i<nbB; i++){ res[nbA+i] = B[i]; }
  for(int i=0; i<nbC; i++){ res[nbA+nbB+i] = C[i]; }
  for(int i=0; i<nbD; i++){ res[nbA+nbB+nbC+i] = D[i]; }
  return res;
}

/*-------------------------------------------------------*
 *              --- function : OrdParam ---
 * sort the parameters by the order of mean
 *-------------------------------------------------------*/
// [[Rcpp::export]]
List OrdParam(NumericVector initPr, NumericMatrix transPr, NumericVector esAvg, NumericVector eSize, NumericVector eProb){
  NumericVector initPr_=clone(initPr), esAvg_=clone(esAvg), eSize_=clone(eSize), eProb_=clone(eProb); 
  NumericMatrix transPr_=clone(transPr);
  Function ord("order");
  IntegerVector ordParam = ord(esAvg_);
  
  for(int q=0; q<initPr_.size(); q++){
    initPr_[q] = initPr[ordParam[q]-1];
    eSize_[q] = eSize[ordParam[q]-1];
    eProb_[q] = eProb[ordParam[q]-1];
    transPr_(q,_) = transPr(ordParam[q]-1,_);
  }
  NumericMatrix transTmp = clone(transPr_);
  for(int q=0; q<initPr_.size(); q++){
    transPr_(_,q) = transTmp(_,ordParam[q]-1);
  }
  return List::create(_["initPr"]=initPr_, _["transPr"]=transPr_, _["esAvg"]=esAvg_, _["eSize"]=eSize_, _["eProb"]=eProb_);
}

/*-------------------------------------------------------*
 *              --- function : FwdBcd---
 * return forward probability and normalized constant
 * return posterior probabilities
 *-------------------------------------------------------*/
// [[Rcpp::export]]
List ForwardBackward_HMM(NumericVector initPr, NumericMatrix transPr, NumericMatrix emisPr){
  short unsigned int nbT = emisPr.ncol(), nbS = emisPr.nrow(); 
  // std::cout<<"init="<<initPr<<"\n";
  // std::cout<<"emisPr=\n"<<emisPr<<"\n";
  /* >>>>>>>>> Forward <<<<<<<<<*/
  NumericMatrix fPr(nbS, nbT); NumericVector vecCst(nbT);
  // for t=1, 
  vecCst[0] = sum(initPr*emisPr(_,0));
  fPr(_,0) = (initPr*emisPr(_,0))/vecCst[0];
  // for t=2,...,T
  for(int t=1; t<nbT; t++){
    NumericVector vecTmp(nbS);
    for(int l=0; l<nbS; l++){
      vecTmp[l] = emisPr(l,t)*sum(fPr(_,t-1)*transPr(_,l));//exp(LogSumExp(log(fPr(_,t-1))+log(transPr(_,l))));//
    }
    vecCst[t] = sum(vecTmp);
    fPr(_,t) = vecTmp/vecCst[t];
  }
  double logNorm = sum(log(vecCst));
  
  /* >>>>>>>>> Backward <<<<<<<<<*/
  NumericMatrix postPr(nbS,nbT), gPr(nbS,nbT);
  //for T
  postPr(_,nbT-1) = fPr(_,nbT-1);
  //for T-1, ..., 1
  for(int t=nbT-2; t>-1; t--){
    for(int l=0; l<nbS; l++){
      // gPr(l,t+1) = exp(LogSumExp(log(fPr(_,t))+log(transPr(_,l))));
      gPr(l,t+1) = sum(fPr(_,t)*transPr(_,l));
    }
    for(int k=0; k<nbS; k++){
      // postPr(k,t) = fPr(k,t)*exp(LogSumExp(log(transPr(k,_))+log(postPr(_,t+1))-log(gPr(_,t+1))));
      postPr(k,t) = fPr(k,t)* sum(transPr(k,_)*postPr(_,t+1)/gPr(_,t+1) );
    }
  }
  return List::create(_["fPr"]=fPr,_["logNorm"]=logNorm,_["postPr"]=postPr,_["gPr"]=gPr);
}

/*----------- function : EmisPr_NB --------------*
 * output: Negative binomial emission probability matrix 
 *---------------------------------------------------*/
// [[Rcpp::export]]
NumericMatrix EmisPr_NB(NumericVector X, NumericVector eSize, NumericVector eProb) {
  int nbT = X.size(), nbS = eProb.size();
  NumericMatrix res(nbS, nbT);
  
  for(int r=0; r<nbS; r++){
    for(int t=0; t<nbT; t++){
      res(r,t) =  R::dnbinom(X(t),eSize(r), eProb(r),false);
      //::Rf_dlnorm(X(t),eParam(r,0), sqrt(eParam(r,1)),false);//
      //exp(-pow(X(t)-eParam(r,0),2.0)/(2*eParam(r,1)))/sqrt(2*PI*eParam(r,1));//PI
    }
  }
  return res;
}

/*----------- function : EmisPr_mixNB --------------*
 * output: mixture negative binomial emission probability matrix 
 *---------------------------------------------------*/
// [[Rcpp::export]]
NumericMatrix EmisPr_mixNB(NumericVector X, IntegerVector component, NumericVector eWeight, 
                           NumericVector eSize, NumericVector eProb) {
  unsigned short int nbT = X.size(), nbK = component.size();
  unsigned short int index=0;
  NumericMatrix res(nbK, nbT);
  
  for(int k=0; k<nbK; k++){
    for(int t=0; t<nbT; t++){
      for(int c=index; c<index + component(k); c++){// problem for start index
        res(k,t) +=  exp(log(eWeight(c))+ R::dnbinom(X(t),eSize(c), eProb(c),true));
      }
    }
    index += component(k);
  }
  return res;
}

/*-------------------------------------------------------*
 *              --- function : EM_NB ---
 * Learning the parameters of HMM with NB distribution
 * return posterior probabilities and parameters
 *-------------------------------------------------------*/
// [[Rcpp::export]]
List EM_NB(NumericVector X, NumericVector initPr, NumericMatrix transPr,
           NumericVector eSize, NumericVector eProb, int iterMax=500, double EPSILON=1.0E-7){
  
  short unsigned int nbT(X.size()), nbS(initPr.size());
  NumericMatrix delta(nbS,nbT);
  NumericVector beta(nbS), lambda(nbS), esAvg(nbS);
  
  NumericMatrix emisPr(nbS, nbT), postPr(nbS,nbT);
  NumericVector param_new = ParamVec(initPr, transPr, eSize, eProb);
  NumericVector param_old = ParamVec(initPr, transPr, eSize, eProb);
  
  int iterStop=iterMax; double logLik=0.0;
  int convergence=0;
  for(int iter=0; iter<iterMax; iter++){
    NumericVector prev_initPr = clone(initPr), prev_eSize=clone(eSize), prev_eProb=clone(eProb), prev_esAvg=clone(esAvg);
    NumericMatrix prev_transPr = clone(transPr);
    // E-step -------------------------------------------------------------------
    // postPr by Forward-Backward
    emisPr = EmisPr_NB(X, eSize, eProb);
    List resFB = ForwardBackward_HMM(initPr, transPr, emisPr);
    NumericMatrix fPr = resFB["fPr"], gPr = resFB["gPr"];
    postPr = as<NumericMatrix>(resFB["postPr"]);
    logLik = resFB["logNorm"];
    // delta
    for(int t=0; t<nbT; t++){
      for(int q=0; q<nbS; q++){
        delta(q,t) = eSize[q]*(Rf_digamma(eSize[q]+X[t]) - Rf_digamma(eSize[q])); 
      }
    }
    // beta
    beta = 1-1/(1-eProb)-1/log(eProb);
    // M-step -------------------------------------------------------------------
    // initPr & transPr
    for(int k=0; k<nbS; k++){
      //update transPr
      NumericMatrix transTmp(nbS,nbS);
      for(int l=0; l<nbS; l++){
        double sumTmp=0.0;
        for(int t=1; t<nbT; t++){
          sumTmp += fPr(k,t-1) * postPr(l,t)/gPr(l,t);
        }
        transTmp(k,l) = transPr(k,l)*sumTmp;
      }
      transPr(k,_) = transTmp(k,_)/sum(transTmp(k,_));
      //update initPr
      initPr(k) = postPr(k,0);
    }
    // lambda
    NumericVector numeratorLambda(nbS), denominatorLambda(nbS);
    for(int k=0; k<nbS; k++){
      numeratorLambda[k] = sum(postPr(k,_)*delta(k,_));
      denominatorLambda[k] = sum(postPr(k,_));
    }
    lambda = numeratorLambda/denominatorLambda;
    // eProb
    NumericVector numeratorProb(nbS), denominatorProb(nbS);
    for(int k=0; k<nbS; k++){
      numeratorProb[k] = beta[k]*sum(postPr(k,_)*delta(k,_));
      denominatorProb[k] = sum(postPr(k,_)*(X-(1-beta[k])*delta(k,_)));
    }
    eProb = numeratorProb/denominatorProb;
    esAvg = eSize*(1-eProb)/eProb;
    eSize = -lambda/log(eProb);
    
    if(any(is_na(eSize)) | any(is_na(eProb))){
      Rf_warning("Convergence error, NA is appeared!");
      initPr = prev_initPr; transPr = prev_transPr;
      esAvg = prev_esAvg; eProb = prev_eProb; eSize = prev_eSize;
      convergence=1;
      break;
    }
    
    List resOrd = OrdParam(initPr, transPr, esAvg, eSize, eProb);
    initPr  = as<NumericVector>(resOrd["initPr"]); 
    transPr = as<NumericMatrix>(resOrd["transPr"]); 
    esAvg   = as<NumericVector>(resOrd["esAvg"]);
    // esAvg = NumericVector::create(25,75,125,175,250,350);
    eSize   = as<NumericVector>(resOrd["eSize"]); 
    eProb   = as<NumericVector>(resOrd["eProb"]);
    
    param_new = ParamVec(initPr, transPr, eSize, eProb);
    double D  = sum(pow(abs(param_new-param_old),2.0));
    param_old = param_new;
    if(D<EPSILON){iterStop = iter+1;  break;}
    std::cout<<iter<<"\t";
  }//end for
  NumericVector esVar = eSize*(1-eProb)/pow(eProb,2);
  
  List parameter; parameter["initPr"]=initPr; parameter["transPr"]=transPr;
  parameter["eProb"]=eProb; parameter["eSize"]=eSize; parameter["esAvg"]=esAvg;  parameter["esVar"]=esVar; 
  return List::create(_["parameter"]=parameter,_["loglik"]=logLik,_["postPr"]=postPr,_["emisPr"]=emisPr,_["convergence"]=convergence);
}  
/*-------------------------------------------------------*
 *              --- function : Viterbi ---
 * Viterbi algorithm to find optimal path 
 *-------------------------------------------------------*/
// [[Rcpp::export]]
SEXP Viterbi_HMM(SEXP initPr_, SEXP transPr_, SEXP emisPr_){
  NumericVector initPr(initPr_);
  NumericMatrix transPr(transPr_),  emisPr(emisPr_);
  int nbT = emisPr.ncol(), nbK = emisPr.nrow();
  NumericMatrix logF(nbK,nbT);
  
  logF(_,0) = log(initPr) + log(emisPr(_,0));
  
  for(int t=1; t<nbT; t++){
    for(int r =0; r<nbK; r++){
      logF(r,t) = log(emisPr(r,t)) + max(logF(_,t-1) + log(transPr(_,r)));
    }
  }
  
  IntegerVector path(nbT);
  path[nbT-1] = which_max(logF(_,nbT-1));
  for(int t=nbT-2; t>-1; t--){
    path[t] = which_max(logF(_,t) + log(transPr(_,path[t+1])));
  }
  path = path+1;
  return path; 
}


/*-------------------------------------------------------*
 *              --- function : InitEM ---
 * initialize the parameters in HMM
 *-------------------------------------------------------*/
// [[Rcpp::export]]
List InitEM(NumericVector X, int nbS){
  //initialize probability
  NumericVector initPr(nbS, 1.0/nbS);
  //transition probability
  NumericMatrix matTmp(nbS,nbS);std::fill(matTmp.begin(),matTmp.end(),0.5/nbS);
  NumericMatrix transPr=C_sum2Mat(NumericMatrix::diag(nbS,0.5),matTmp);
  //emission parameters
  Function      kmeans("kmeans");
  List          clust = kmeans(X,nbS);
  NumericVector mm = clust[1]; std::sort(mm.begin(),mm.end());
  double totwithinss = clust[4]; 
  // double totvar =  totwithinss / X.size();
  NumericVector esAvg(mm); NumericVector esVar=rep(totwithinss,nbS);
  NumericVector eSize=pow(esAvg,2)/(esVar-esAvg); NumericVector eProb=esAvg/esVar;
  
  return List::create(_["initPr"]=initPr, _["transPr"]=transPr, _["eSize"]=eSize, _["eProb"]=eProb);
}

/*-------------------------------------------------------*
 *              --- function : HMM_NB_CF ---
 * HMM to search the optimal path : closed form
 *-------------------------------------------------------*/
// [[Rcpp::export]]
List HMM_NB_CF(NumericVector X, int nbS, bool viterbi=true, Nullable<NumericVector> initPr=R_NilValue, Nullable<NumericMatrix> transPr=R_NilValue,
                       Nullable<NumericVector> eSize=R_NilValue, Nullable<NumericVector> eProb=R_NilValue,
                       int iterMax=500, double EPSILON=1.0E-7){
  
  NumericVector initPrIn(nbS), eSizeIn(nbS), eProbIn(nbS); NumericMatrix transPrIn(nbS,nbS);
  List resEM;
  if(initPr.isNotNull() && transPr.isNotNull() && eSize.isNotNull() && eProb.isNotNull()){
    initPrIn=as<NumericVector>(clone(initPr)); transPrIn=as<NumericMatrix>(clone(transPr)); 
    eSizeIn=as<NumericVector>(clone(eSize)); eProbIn=as<NumericVector>(clone(eProb));  
  }else{
    List initEM = InitEM(X, nbS);
    initPrIn=as<NumericVector>(initEM["initPr"]); transPrIn=as<NumericMatrix>(initEM["transPr"]); 
    eSizeIn=as<NumericVector>(initEM["eSize"]); eProbIn=as<NumericVector>(initEM["eProb"]);
  }

  resEM = as<List>(EM_NB(X, initPrIn, transPrIn, eSizeIn, eProbIn, iterMax, EPSILON));
  
  List paramOut = as<List>(resEM["parameter"]);
  NumericMatrix  postPrOut = as<NumericMatrix>(resEM["postPr"]);
  NumericMatrix  emisPrOut = as<NumericMatrix>(resEM["emisPr"]);
  double loglikOut  = as<double>(resEM["loglik"]);
  bool convergOut   = as<bool>(resEM["convergence"]);
  
  NumericVector initPrOut = paramOut["initPr"];
  NumericMatrix transPrOut = paramOut["transPr"];
    
  int nbT = X.size(); IntegerVector path(nbT);
  if(viterbi){
    path = Viterbi_HMM(initPrOut, transPrOut, emisPrOut);
  }else{
    for(int i=0; i<nbT; i++){path[i]=which_max(postPrOut(_,i))+1;}
  }  

  return List::create(_["parameter"]=paramOut, _["postPr"]=postPrOut, _["logLik"]=loglikOut, 
                      _["convergence"]=convergOut, _["status"]=path); 
}  
















