#include <Rcpp.h>
using namespace Rcpp;

NumericVector rowSums(NumericMatrix a){
  NumericVector ans(a.nrow());
  for(int i=0;i<a.nrow();i++){
    ans[i]=sum(a(i,_));
  }
  return ans;
}
NumericVector reorder(NumericVector vec, NumericVector ord){ 
  int n(vec.size());
  NumericVector reord(n);
  for(int i=0;i<n;i++) reord[i]=vec[ord[i]-1];
  return reord;
}
List initHMM(NumericVector X, int nbS, String meth_init="kmeans"){
  // NumericVector init (nbS,1/nbS*1.0);
  NumericVector init (nbS,1.0/nbS);
  NumericVector transv (nbS*nbS,0.5/(nbS));
  NumericMatrix trans (nbS,nbS,transv.begin());
  trans.fill_diag(0.5+0.5/(nbS));
  NumericVector esAvg;
  NumericVector esVar;
  NumericMatrix eParam(nbS,2);
  if(meth_init == "kmeans") {
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

      eParam(_,0)=esAvg;
      eParam(_,1)=esVar;

  }
  return List::create(_["init"]=init,_["trans"]=trans,_["eParam"]=eParam);
}
NumericMatrix EmisPr(NumericVector X, NumericMatrix eParam){
  // eParam : a matrix of size K*2 storing the mean and variance
  int nbS   = eParam.nrow();
  int nbT=X.size();
  NumericVector esAvg=eParam(_,0);
  NumericVector esVar=eParam(_,1);
  NumericVector eSize = pow(esAvg,2)/(esVar-esAvg);
  NumericVector eProb = esAvg/esVar;
  NumericMatrix res(nbS, nbT);

  for(int r=0; r<nbS; r++){
    for(int t=0; t<nbT; t++){
      res(r,t) =  R::dnbinom(X(t),eSize(r), eProb(r),false);
    }
  }
  return res;
}
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
      gPr(l,t+1) = sum(fPr(_,t)*transPr(_,l));
    }
    for(int k=0; k<nbS; k++){
      postPr(k,t) = fPr(k,t)* sum(transPr(k,_)*postPr(_,t+1)/gPr(_,t+1) );
    }
  }
  return List::create(_["fPr"]=fPr,_["logNorm"]=logNorm,_["postPr"]=postPr,_["gPr"]=gPr);
}

List Estep(NumericVector init, NumericMatrix trans, NumericMatrix emisPr,NumericMatrix eParam){
    List resFB  = ForwardBackward_HMM(init, trans, emisPr);
    NumericMatrix Fpr   = resFB["fPr"];
    double logLik   = resFB["logNorm"];
    NumericMatrix postPr   = resFB["postPr"];
    NumericMatrix Gpr   = resFB["gPr"];
    return List::create(_["Fpr"]=Fpr,_["Gpr"]=Gpr,_["postPr"]=postPr,_["logLik"]=logLik,_["eParam"]=eParam);
}

// objFunc <- function(par){-sum(dnbinom(X, size=par[1], mu=par[2], log=TRUE) * postPr[k,])}
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
List Mstep(NumericMatrix Fpr, NumericMatrix Gpr, NumericMatrix postPr,NumericMatrix trans,NumericVector X,NumericMatrix eParam){

  int nbT = Fpr.ncol(); 
  int nbS = Fpr.nrow();
  NumericMatrix Pi_tmp(nbS, nbS);
  // NumericMatrix eParam(nbS,2);

  for(int k=0;k<nbS;k++){
    for(int l=0;l<nbS;l++){
      NumericMatrix f=Fpr(Range(k,k),Range(0,nbT-2));
      NumericMatrix p=postPr(Range(l,l),Range(1,nbT-1));
      NumericMatrix g=Gpr(Range(l,l),Range(1,nbT-1));
      Pi_tmp(k,l) = trans(k,l)*sum( f*p/g );
    }
    // double Mean_tmp = sum(postPr(k,_) * X) / sum(postPr(k,_));
    // double Var_tmp  = sum(postPr(k,_) * pow((X - Mean_tmp),2)) / sum(postPr(k,_));
    // 
    // eParam(k,0)=Mean_tmp;
    // eParam(k,1)=Var_tmp;
  }
  
  List res;
  Rcpp::Function optim("optim");
  for(int k=0;k<nbS;k++){
    NumericVector par=eParam(k,_);
    NumericVector postPrv=postPr(k,_);
    // objFunc <- function(par){-sum(dnbinom(X, size=par[1], mu=par[2], log=TRUE) * postPr[k,])}
     res=optim(Rcpp::_["par"]    = par,
                   Rcpp::_["fn"]     = Rcpp::InternalFunction(&objFunc),
                   Rcpp::_["method"] = "BFGS",
                   Rcpp::_["X"]=X,
                   Rcpp::_["postPrv"]=postPrv);
    // eSize[k] = eParam[1]; esAvg[k] = eParam[2]
    NumericVector newpar=res["par"];
    eParam(k,_)=newpar;

  }

  // NumericMatrix trans_m = Pi_tmp/rowSums(Pi_tmp);
  NumericMatrix trans_m(nbS, nbS);
  for(int i=0;i<nbS;i++){
    NumericVector tmp=rowSums(Pi_tmp);
    trans_m(i,_)=Pi_tmp(i,_)/tmp[i];
  }
  // return List::create(_["hi"]=1);
// # ordMean = order(eParam[,1])
// # init = postPr[,1][ordMean]
// # eParam = eParam[ordMean,]
// # trans = trans[ordMean, ordMean]
  NumericVector init = postPr(_,1);
  return List::create(_["init"]=Fpr,_["init"]=Gpr,_["trans"]=trans_m,_["eParam"]=eParam);
}
List EM(NumericVector X, int nbS, int iterMax, double EPSILON){
  List resInitHMM = initHMM(X,nbS);
  // return resInitHMM;
  NumericVector init    = resInitHMM["init"];
  NumericMatrix trans   = resInitHMM["trans"];
  NumericMatrix eParam  = resInitHMM["eParam"];
  // NumericVector Param   = (init, as.vector(trans), as.vector(eParam));
  NumericVector Param;
  for(int i=0;i<init.size();i++)Param.push_back(init[i]);
  for(int i=0;i<trans.size();i++)Param.push_back(trans[i]);
  for(int i=0;i<eParam.size();i++)Param.push_back(eParam[i]);
  NumericMatrix Fpr;
  NumericMatrix Gpr;
  NumericMatrix postPr;
  double logLik;
  int iter;
  NumericMatrix emisPr;

  for(iter=0;iter<iterMax;iter++){
     emisPr = EmisPr(X, eParam);

    List resEstep = Estep(init, trans, emisPr,eParam);
    NumericMatrix f=resEstep["Fpr"];
    NumericMatrix g=resEstep["Gpr"];
    NumericMatrix p=resEstep["postPr"];
     Fpr  = f;
     Gpr  = g;
     postPr  = p;
     logLik  = resEstep["logLik"];

    List resMstep = Mstep(Fpr, Gpr, postPr, trans, X,eParam);

    NumericVector init_new   = resMstep["init"];
    NumericMatrix trans_new  = resMstep["trans"];
    NumericMatrix eParam_new = resMstep["eParam"];

    // NumericVector Param_new (init_new, as.vector(trans_new), as.vector(eParam_new));
    NumericVector Param_new;

    for(int i=0;i<init_new.size();i++)Param_new.push_back(init_new[i]);
    for(int i=0;i<trans_new.size();i++)Param_new.push_back(trans_new[i]);
    for(int i=0;i<eParam_new.size();i++)Param_new.push_back(eParam_new[i]);
    double D = sum(pow((Param_new - Param),2));
    Param  = Param_new;
    init   = init_new;
    trans  = trans_new;
    eParam = eParam_new;
    Rcpp::Function order("order");
    NumericVector ord=order(eParam(_,0));
    eParam(0,_) = reorder(eParam(0,_),ord);
    eParam(1,_) = reorder(eParam(1,_),ord);

    if(D<EPSILON){break;};
  }
 
  return List::create(_["initPr"]=init,_["transPr"]=trans,_["eParam"]=eParam,
                      _["postPr"]=postPr,_["logLik"]=logLik,_["emisPr"]=emisPr,_["iter"]=iter+1);
}

NumericVector Viterbi(NumericVector initPr,NumericMatrix transPr,NumericMatrix emisPr){
  int nbT = emisPr.ncol();
  int nbS = emisPr.nrow();
  NumericMatrix logF(nbS,nbT);
  logF(_,0) = log(initPr) + log(emisPr(_,0));
  
  for (int tt=1;tt<nbT;tt++) {
    for (int k=0;k<nbS;k++) {
      logF(k,tt) = log(emisPr(k,tt)) +  max(logF(_,tt-1) + log(transPr(_,k)));
    }
  }
  NumericVector path(nbT);
  path[nbT-1] = which_max(logF(_,nbT-1));
  
  for (int tt=nbT-2;tt>=0;tt--) {
    path[tt] = which_max(logF(_,tt) + log(transPr(_,path[tt+1])));
  }
  return path+1;
}

// [[Rcpp::export]]
List HMM_NB_NS(NumericVector X, int nbS, bool viterbi=true,int iterMax=5000, double EPSILON=1e-7){
  
  List resEM= EM(X, nbS, iterMax, EPSILON);
  NumericMatrix postPr = resEM["postPr"];
  NumericVector initPr = resEM["initPr"];
  NumericMatrix transPr = resEM["transPr"];
  NumericMatrix eParam  = resEM["eParam"];
  NumericMatrix emisPr  =resEM["emisPr"];
  int iterStop = resEM["iter"];
  double logLik  = resEM["logLik"];
  int nbT=X.size();
  NumericVector path(nbT);
  if(viterbi){
    path= Viterbi(initPr,transPr,emisPr);
  }else{
    // NumericVector path= apply(postPr,2,which.max);
    
    for(int i=0;i<nbT;i++){
      path[i]=which_max(postPr(_,i))+1;
    }
  }
  return List::create(_["initPr"]=initPr,_["transPr"]=transPr,_["eParam"]=eParam,_["postPr"]=postPr,_["logLik"]=logLik,_["emisPr"]=emisPr,_["iterStop"]=iterStop,_["path"]=path);
}
