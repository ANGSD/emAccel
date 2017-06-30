#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>


typedef struct{
  size_t l;
  double *F_em1;
  double *F_em2;
  double *F_diff1;
  double *F_diff2;
  double *F_diff3;
  double *F_tmp;
}emaccl_internal;


typedef struct{
  emaccl_internal *internal;//keeps track of internal variable and allocations
  double (*llh)(double *,void *);//llh function. Log likelihood. This should be >0. and we seek the parameters that minimizes this function.
  double (*emstep)(double *,void*,double*);//function that returns the next guess.
  int verbose;//should we printout stuff
  double ttol;//tolerance
  double l;//dimension of parameterspase
  double type;//type=0 normal em, type=1 accelerated em.
}empars;


typedef struct{
  int n;
  double *d;
}testp;


double llh(double f,double *d,int n){
  double llh=0;
  for(int i=0;i<n;i++)
    llh += log( f*exp(d[i*2])+(1-f)*exp(d[i*2+1]));
  return -llh;
}

int emstep(double *pre,double *post,double *d,int n){
  double f2=0;
   for(int i=0;i<n;i++)
     f2 += pre[0]*exp(d[i*2])/(pre[0]*exp(d[i*2])+pre[1]*exp(d[2*i+1]));
   post[0] = f2/(1.0*n);
   post[1] = 1-post[0];
   //   fprintf(stderr,"diff1:%e\n",1-post[0]-post[1]);
   return 0;
}


emaccl_internal *emaccl_internal_alloc(int ndim){
  emaccl_internal *ret = malloc(sizeof(emaccl_internal)); 
  ret->l=ndim;
  ret->F_em1 = malloc(ret->l*sizeof(double));
  ret->F_em2 = malloc(ret->l*sizeof(double));
  ret->F_diff1 = malloc(ret->l*sizeof(double));
  ret->F_diff2 = malloc(ret->l*sizeof(double));
  ret->F_diff3 = malloc(ret->l*sizeof(double));
  ret->F_tmp = malloc(ret->l*sizeof(double));
  return ret;
}

#if 0

int em1(double *sfs,double  **emis,double tole,int maxIter,int len,int verbose){
  double oldLik,lik;
  oldLik = loglike(sfs,emis,len);
  if(verbose)
    fprintf(stderr,"startlik=%f %f %f %f\n",oldLik,sfs[0],sfs[1],sfs[2]);
  fflush(stderr);

  double tmp[3];
  int it;
  for(it=0;it<maxIter;it++) {
    emStep1(sfs,emis,tmp,len);
    for(int i=0;i<3;i++)
      sfs[i]= tmp[i];
    lik = loglike(sfs,emis,len);

    if(verbose)
      fprintf(stderr,"[%d] lik=%f diff=%g\n",it,lik,fabs(lik-oldLik));

    if(fabs(lik-oldLik)<tole){
     
      oldLik=lik;
      break;
    }
    oldLik=lik;
  }
  return it;
}
#endif

#if 0
int emAccel(double *F,double **emis,double *F_new,int len){
  //  fprintf(stderr,"calling emaccel \n");
   double ttol=0.0000001;

  //  fprintf(stderr,"tol:%f\n",tol);
  //maybe these should be usersettable?
  double stepMin =1;
  const   double stepMax0 = 1;
  static double stepMax=stepMax0;
  double mstep=4;
  //  double objfnInc=1;


  double F_em1[3];
  double F_diff1[3];
  double F_em2[3];
  double F_diff2[3];
  double F_diff3[3];
  double F_tmp[3];

  emStep1(F,emis,F_em1,len);
  // stayin(F_em1);
  minus(F_em1,F,F_diff1);
  double sr2 = sumSquare(F_diff1);
  
  if(sqrt(sr2)<ttol){
    //    fprintf(stderr,"sr2 break:%f\n",sr2);
    return 0;
    //break;
  }
  emStep1(F_em1,emis,F_em2,len);
  //  stayin(F_em2);
  minus(F_em2,F_em1, F_diff2);

  double sq2 = sumSquare(F_diff2);
  if(sqrt(sq2)<ttol){
    //fprintf(stderr,"sq2\n");
    return 0;
    //    break;
  }


  minus(F_diff2,F_diff1, F_diff3);
  
  double sv2 = sumSquare(F_diff3);
  
  double alpha = sqrt(sr2/sv2);
  alpha = std::max(stepMin,std::min(stepMax,alpha));
  for(size_t i=0;i<3;i++)
    F_new[i] = F[i]+2*alpha*F_diff1[i]+alpha*alpha*F_diff3[i];
  //  fprintf(stderr,"F_new (this the linear jump: (%f,%f,%f)\n",F_new[0],F_new[1],F_new[2]);
  int outofparspace =0;
  for(int i=0;i<3;i++){
    if(F_new[i]<0||F_new[i]>1)
      outofparspace++;
  }
  if(outofparspace){
    //    fprintf(stderr,"outofparspace will use second emstep as jump\n");
    for(int i=0;i<3;i++)
      F_new[i] = F_em2[i];
    return 1;
  }
    
  
  if (fabs(alpha - 1) > 0.01){
    emStep1(F_new,emis,F_tmp,len);
    //    stayin(F_tmp);
    for(int i=0;i<3;i++)
      std::swap(F_new[i],F_tmp[i]);
  }

  //  double lnew = 1;
  if ((alpha - stepMax) > -0.001) {
    stepMax = mstep*stepMax;
  }
  //  print(stderr,3,F_new);
  //  fprintf(stderr,"alpha %f stepMax %f\n",alpha,stepMax);

  // fprintf(stderr,"calling emaccel \n");
  // fprintf(stderr,"F_new:(%f,%f,%f)\n",F_new[0],F_new[1],F_new[2]);
  return 1;

  
}

#endif

void mean(int *obs,int n){
 double s=0;
  for(int i=0;i<n;i++)
    s+=obs[i];
  fprintf(stderr,"est:%f\n",s/(1.0*n));
}


double *simdata(double f,int n,double errate,int dep){
  int *obs=malloc(n*sizeof(int));
  for(int i=0;i<n;i++)
    obs[i]=drand48()<=f?1:0;
  mean(obs,n);

  double *ret =(double *)malloc(2*n*sizeof(double));
   for(int i=0;i<n;i++){
    double *persite=ret+2*i;
    persite[0]=persite[1]=0;

    int md = lrand48() % dep;
    //    fprintf(stderr,"%d/%d\n",md,dep );
    for(int d=0;d<md;d++){
      int val = obs[i];
      if(drand48()<=errate)
	val=val?0:1;
      
      if(val){
	persite[0]+=log(errate);
	persite[1]+=log(1-errate);
      }else{
	persite[1]+=log(errate);
	persite[0]+=log(1-errate);
      }
      double mmax=persite[0]>persite[1]?persite[0]:persite[1];
      persite[0]-=mmax;persite[1]-=mmax;
    }
  }
 
  return ret;;
}

int main(){
  double f=0.000018;
  double err = 0.025;
  int dep =3;
  int n=10000000;
  double *d=simdata(f,n,err,dep);
  for(int i=0;0&&i<n;i++)
    fprintf(stdout,"%f\t%f\n",d[i*2],d[i*2+1]);
  fflush(stdout);
  double start[2] = {0.5,0.5};
  double oldlik  =llh(start[0],d,n);
  fprintf(stderr,"(%f,%f):llh1:%f\n",start[0],start[1],oldlik);
  double tmp[2];
  for(int i=0;i<10000;i++){
    emstep(start,tmp,d,n);
    double newllh = llh(tmp[0],d,n);
    fprintf(stderr,"[%d]:(%f,%f):llh2:%f diff:%e\n",i,tmp[0],tmp[1],newllh,newllh-oldlik);
    assert(newllh<=oldlik);
    start[0]=tmp[0];start[1]=tmp[1];
    if(newllh-oldlik==0)
      break;

    oldlik=newllh;
  }
      
  return 0;
}
