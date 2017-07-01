#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <time.h>

typedef struct{
  size_t l;
  double *F_em1;
  double *F_em2;
  double *F_diff1;
  double *F_diff2;
  double *F_diff3;
  double *F_tmp;
  double stepMax;
}emaccl_internal;


typedef struct{
  emaccl_internal *internal;//keeps track of internal variable and allocations
  double (*llhFP)(double *,void *);//llh function. Log likelihood. This should be >0. and we seek the parameters that minimizes this function.
  int (*emstepFP)(double *,void*,double*);//function that returns the next guess.
  int verbose;//should we printout stuff
  double ttol;//tolerance used in accelerated em
  double tole;//tolerance used for difference in llh values
  double l;//dimension of parameterspase
  double type;//type=0 normal em, type=1 accelerated em.
  double maxIter;
}emaccl_pars;


emaccl_internal *emaccl_internal_alloc(int ndim){
  emaccl_internal *ret = malloc(sizeof(emaccl_internal)); 
  ret->l=ndim;
  ret->F_em1 = malloc(ret->l*sizeof(double));
  ret->F_em2 = malloc(ret->l*sizeof(double));
  ret->F_diff1 = malloc(ret->l*sizeof(double));
  ret->F_diff2 = malloc(ret->l*sizeof(double));
  ret->F_diff3 = malloc(ret->l*sizeof(double));
  ret->F_tmp = malloc(ret->l*sizeof(double));
  ret->stepMax=1;
  return ret;
}

emaccl_pars *emaccl_pars_alloc(int ndim){
  emaccl_pars *ep = malloc(sizeof(emaccl_pars));
  ep->internal = emaccl_internal_alloc(ndim);
  ep->ttol=1e-9;
  ep->tole=1e-6;
  ep->l=ndim;
  ep->type=1;
  ep->maxIter = 100;
  ep->llhFP=NULL;
  ep->emstepFP=NULL;
  ep->verbose=100;
  return ep;
}

void minus(double *fst,double *sec,double *res,int l){
  for(int i=0;i<l;i++)
    res[i] = fst[i]-sec[i];
}

double sumSquare(double *mat, int l){
  double tmp=0;
  for(size_t i=0;i<l;i++){
    //    fprintf(stderr,"%f \n",mat[i]);
    tmp += mat[i]*mat[i];
  }

  return tmp;
}


int emAccel(double *F,void *vp,double *F_new,emaccl_pars *ep){
  //  fprintf(stderr,"tol:%f\n",tol);

  double stepMin =1;
  double *stepMax=&ep->internal->stepMax;
  double mstep=4;
  //  double objfnInc=1;


  double *F_em1=ep->internal->F_em1;
  double *F_em2=ep->internal->F_em2;
  double *F_diff1=ep->internal->F_diff1;
  double *F_diff2=ep->internal->F_diff2;
  double *F_diff3=ep->internal->F_diff3;
  double *F_tmp=ep->internal->F_tmp;

  ep->emstepFP(F,vp,F_em1);
  minus(F_em1,F,F_diff1,ep->l);
  double sr2 = sumSquare(F_diff1,ep->l);
  
  if(sqrt(sr2)<ep->ttol){
    fprintf(stderr,"sr2 break:%f\n",sr2);
    return 0;
  }

  ep->emstepFP(F_em1,vp,F_em2);
  minus(F_em2,F_em1, F_diff2,ep->l);

  double sq2 = sumSquare(F_diff2,ep->l);
  if(sqrt(sq2)<ep->ttol){
    fprintf(stderr,"sq2\n");
    return 0;
  }

  minus(F_diff2,F_diff1, F_diff3,ep->l);
  double sv2 = sumSquare(F_diff3,ep->l);
  
  double alpha = sqrt(sr2/sv2);
  double mmin=*stepMax<alpha?*stepMax:alpha;
  alpha =stepMin>mmin?stepMin:mmin;// ;std::max(stepMin,std::min(*stepMax,alpha));
  for(size_t i=0;i<ep->l;i++)
    F_new[i] = F[i]+2*alpha*F_diff1[i]+alpha*alpha*F_diff3[i];
  //  fprintf(stderr,"F_new (this the linear jump: (%f,%f)\n",F_new[0],F_new[1]);
  int outofparspace =0;
  for(int i=0;i<ep->l;i++){
    if(F_new[i]<0||F_new[i]>1)
      outofparspace++;
  }
  if(outofparspace){
    fprintf(stderr,"outofparspace will use second emstep as jump instead\n");
    for(int i=0;i<ep->l;i++)
      F_new[i] = F_em2[i];
    return 1;
  }
    
  
  if (fabs(alpha - 1) > 0.01){
    ep->emstepFP(F_new,vp,F_tmp);
    for(int i=0;i<ep->l;i++){
      double tmp = F_new[i];
      F_new[i]=F_tmp[i];
      F_tmp[i]=tmp;
    }
  }

  //  double lnew = 1;
  if ((alpha - *stepMax) > -0.001) 
    *stepMax = mstep**stepMax;

  
  //  print(stderr,3,F_new);
  if(ep->verbose>1)
    fprintf(stderr,"\t-> alpha %f stepMax %f\n",alpha,*stepMax);

  return 1;
}



int em1(double *sfs,void *vpp,void *vp){
  emaccl_pars *ep = vpp;
  assert(ep->llhFP);assert(ep->emstepFP);
  double oldLik,lik;
  oldLik = ep->llhFP(sfs,vp);
  if(ep->verbose){
    fprintf(stderr,"startlik=%f",oldLik);
    for(int i=0;i<ep->l;i++)
      fprintf(stderr,"\t%f",sfs[i]);
    fprintf(stderr,"\n");
    fflush(stderr);
  }


  double *tmp =malloc(ep->l*sizeof(double));
  int it;
  for(it=0;it<ep->maxIter;it++) {
    if(ep->type==0)
      ep->emstepFP(sfs,vp,tmp);
    else if(ep->type==1)
      emAccel(sfs,vp,tmp,ep);
    lik = ep->llhFP(tmp,vp);
    
    double pdif=0;
    for(int i=0;i<ep->l;i++)
      pdif +=(tmp[i]-sfs[i])*(tmp[i]-sfs[i]);
    if(ep->verbose){
      fprintf(stderr,"[%d] lik=%f diff_llh=%g sqrt(diff_pars):%e ",it,lik,lik-oldLik,sqrt(pdif));
      for(int i=0;i<ep->l;i++)
	fprintf(stderr,"\t %e",tmp[i]-sfs[i]);
      fprintf(stderr,"\n");
    }
    if(lik>oldLik){
      fprintf(stderr,"\t-> New likelihood from EM-step is worse, assuming optima has been achieved and we are observing an issue of numerical stability\n");
      break;
    }else{
      for(int i=0;i<ep->l;i++)
	sfs[i]= tmp[i];
      assert(lik<=oldLik);
      if(lik-oldLik==0||fabs(lik-oldLik)<ep->tole){
	if(ep->verbose)
	  fprintf(stderr,"Breaking in emalgo, since difference in llh are smaller than tolerance\n");
	oldLik=lik;
	break;
      }
      oldLik=lik;
    }
  }
  return it;
}


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

typedef struct{
  int n;
  double *d;
}testp;


double llh(double *f,void *vp){
  testp *p =(testp*) vp;
  double llh=0;
  for(int i=0;i<p->n;i++)
    llh += log( f[0]*exp(p->d[i*2])+(f[1])*exp(p->d[i*2+1]));
  return -llh;
}

int emstep(double *pre,void *vp,double *post){
  testp *p =(testp*) vp;
  double f2=0;
   for(int i=0;i<p->n;i++)
     f2 += pre[0]*exp(p->d[i*2])/(pre[0]*exp(p->d[i*2])+pre[1]*exp(p->d[2*i+1]));
   post[0] = f2/(1.0*p->n);
   post[1] = 1-post[0];
   //   fprintf(stderr,"diff1:%e\n",1-post[0]-post[1]);
   return 0;
}




int main(){
  clock_t t=clock();
  time_t t2=time(NULL);
  double f=0.000018;
  double err = 0.025;
  int dep =3;
  int n=10000000;
  double *d=simdata(f,n,err,dep);

  double start[2] = {0.5,0.5};
  testp tp;
  tp.d=d;tp.n=n;
  double oldlik  =llh(start,&tp);
  fprintf(stderr,"(%f,%f):llh1:%f\n",start[0],start[1],oldlik);

  emaccl_pars *ep=emaccl_pars_alloc(2);
  ep->verbose=1;ep->maxIter=2000;
  ep->llhFP=llh;
  ep->emstepFP=emstep;
  ep->type =1;
  em1(start,ep,&tp);

  fprintf(stderr, "\t[ALL done] cpu-time used =  %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
  fprintf(stderr, "\t[ALL done] walltime used =  %.2f sec\n", (float)(time(NULL) - t2));  

  t=clock();
  t2=time(NULL);

  start[0]=start[1]=0.5;
  ep->type =0;
  em1(start,ep,&tp);
  fprintf(stderr, "\t[ALL done] cpu-time used =  %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
  fprintf(stderr, "\t[ALL done] walltime used =  %.2f sec\n", (float)(time(NULL) - t2));  
      
  return 0;
}
