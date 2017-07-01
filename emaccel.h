
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

#ifdef __cplusplus
extern "C" {
#endif

emaccl_pars *emaccl_pars_alloc(int ndim);
void emaccl_pars_free(emaccl_pars *ep);
int em1(double *sfs,void *vpp,void *vp);

#ifdef __cplusplus
}
#endif

