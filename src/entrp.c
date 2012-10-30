#include "heap.h"
#include "entrp.h"
#include <stdlib.h>

/* #define VERBOSE 1 */

void entrp(double *A,
	   int *lda,
	   int *na,
	   int *nf,
	   int *ne,
	   int *ns, 
	   int *S,
	   double *opt,
	   int *S_Work,
	   int *F, 
	   int *E,
	   int *ind,
	   int *ind1,
	   double *As,
	   int *ldas,
	   double *Bs,
	   int *ldbs,
	   double *Cs,
	   int *ldcs,
	   double *Inv,
	   int *ldinv,
	   double *W,
	   double* WORK,
	   int *LWORK,
	   int *IWORK,
	   double *tol,
	   int *maxcount,
	   int *iter,
	   int *verbose
	   )

{ 
    double diff;
    double LB,UB,det,RCOND,UB_Work; 
    heap_type *subproblems;
    heap_element_type *act_subproblem, *act_subproblem_Work;
    int i,j,k,ni,cardf,carde,ierr,lerr;
    *iter=0;lerr=0;
	LB=*opt;

#ifdef VERBOSE
    if(*verbose!=0)
	{
	    printf("Init: \n");
	    printf("LB:%e\n",LB);
	}
#endif

    for(i=0;i<=(*nf)-1;i++)
	ind[i]=i+1;
    
    for(i=0;i<=(*ne)-1;i++)
	{
	    E[i]=1;
	    F[i]=0;
	}

    cardf=*nf;
    carde=*ne;
    if (cardf!=0)
	{
	    UB_Work=F77_CALL(upbnd)(A,
				      lda,
				      na,
				      F,
				      E,
				      ne,
				      ns,
				      As,
				      ldas,
				      Bs,
				      ldbs,
				      Cs,
				      ldcs,
				      Inv,
				      ldinv,
				      ind1,
				      &ierr,
				      W,
				      WORK,
				      LWORK,
				      IWORK); 
	    UB=UB_Work;
	}
    if (cardf==0)
	UB=2*LB;
    
#ifdef VERBOSE
    if(*verbose!=0)
	{
	    printf("\nbound\n"); 
	    printf("UB:%e\n",UB);
	}
#endif
    if(! (subproblems=heap_init())){
	/* error("memory exhausted"); exit(1); */
        lerr=1;
    }
    if(! (act_subproblem=heap_element_init(F,E,*ne,UB,0))){
	/* error("memory exhausted"); exit(1); */
        lerr=1;
    }


    /* heap_element_print(act_subproblem);  */
    heap_insert(subproblems, act_subproblem);    
    /* heap_traverse(subproblems); */
    
    diff=UB-LB;
#ifdef VERBOSE
    if(*verbose!=0)
	{
	    printf("\nbegin of the loop!! tolerance: %e\n",*tol);  
	}
#endif
    while(diff>*tol && subproblems->count>0)
	{
	    (*iter)++; 
	    
#ifdef VERBOSE
	    if(*verbose!=0)
	 	{
		    printf("iter: %i, heap: %i\n",*iter,subproblems->count); 
		}
#endif
	    /* heap_traverse(subproblems);  */
	    
	    act_subproblem=heap_drop(subproblems);
	    if(act_subproblem!=NULL)
		{
		    /* printf("dropped element:\n"); */
		    /* heap_element_print(act_subproblem); */

		    carde=0;
		    cardf=*nf;
		    for(j=0;j<=(*ne)-1;j++)
			{
			    carde=carde+act_subproblem->E[j]*act_subproblem->E[j];
			    cardf=cardf+act_subproblem->F[j]*act_subproblem->F[j];
			}
		    /* if(*verbose!=0)*/
		    /*   {*/
		    /*	    F77_CALL(ivecpr)(act_subproblem->E,ne); */
		    /*	    F77_CALL(ivecpr)(act_subproblem->F,ne); */
		    /*	}*/


		    for(j=0;j<=(*ne)-1;j++)
			S_Work[j]=0;

		    i=act_subproblem->ibranch;

		    /***** loop *****/
		    /*		    while(i<(*ne)-1 && diff>*tol)
				    { */
		    do
			i++;
		    while(act_subproblem->E[i]==0 && i<(*ne)-1);
		    if(act_subproblem->E[i]>0)
			{
			    
			    if(! (act_subproblem_Work=heap_element_init(NULL,NULL,*ne,0,0))){
				/* error("memory exhausted"); exit(1);*/  
                                lerr=1;
                            }
			    heap_element_copy(act_subproblem_Work,act_subproblem);
			    
			    act_subproblem_Work->E[i]=0;
			    act_subproblem_Work->ibranch=i;
			    
			    if(cardf+carde-1>(*ns)+(*nf))
				{
				    UB_Work=F77_CALL(upbnd)(A,
							      lda,
							      na,
							      act_subproblem_Work->F,
							      act_subproblem_Work->E,
							      ne,
							      ns,
							      As,
							      ldas,
							      Bs,
							      ldbs,
							      Cs,
							      ldcs,
							      Inv,
							      ldinv,
							      ind1, 
							      &ierr,
							      W,
							      WORK,
							      LWORK,
							      IWORK); 
				    act_subproblem_Work->b=UB_Work;
				    /* fathoming by lower bounds */
				    if(UB_Work>=LB) 
					heap_insert(subproblems, act_subproblem_Work);
				    /*	    printf("LB: %e UB: %e\n",LB, UB_Work); */
				}
			    if(cardf+carde-1==(*ns)+(*nf))
				{
				    for(j=0;j<=(*ne)-1;j++)
					S_Work[j]=(act_subproblem_Work->F[j]==1 ||
						   act_subproblem_Work->E[j]!=0 ?  j+(*nf)+1 : 0);
				    
				    ni=*nf;
				    for(j=0;j<=(*ne)-1;j++)
					{
					    if(S_Work[j]!=0)
						{
						    ind[ni]=S_Work[j];
						    ni=ni+1;
						}
					}
				    
				    /* calculation of the determinant */
				    F77_CALL(psubm)(A,lda,na,As,ldas,ind,&ni);
				    F77_CALL(chol1)(As,ldas,&ni,Bs,ldbs,&ierr);
				    det=F77_CALL(chdet)(Bs,ldbs,&ni); 
				    if (det>LB) 
					{
					    LB=det;
					    for(j=0;j<=(*ne)-1;j++)
						S[j]=S_Work[j];
					}
				}
							    
			    if(! (act_subproblem_Work=heap_element_init(NULL,NULL,*ne,0,0))){
				/* error("memory exhausted"); exit(1);*/  
                                lerr=1;
                            }
			    
			    heap_element_copy(act_subproblem_Work,act_subproblem);
			    act_subproblem_Work->E[i]=0;
			    act_subproblem_Work->F[i]=1;
			    act_subproblem_Work->ibranch=i;
			    
			    
			    if(cardf+1<(*ns)+(*nf))
				{
				    UB_Work=F77_CALL(upbnd)(A,
							      lda,
							      na,
							      act_subproblem_Work->F,
							      act_subproblem_Work->E,
							      ne,
							      ns,
							      As,
							      ldas,
							      Bs,
							      ldbs,
							      Cs,
							      ldcs,
							      Inv,
							      ldinv,
							      ind1,
							      &ierr,
							      W,
							      WORK,
							      LWORK,
							      IWORK); 
				    act_subproblem_Work->b=UB_Work;
				    /* fathoming by lower bounds */
				    /*  printf("LB: %e UB: %e\n",LB, UB_Work); */
				    if(UB_Work>=LB)
					heap_insert(subproblems, act_subproblem_Work);
				}
			    if(cardf+1==(*ns)+(*nf))
				{
				    for(j=0;j<=(*ne)-1;j++)
					S_Work[j]=(act_subproblem_Work->F[j]==1 ? j+(*nf)+1 : 0);
				    
				    ni=*nf;
				    for(j=0;j<=(*ne)-1;j++)
					{
					    if(S_Work[j]!=0)
						{
						    ind[ni]=S_Work[j];
						    ni=ni+1;
						}
					}
				    /* calculation of the determinant */
				    F77_CALL(psubm)(A,lda,na,As,ldas,ind,&ni);
				    F77_CALL(chol1)(As,ldas,&ni,Bs,ldbs,&ierr);
				    det=F77_CALL(chdet)(Bs,ldbs,&ni);
				    if (det>LB) 
					{
					    LB=det;
					    for(j=0;j<=(*ne)-1;j++)
						S[j]=S_Work[j];
					}
				}
			    UB=subproblems->top->down->b;
			    diff=UB-LB;
			    
#ifdef VERBOSE
			    if (subproblems->count!=0 && *verbose!=0)
				printf("diff:%e LB: %e UB: %e\n",diff, LB, UB);
#endif
			    
			    
			} 
		    
		    /*			}*/
		    /***** End While *****/
		}
	}
#ifdef VERBOSE
    if(*verbose!=0)
	{
	    printf("\nEnd\n");
	    /* printf("diff:%e LB: %e UB: %e\n",diff, LB, UB); */
	}
#endif
    
    /* LB is numerically better than UB !! */
    *opt=LB;
    *maxcount=subproblems->maxcount;
    
}
