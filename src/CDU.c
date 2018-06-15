//
//  CDU.c
//  
//
//  Created by Brad Price on 3/14/16.
//
//

#include "CDU.h"

void CDU(double *beta, double *xx, double *xy, double *y_clusters, double *delta, double *gamma_y, double *eps, int *miter, int *r, int *p, double *beta0, double *set,double *First,double *Second, double *mine){
    
    int nr=r[0];
    int np=p[0];
    int m,c,j,h,k,s,l;
    double diff=1;
    double old=0;
    int iter=0;
    /*double First;
    double Second;*/
    
    for(m=0; m<nr; m++){
        set[m]=0;
        for(c=0; c<nr; c++){
            if(y_clusters[c]==y_clusters[m]){
                set[m]=set[m]+1;
            }
        }
    }

    while(diff>eps[0] && iter<miter[0]){
        diff=0;
        old=0;
        iter=iter+1;
        for(j=0; j<nr; j++){
            for(h=0; h<np; h++){
                First[0]=0;

                for(k=0;k<np;k++){
                    if(k != h){
                        First[0]=First[0]+xx[h*np+k]*beta[j*np+k];
                    }
                }
                
                Second[0]=0;
                for(s=0; s<nr;s++){
                    if(y_clusters[s]==y_clusters[j] && s !=j){
                        for(l=0; l<np; l++){
                            Second[0]=Second[0]+xx[h*np+l]*beta[s*np+l];
                        }
                    }
                }
                mine[j*np+h]=beta[j*np+h];
                /* This is where the code breaks */
                beta0[j*np+h]=copysign(fmax(0,fabs(xy[j*np+h]-(1+(gamma_y[0]*(set[j]-1)/set[j]))*First[0]+(gamma_y[0]/set[j])*Second[0])-delta[0]/2),(xy[j*np+h]-(1+(gamma_y[0]*(set[j]-1)/set[j]))*First[0]+(gamma_y[0]/set[j])*Second[0]));
                beta0[j*np+h]=beta0[j*np+h]/((1+(gamma_y[0]*(set[j]-1)/set[j]))*xx[h*np+h]);
                /*up[j*np+h]=xy[j*np+h]-(1+(gamma_y[0]*(set[j]-1)/set[j]))*First[0]+(gamma_y[0]/set[j])*Second[0];
                beta[j*np+h]=copysign(fmax(0,(fabs(up)-delta[0]/2)),up)/(1+(gamma_y[0]*(set[j]-1)/set[j])*xx[h*np+h]);*/
                beta[j*np+h]=beta0[j*np+h];
                old=old+pow(mine[j*np+h],2);
                diff=diff+pow(mine[j*np+h]-beta[j*np+h],2);
            }
        }
        diff=diff/(old+.0000001);
    }
    /*up[0]=up2;*/
}
