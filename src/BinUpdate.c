//
//  BinUpdate.c
//  
//
//  Created by Brad Price on 3/28/16.
//
//

#include "BinUpdate.h"

void BinUp(double *TmV ,double *Tm, double *V, double *C, double *R, double *Beta, int* r, int* p, int* n, double *y_clusters, double *sl,
           double *gamma_y, double *delta, double *iter, double *eps,double *third, double *fourth, double *sixth, double *first, double* second, double* fifth){
    
    int nr=r[0];
    int np=p[0];
    //int nn=n[0];
    //double diff=1;
    //double old=0;
    //int iter=0;
    int j,l,k,s;
    double First;
    double Second;
    double Third;
    double Fourth;
    double Fifth;
    double Sixth;
    //double Top;
    //double Bottom;
    double New=10;
    double it=0;
    double Est;
    
  while(New>eps[0] && it<iter[0]){
        New=0;
        for(l=0; l<nr; l++){
            for(k=0; k<np; k++){
                First=0;
                Second=0;
                Third=0;
                Fourth=0;
                Fifth=0;
                Sixth=0;
                Est=0;
               // if(k==0){
                 
               //  First=TmV[(l*np)+k];
                 
                 //   for(i=0; i<nn; i++){
                //      First=First+Tm[(l*nn*np)+(k*nn)+i]*V[(l*nn)+i];
                  //  }
                //    for(j=0; j<np; j++){
                //        if(j!=k){
                //            Second=Second+C[(l*np*np)+(k*np)+j]*Beta[(l*np)+j];
                //        }
                //    }
                    
                //    Fifth=C[(l*np*np)+(k*np)+k];
                //    Est=(First-Second)/Fifth;
                //    New=New+(pow(Beta[(l*np)+k]-Est,2)/(np*nr));
                //    Beta[l*np+k]=Est;
                //}
                //if(k!=0){
                    
                    First=TmV[(l*np)+k];
                    //for(i=0; i<nn; i++){
                    //    First=First+Tm[(l*nn*np)+(k*nn)+i]*V[(l*nn)+i];
                    //}
                    for(j=0; j<np; j++){
                        if(j!=k){
                            Second=Second+C[(l*np*np)+(k*np)+j]*Beta[(l*np)+j];
                            Third=Third+R[(k*np)+j]*Beta[(l*np)+j];
                        }
                        for(s=0; s<nr;s++){
                            if(y_clusters[s]==y_clusters[l] && s!=l){
                                Fourth=Fourth+R[(k*np)+j]*Beta[(s*np)+j];
                            }
                        }
                        
                    }
                    
                    
                    
                    Fifth=C[(l*np*np)+(k*np)+k];
                    Sixth=R[k*np+k];
                if(k==0){
                    Est=(First-Second-(gamma_y[0]*(sl[l]-1)/sl[l])*Third+(gamma_y[0]/sl[l])*Fourth);
                }
                if(k!=0){
                    Est=copysign(fmax(0,fabs(First-Second-(gamma_y[0]*(sl[l]-1)/sl[l])*Third+(gamma_y[0]/sl[l])*Fourth)-delta[0]/2),
                                 (First-Second-(gamma_y[0]*(sl[l]-1)/sl[l])*Third+(gamma_y[0]/sl[l])*Fourth));
                }
                    Est=Est/(Fifth+(gamma_y[0]*(sl[l]-1)/sl[l])*Sixth);
                    New=New+pow((Beta[l*np+k]-Est)/pow(nr*np,.5),2);
                    Beta[l*np+k]=Est;
                // }
                first[l*np+k]=First;
                second[l*np+k]=Second;
                third[l*np+k]=Third;
                fourth[l*np+k]=Fourth;
                fifth[l*np+k]=Fifth;
                sixth[l*np+k]=Sixth;
            }
            
        }
        it=it+1;
        }
    
}
