#include "mex.h"
#include "ingredients.cpp"
#include "twomedian.cpp"
#include "qmedian.cpp"
#include <time.h>

/*
 *using coordinate descent method to directly solve anisotropic TV image denoising
 *with the formulation,
 *       0.5||x-y||^2+lambda||Dx||_1.
 *e.g. x=CoorDenoiseC(y,lambda,Tol)
 *y is the degraded image
 */

void mexFunction(int nlhs,mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int i,j,k,l,iter,maxIter,ii,   count=0;
    int i0,j0,i1,j1,i2,j2,i3,j3,i4,j4;
    int m,n;
    double *y,*x,lambda,dchange,tmp[9];
    bool *tilde_x;
//     double *theta,*sintheta,*costheta,*a,*e,*sa,dL=0,Axsp;
//     int L,kk;
    double permu,x_pre;
    
    int outsize[2];
    
    m=(int)mxGetM(prhs[0]);
    n=(int)mxGetN(prhs[0]);
    outsize[0]=m;
    outsize[1]=n;
    
    //get inputs arguments
    y=(double *)mxGetData(prhs[0]);
    lambda=(double)mxGetScalar(prhs[1]);
    maxIter=(int)mxGetScalar(prhs[2]);
    permu=(double)mxGetScalar(prhs[3]);
//     double *ini=(double *)mxGetData(prhs[5]);
//     int fla=(int)mxGetScalar(prhs[6]);
    
    //set outputs auguments
    plhs[0] = mxCreateNumericArray(2, outsize, mxDOUBLE_CLASS, mxREAL);//output matrix
    x = (double *)mxGetPr(plhs[0]);
    memcpy(x,y,(m*n*sizeof(double)));//initilize x=y
    
//     consts=(double **)malloc(m*n*sizeof(double*));
//     for(i=0;i<m*n;++i)
//     {
//         consts[i]=(double *)malloc(9*sizeof(double));
//
//         consts[i][4]=y[i]-4*lambda;
//         consts[i][5]=y[i]-2*lambda;
//         consts[i][6]=y[i];
//         consts[i][7]=y[i]+2*lambda;
//         consts[i][8]=y[i]+4*lambda;
//
//     }
    
    
    //strore solutions of each subproblem
    tilde_x=(bool *)malloc((m*n)*sizeof(bool));
    memset(tilde_x,0,m*n*sizeof(bool));
    
//     for(i=0;i<m;++i)
//     {
//         for(j=0;j<n;++j)
//         {
//             SwSubscript(i-1,j,m,n,&i1,&j1);
//             SwSubscript(i,j-1,m,n,&i2,&j2);
//             SwSubscript(i+1,j,m,n,&i3,&j3);
//             SwSubscript(i,j+1,m,n,&i4,&j4);
//             tmp[1]=x[i1*n+j1];tmp[2]=x[i2*n+j2];tmp[3]=x[i3*n+j3];tmp[4]=x[i4*n+j4];
//             tilde_x[i*n+j]=SolveEleSub(y[i*n+j],lambda,tmp,4);
//         }
//     }
    
    for(iter=0;iter<maxIter;++iter)
    {
//         if (fla==1)
        permu/=1.15;
        
        
        
        for(ii=0;ii<m*n;++ii)
        {
            j=(int)(ii/m);
            i=ii-j*m;
            
//             if(iter==0)
//             {
//                 consts[j*m+i]=(double *)malloc(9*sizeof(double));
//
//         consts[j*m+i][4]=y[j*m+i]-4*lambda;
//         consts[j*m+i][5]=y[j*m+i]-2*lambda;
//         consts[j*m+i][6]=y[j*m+i];
//         consts[j*m+i][7]=y[j*m+i]+2*lambda;
//         consts[j*m+i][8]=y[j*m+i]+4*lambda;
//             }
            
            
            // if(iter==0)// calculate the entire solution to tilede_x.
            if(tilde_x[j*m+i]==0)
            {count++;
             if(i<1 || i>m-2 || j<1 || j>n-2)
             {
                 SwSubscript(i-1,j,m,n,&i1,&j1);
                 SwSubscript(i,j-1,m,n,&i2,&j2);
                 SwSubscript(i+1,j,m,n,&i3,&j3);
                 SwSubscript(i,j+1,m,n,&i4,&j4);
             }
             else{
                 i1=i-1;j1=j;
                 i2=i;j2=j-1;
                 i3=i+1;j3=j;
                 i4=i;j4=j+1;
             }
//                 tmp[1]=x[j1*m+i1];tmp[2]=x[j2*m+i2];tmp[3]=x[j3*m+i3];tmp[4]=x[j4*m+i4];
//                 tilde_x[j*m+i]=SolveEleSub(y[j*m+i],lambda,tmp,4);
             
             // solve subproblem. vectors a and e are first constructed
             
             
             tmp[0]=x[j1*m+i1];tmp[1]=x[j2*m+i2];tmp[2]=x[j3*m+i3];tmp[3]=x[j4*m+i4];
             tmp[4]=y[j*m+i]-4*lambda;
             tmp[5]=y[j*m+i]-2*lambda;
             tmp[6]=y[j*m+i];
             tmp[7]=y[j*m+i]+2*lambda;
             tmp[8]=y[j*m+i]+4*lambda;
//              consts[j*m+i][0]=x[j1*m+i1];
//              consts[j*m+i][1]=x[j2*m+i2];
//              consts[j*m+i][2]=x[j3*m+i3];
//              consts[j*m+i][3]=x[j4*m+i4];
//
//              if(iter==1)
//              {consts[j*m+i][4]=y[j*m+i]-4*lambda;consts[j*m+i][5]=y[j*m+i]-2*lambda;consts[j*m+i][6]=y[j*m+i];consts[j*m+i][7]=y[j*m+i]+2*lambda;consts[j*m+i][8]=y[j*m+i]+4*lambda;}
             
             
             
             //mexPrintf("%.3lf %.3lf %.3lf %.3lf\n",tmp[0],tmp[1],tmp[2],tmp[3]);
             
             
//              system("pause");
             
             x_pre=x[j*m+i];
             
             x[j*m+i]=opt_med9(tmp);
             
             tilde_x[j*m+i]=1;
             if(fabs(x_pre-x[j*m+i])>1e-5)
             {
                 double sr=(double)(1*rand()/(RAND_MAX+1.0))-0.5;
//                      cout<<sr<<endl;
                 x[j*m+i]+=sr*permu;
                 
                 tilde_x[j1*m+i1]=0;tilde_x[j2*m+i2]=0;tilde_x[j3*m+i3]=0;tilde_x[j4*m+i4]=0;
             }
             
             
            }
        }
    }
    free(tilde_x);
//     mexPrintf("%d\n",count);
}


