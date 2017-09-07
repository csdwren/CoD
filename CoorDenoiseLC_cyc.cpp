#include "mex.h"
#include "ingredients.cpp"
#include <time.h>

/*
 *using coordinate descent method to directly solve both anisotropic and 
 * isotropic TV image denoising with the formulation,
 *       0.5||x-y||^2+lambda||Dx||
 *e.g. x=CoorDenoiseC_cyc(y,lambda,maxIter,L,perturbation)
 *
 */

void mexFunction(int nlhs,mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int i,j,k,l,iter,maxIter,ii,   count=0;
    int i0,j0,i1,j1,i2,j2,i3,j3,i4,j4;
    int m,n;
    double *y,*x,lambda,dchange,tmp[5];
    bool *tilde_x;
    double *theta,*sintheta,*costheta,*a,*e,*sa,dL=0,Axsp;
    int L,kk;
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
    L=(int)mxGetScalar(prhs[3]);
    permu=(double)mxGetScalar(prhs[4]);
    //double *ini=(double *)mxGetData(prhs[5]);
//     int fla=(int)mxGetScalar(prhs[6]);
    
    //set outputs auguments
    plhs[0] = mxCreateNumericArray(2, outsize, mxDOUBLE_CLASS, mxREAL);//output matrix
    x = (double *)mxGetPr(plhs[0]);
    memcpy(x,y,(m*n*sizeof(double)));//initilize x=y
    
    //theta
    theta=(double*)malloc(L*sizeof(double));
    sintheta=(double*)malloc(L*sizeof(double));
    costheta=(double*)malloc(L*sizeof(double));
    for(i=0;i<L;++i)
    {
        theta[i]=PI*i/2/L;
        sintheta[i]=sin(theta[i]);
        costheta[i]=cos(theta[i]);
        dL+=costheta[i]+sintheta[i];
    }
    
    // coeficients vector for subproblem
    a=(double*)malloc((4*L+1)*sizeof(double));
    e=(double*)malloc((4*L+1)*sizeof(double));
    sa=(double*)malloc((4*L+1)*sizeof(double));
    
    //strore solutions of each subproblem
    tilde_x=(bool *)malloc((m*n)*sizeof(bool));
    memset(tilde_x,0,m*n*sizeof(bool));
    
    double minper=permu/100000000;
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
        permu/=1.18;
     //   if(permu<minper) permu=minper;
        
        
        for(ii=0;ii<m*n;++ii)
        {
            j=(int)(ii/m);
            i=ii-j*m;
            
            
            // if(iter==0)// calculate the entire solution to tilede_x.
            if(tilde_x[j*m+i]==0)
            {count++;
//                 if(i<1 || i>m-2 || j<1 || j>n-2)
                {
                SwSubscript(i-1,j,m,n,&i1,&j1);
                SwSubscript(i,j-1,m,n,&i2,&j2);
                SwSubscript(i+1,j,m,n,&i3,&j3);
                SwSubscript(i,j+1,m,n,&i4,&j4);
                }
//                 else{
//                     i1=i-1;j1=j;
//                     i2=i;j2=j-1;
//                     i3=i+1;j3=j;
//                     i4=i;j4=j+1;
//                 }
//                 tmp[1]=x[j1*m+i1];tmp[2]=x[j2*m+i2];tmp[3]=x[j3*m+i3];tmp[4]=x[j4*m+i4];
//                 tilde_x[j*m+i]=SolveEleSub(y[j*m+i],lambda,tmp,4);
                
                // solve subproblem. vectors a and e are first constructed
                for(k=0;k<L;++k)
                {
                    a[k*4+1]=costheta[k]+sintheta[k];
                    a[k*4+2]=a[k*4+1];
                    e[k*4+1]=x[j2*m+i2]*costheta[k]+x[j1*m+i1]*sintheta[k];
                    e[k*4+2]=x[j4*m+i4]*costheta[k]+x[j3*m+i3]*sintheta[k];
                    
                    Axsp=costheta[k]-sintheta[k];
                    if(Axsp>0)
                    {
                        a[k*4+3]=costheta[k]-sintheta[k];
                        a[k*4+4]=a[k*4+3];
                        e[k*4+3]=x[j1*m+i1]*costheta[k]-x[j2*m+i2]*sintheta[k];
                        e[k*4+4]=x[j3*m+i3]*costheta[k]-x[j4*m+i4]*sintheta[k];
                    }
                    else
                    {
                        a[k*4+3]=-costheta[k]+sintheta[k];
                        a[k*4+4]=a[k*4+3];
                        e[k*4+3]=-x[j1*m+i1]*costheta[k]+x[j2*m+i2]*sintheta[k];
                        e[k*4+4]=-x[j3*m+i3]*costheta[k]+x[j4*m+i4]*sintheta[k];
                    }
                }
//             tmp[1]=x[j1*m+i1];tmp[2]=x[j2*m+i2];tmp[3]=x[j3*m+i3];tmp[4]=x[j4*m+i4];
                x_pre=x[j*m+i];
                
                x[j*m+i]=SolveEleSub(y[j*m+i],lambda/dL,e,4*L,a,sa);
                
                tilde_x[j*m+i]=1;
                if(fabs(x_pre-x[j*m+i])>1e-4 && iter<maxIter-1)
                {
                    double sr=(double)(1*rand()/(RAND_MAX+1.0))-0.5;
                  //  if(iter%2==1) sr=1e-3;
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


