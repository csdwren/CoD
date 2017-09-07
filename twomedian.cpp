#include <iostream>
#include <time.h>
#include "mex.h"

using namespace std;

int n;

void INSERTION_SORT(double A[],int r)//
{
    int key;
    for (int j=1;j<=r;j++)
    {
        key=A[j-1];
        int i=j-1;
        while (i>0&&A[i-1]>key)//
        {
            A[i]=A[i-1];
            i=i-1;
        }
        A[i]=key;
    }
}

double Two_groups_array_Median(double A[],double B[],int beginA,int endA,int beginB,int endB)
{
//     int n=endA+1;
   int p=(beginA+endA)/2;int r=(beginB+endB)/2;
   int t=p+r+2,flag=0;//
   if (A[n-1]<B[0])//
   {
        return A[n-1];
   }
   else if (A[0]>B[n-1])
   {
        return B[n-1];
   }
   if (beginA>=endA||beginB>=endB)//
   {
       if (n-t>0)
       {
           for (int i=1;i<=n-t;i++)//
           {
               ++p;++r;
               if (A[p]<B[r])
               {
                   --r;
                   flag=1;
               }
               else
               {
                   --p;
                   flag=0;
               }
           }
           if (flag)//
           {
               return A[p];
           }
           else
           {
               return B[r];
           }
       }
       else//
       {
           if (A[p]<B[r])
           {
               return B[r];
           }
           else
           {
               return A[p];
           }
       }
   }
   if (A[p]>B[r])//
   {
      return Two_groups_array_Median(A,B,beginA,p,r,endB);//
   }
   else //
   {
      return Two_groups_array_Median(A,B,p,endA,beginB,r);
   }
}

// void mexFunction(int nlhs,mxArray *plhs[], int nrhs, const mxArray *prhs[])
// {
//     
//     int outsize[2],m;
//     outsize[0]=1;outsize[1]=1;
//     
//     m=(int)mxGetM(prhs[0]);
//     n=(int)mxGetN(prhs[0]);
// 
//     if(m>n) n=m;
//     
//     //get inputs arguments
//     double *a=(double *)mxGetData(prhs[0]);
//     double *b=(double *)mxGetData(prhs[1]);
//     
//     
// //     double *tttt=(double*)mxGetData(prhs[5]);
//     
//     //set outputs auguments
//     plhs[0] = mxCreateNumericArray(2, outsize, mxDOUBLE_CLASS, mxREAL);//output matrix
//     double *x = (double *)mxGetPr(plhs[0]);
//     
// //     memset(x,0,m*n*sizeof(double));
//     x[0]=Two_groups_array_Median(a,b,0,n-1,0,n-1);
// 
// }