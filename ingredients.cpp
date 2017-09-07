
/*
 * Contain all the ingredients used in coordinate descent method for TV image restoration.
 */
 
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <algorithm> 
 
using namespace std;
 
#define Left(i)        ((i) << 1)
#define Right(i)    (((i) << 1) + 1)
 
#define PI 3.14159
 
/*
 * used for heapsort
 */
void swap(double *x, double *y)
{
    double temp;
    
    temp = *x;
    *x = *y;
    *y = temp;
}
 
void swap(int *x, int *y)
{
    int temp;
    
    temp = *x;
    *x = *y;
    *y = temp;
}
 
/*
 * maxheapify: used to let the value at a[i] "double down" in the
 * max-heap so that the subtree rooted at index i becomes a max-heap.
 */
void maxheapify(double a[], int i, int heapsize)
{
    int largest, left, right;
    
    left = Left(i);
    right = Right(i);
    
    if (left <= heapsize && a[left] > a[i])
        largest = left;
    else
        largest = i;
    if (right <= heapsize && a[right] > a[largest])
        largest = right;
    
    if (largest != i)
    {
        swap(&a[i], &a[largest]);
        maxheapify(a, largest, heapsize);    /* recurrsion */
    }
    
    return;        /* return when largest == i */
}
 
/*
 * heapsort: contains two procedures, one is to build the max-heap,
 * the other is to delete max to gain the sorted array.
 */
 
void heapsort(double a[], int heapsize)
{
    int i;
    
    for (i = heapsize / 2; i >= 1; i--)        /* build max-heap */
        maxheapify(a, i, heapsize);
    
    for (i = heapsize; i >= 2; i--)            /* delete max */
    {
        swap(&a[1], &a[i]);
        heapsize--;
        maxheapify(a, 1, heapsize);
    }
    
}
 
/*
 * used for solving isotropic TV. Sort c while sorting a,,,,
 */
void maxheapify(double a[], int i, int heapsize, double *c)
{
    int largest, left, right;
    
    left = Left(i);
    right = Right(i);
    
    if (left <= heapsize && a[left] > a[i])
        largest = left;
    else
        largest = i;
    if (right <= heapsize && a[right] > a[largest])
        largest = right;
    
    if (largest != i)
    {
        swap(&a[i], &a[largest]);
        swap(&c[i], &c[largest]);
        maxheapify(a, largest, heapsize);    /* recurrsion */
    }
    
    return;        /* return when largest == i */
}
 
void heapsort(double *a, int heapsize, double *c)
{
    int i;
    
    for (i = heapsize / 2; i >= 1; i--)        /* build max-heap */
        maxheapify(a, i, heapsize, c);
    
    for (i = heapsize; i >= 2; i--)            /* delete max */
    {
        swap(&a[1], &a[i]);
        swap(&c[1], &c[i]);
        heapsize--;
        maxheapify(a, 1, heapsize, c);
    }
    
}
 
/*
 *bubble sort.. a[0] is not used, and n=length(a)-1
 */
void bubblesort(double a[], int n)
{
    int i,j;
    
    for(j=1;j<n;++j)
    {
        for(i=j;i<=n;++i)
        {
            if(a[j]>a[i])
            {
                swap(&a[i],&a[j]);
            }
        }
    }
    
}
 
void bubblesort(double a[], int n, double *c)
{
    int i,j;
    
    for(j=1;j<n;++j)
    {
        for(i=j;i<=n;++i)
        {
            if(a[j]>a[i])
            {
                swap(&a[i],&a[j]);
                swap(&c[i],&c[j]);
            }
        }
    }
    
}
 
void bubblesort(double a[], int n, int *c)
{
    int i,j;
    
    for(j=1;j<n;++j)
    {
        for(i=j;i<=n;++i)
        {
            if(a[j]>a[i])
            {
                swap(&a[i],&a[j]);
                swap(&c[i],&c[j]);
            }
        }
    }
    
}
 
 
 
 
/*
 * @function randindex
 * @remark ?：??：?????[low,high]???????：???
 */
int randindex ( int low, int high )
{
    return rand () % ( high - low ) + low;
}
 
 
void quicksort ( double * array_to_sort, int low, int high )
{
    if ( low >= high ) return;
 
    /* 
     * ?? [l + 1,h] ???????????????：????? index??
     * ???????????：＜???????????????????：? index??????????????????
     * ????????????：＜ pivot ??????????：????：＜?????????????????????????：???????：＜
     */
    int index = randindex ( low + 1, low );
    swap ( array_to_sort + low, array_to_sort + index ); 
 
    int l = low, h = high;
    double pivot = array_to_sort[l];
 
    while ( l < h )
    {
        while ( l < h && array_to_sort[h] >= pivot ) h--;
        array_to_sort[l] = array_to_sort[h];
        while ( l < h && array_to_sort[l] <= pivot ) l++;
        array_to_sort[h] = array_to_sort[l];
    }
    array_to_sort[l] = pivot;
 
    quicksort ( array_to_sort, low, l - 1 );
    quicksort ( array_to_sort, l + 1, high );
}
 
void quicksort ( double * array_to_sort, int low, int high , double *c)
{
    if ( low >= high ) return;
 
    /* 
     * ?? [l + 1,h] ???????????????：????? index??
     * ???????????：＜???????????????????：? index??????????????????
     * ????????????：＜ pivot ??????????：????：＜?????????????????????????：???????：＜
     */
    int index = randindex ( low + 1, low );
    swap ( array_to_sort + low, array_to_sort + index ); 
    swap(c+low,c+index);
 
    int l = low, h = high;
    double pivot = array_to_sort[l];
    double tt=c[l];
 
    while ( l < h )
    {
        while ( l < h && array_to_sort[h] >= pivot ) h--;
        array_to_sort[l] = array_to_sort[h];
        c[l]=c[h];
        while ( l < h && array_to_sort[l] <= pivot ) l++;
        array_to_sort[h] = array_to_sort[l];
        c[h]=c[l];
    }
    array_to_sort[l] = pivot;
    c[l]=tt;
 
    quicksort ( array_to_sort, low, l - 1 ,c);
    quicksort ( array_to_sort, l + 1, high ,c);
}
 
 
 
/*
 *handle the overflow of matrix with the periodic boundary condition.
 */
void SwSubscript(int i,int j,int m,int n,int *dsti,int *dstj)
{//
    *dsti=i;
    *dstj=j;
    
    if(i<0)
        *dsti=m+i;
    else if(i>=m)
        *dsti=i-m;
    if(j<0)
        *dstj=n+j;
    else if(j>=n)
        *dstj=j-n;
}
 
/*
 *used to solve elementary subproblem with the form  0.5*(x-beta)^2+miu*||x-e||_l1
 *the input argument, e is a n+1 vector, e[0] is not used for the function 'heapsort'
 */
double SolveEleSub(double beta,double miu,double *e,int n)
{
    double x,x1,x2;
    int i;
    
    //sort vector e
    bubblesort(e,n);//e[0] is not used
//     quicksort(e,1,n);
    
    
    //start the solving procedure
    // when x < e(1)
    x=n*miu+beta;
    if(x<e[1])
        return x;
    
// when x > e(n)
    x=-n*miu+beta;
    if( x>e[n])
        return x;
    
// when e(i)<x<e(i+1)
    for (i=1;i<(n);++i)
    {
        //x=-(i-(d-i))*miu/2+beta;
        x=(n-2*i)*miu+beta;
        if( x>e[i] && x<e[i+1])
            return x;
    }
    
// when x=e(i)
    for( i=1;i<=n;++i)
    {
        x=e[i];
        // x1=-((i-1)-(d-i)+1)*miu/2+beta;
        // x2=-((i-1)-(d-i)-1)*miu/2+beta;
        x1=(n-2*i)*miu+beta;
        x2=(n-2*i+2)*miu+beta;
        if( (x>x1||fabs(x-x1)<1e-15) && (x<x2||fabs(x-x2)<1e-15))
//             if( ((x-x1)>=0) && ((x2-x)<=0))
            return x;
    }
}
 
double SolveEleSub(double beta,double miu,double *e,int n, int * s)
{
    double x,x1,x2;
    int i;
    
    //sort vector e
    bubblesort(e,n);//e[0] is not used
    
    //start the solving procedure
    // when x < e(1)
    x=n*miu+beta;
    if(x<e[s[1]])
        return x;
    
// when x > e(n)
    x=-n*miu+beta;
    if( x>e[s[n]])
        return x;
    
// when e(i)<x<e(i+1)
    for (i=1;i<(n);++i)
    {
        //x=-(i-(d-i))*miu/2+beta;
        x=(n-2*i)*miu+beta;
        if( x>e[s[i]] && x<e[s[i+1]])
            return x;
    }
    
// when x=e(i)
    for( i=1;i<=n;++i)
    {
        x=e[s[i]];
        // x1=-((i-1)-(d-i)+1)*miu/2+beta;
        // x2=-((i-1)-(d-i)-1)*miu/2+beta;
        x1=(n-2*i)*miu+beta;
        x2=(n-2*i+2)*miu+beta;
        if( (x>x1||fabs(x-x1)<1e-15) && (x<x2||fabs(x-x2)<1e-15))
            return x;
    }
}
/*
 *used to solve elementary subproblem with the form  0.5*(x-beta)^2+miu*sum_i||ai*x-ei||_l1
 *the input argument, e is a n+1 vector, e[0] is not used, i.e., n=length(e)-1
 */
double SolveEleSub(double beta,double miu,double *e,int n, double *a, double *sa)
{
    double x,x1,x2;
    int i;
    
//     sa=(double*)malloc((n+1)*sizeof(double));
    
    sa[0]=0;
    for(i=1;i<=n;++i)
    {
        e[i]/=a[i];
        sa[0]-=a[i];
    }
    sa[n]=-sa[0];
    
    //sort vector e/a and a
    bubblesort(e,n,a);//e[0] is not used
//     quicksort(e,1,n,a);
    
    //start the solving procedure
    // when x < e(1)
//     x=n*miu+beta;
    x=beta-sa[0]*miu;
    if(x<e[1])
    {
//         free(sa);
        return x;
    }
    
// when x > e(n)
    x=beta-sa[n]*miu;
    if(x>e[n])
    {
//         free(sa);
        return x;
    }
    
// when e(i)<x<e(i+1)
    for (i=1;i<n;++i)
    {
        sa[i]=sa[i-1]+2*a[i]; //precompute -sum_i{a}+sum_(n-i){a}
        
        //x=-(i-(d-i))*miu/2+beta;
//         x=(n-2*i)*miu+beta;
        x=beta-sa[i]*miu;
        
        if( x>e[i] && x<e[i+1])
        {
//             free(sa);
            return x;
        }
    }
    
// when x=e(i)
    for(i=1;i<=n;++i)
    {
        x=e[i];
        // x1=-((i-1)-(d-i)+1)*miu/2+beta;
        // x2=-((i-1)-(d-i)-1)*miu/2+beta;
        x1=beta-sa[i]*miu;
        x2=beta-sa[i-1]*miu;
        if( (x>x1||fabs(x-x1)<1e-15) && (x<x2||fabs(x-x2)<1e-15))// x1 <= x <= x2
        {
//             free(sa);
            return x;
        }
    }
    
}
 
 
/*
 *used to solve elementary subproblem with the form  sum_i||ai*x-ei||_l1
 *the input argument, e is a n+1 vector, e[0] is not used, i.e., n=length(e)-1
 */
double SolveEleSubM(double *e,int n, double *a, double *sa)
{
    double x,x1,x2;
    int i;
    
//     double *sa=(double*)malloc((n+1)*sizeof(double));
    
    sa[0]=0;
    for(i=1;i<=n;++i)
    {
//         e[i]/=a[i];
        sa[0]-=a[i];
    }
//     sa[n]=-sa[0];
    
    //sort vector e/a and a
    bubblesort(e,n,a);//e[0] is not used
    
    //start the solving procedure
    for (i=1;i<n;++i)
    {
        sa[i]=sa[i-1]+2*a[i]; //precompute -sum_i{a}+sum_(n-i){a}
//         mexPrintf("%d\n",i);
        if( !(sa[i]<0))
        {
            return e[i];
        }
    }
    
    return e[n];
    
}
 
 
 
 
 
 
/*
 * compute Local Adjustment Table (LAT) in advance
 *
 */
void LAT(double *A,int d, double *adpcoef)
{
    int k,l,r,c,m;
    double coef_x;
    
    m=2*d-1;
    
//     ii=0;//calculate LAT with 4 patrs.
    for(k=-d+1;k<=0;++k)//left up
    {
        for(l=-d+1;l<=0;++l)
        {
            coef_x=0;
            for(r=0;r<d+k;++r)
            {
                for(c=0;c<d+l;++c)
                {
                    coef_x+=A[(d-1-c)*d+(d-1-r)]*A[(d-1+l-c)*d+(d-1+k-r)];
                }
            }
//             if(k==0 && l==0)
//             {coef_x=0;}
            adpcoef[(d-1+l)*m+(d-1+k)]=coef_x;
        }
    }
    adpcoef[(d-1)*m+d-1]=0;
    
//     for(k=-d+1;k<=0;++k)//right up
//     {
//         for(l=1;l<=d-1;++l)
//         {
//             coef_x=0;
//             for(r=0;r<d+k;++r)
//             {
//                 for(c=0;c<d-l;++c)
//                 {
//                     coef_x+=A[(c)*d+(d-1-r)]*A[(c)*d+(r)];
//                 }
//             }
//
//             adpcoef[(d-1+l)*(2*d-1)+(d-1+k)]=coef_x;
//         }
//     }
    
    for(k=0;k<d;++k)//right up
    {
        for(l=d;l<m;++l)
        {
            adpcoef[l*m+k]=adpcoef[(d-2-(l-d))*m+k];
        }
    }
    
    for(k=d;k<m;++k)//left down
    {
        for(l=0;l<d;++l)
        {
            adpcoef[l*m+k]=adpcoef[l*m+(d-2-(k-d))];
        }
    }
    
    for(k=d;k<m;++k)//right down
    {
        for(l=d;l<m;++l)
        {
            adpcoef[l*m+k]=adpcoef[(d-2-l+d)*m+(d-2-(k-d))];
        }
    }
    
    
}
 
/*
 * Generate random integer i with probabilities Prob[i=k] = L[k],
 * with complexity O(log n)
 **/
double ** genP(double *L, int n, int m)
{//m=log2(n)
    
    int i,k;
    double **P;
    
    P=(double **)malloc((m+1)*sizeof(double *));
    P[0]=L;
    for(k=1; k<=m; ++k)
    {
        P[k]=(double*)malloc((int)pow(2,double(m-k))*sizeof(double));
    }
    
    for(k=1; k<=m; ++k)
    {
        for(i=0; i<pow(2,(double)(m-k)); ++i)
        {
            P[k][i]=P[k-1][2*i]+P[k-1][2*i+1];
        }
    }
    
    return P;
}
 
void modifyP(double **P, int m, double a,int i)
{
    //modify P, position: i, value: a
    int k,pos=i,prepos;
    double minu=a-P[0][i];
    
    P[0][i]=a;
    for(k=1;k<=m;++k)
    {
        pos/=2;
        P[k][pos]+=minu;
    }
}
 
int randCounter(double **P,int m)
{
    int i,k;
    
//     srand(time(NULL));
    
    i=0;
    
    for(k=m; k>=1; --k)
    {
        if(rand()%1000 < 1000*(P[k-1][2*i]/P[k][i]))
            i=2*i;
        else
            i=2*i+1;
    }
    
    return i;
}
 
void freeP(double **P, int m)
{
    for(int i=1;i<=m;++i)
        free(P[i]);
    
    free(P);
}
 
 
 


