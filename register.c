#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
void dge0(double *A,double *B,double *C,int n)
{
    int i,j,k;
    for(i=0;i<n;i++)
    {
        for(j=0;j<n;j++) 
        {
            for(k=0;k<n;k++)
            {
                C[i*n+j]+=A[i*n+k]*B[k*n+j];
            }
        }
    }
}

void dge1(double *A,double *B,double *C,int n)
{
    int i,j,k;
    for (i=0; i<n; i++) 
    {
        for (j=0; j<n; j++) 
        {
            register double r = C[i*n+j] ;
            for (k=0; k<n; k++)
            {
                r += A[i*n+k] * B[k*n+j];
                C[i*n+j] = r;
            }
        }
    }
}

void dge2(double *a, double *b, double *c,int n)
{
    int i, j, k;
    for (i = 0; i < n; i+=2)
    {
	for (j = 0; j < n; j+=2)
        {
             for (k = 0; k < n; k+=2)
             {
	         c[i*n + j]         = a[i*n + k]*b[k*n + j] + a[i*n + k+1]*b[(k+1)*n + j] + c[i*n + j];
                 c[(i+1)*n + j]     = a[(i+1)*n + k]*b[k*n + j] + a[(i+1)*n + k+1]*b[(k+1)*n + j] + c[(i+1)*n + j];
                 c[i*n + (j+1)]     = a[i*n + k]*b[k*n + (j+1)] + a[i*n + k+1]*b[(k+1)*n + (j+1)] + c[i*n + (j+1)];
                 c[(i+1)*n + (j+1)] = a[(i+1)*n + k]*b[k*n + (j+1)] + a[(i+1)*n + k+1]*b[(k+1)*n + (j+1)] 
                                       + c[(i+1)*n+ (j+1)];
             }
        }
    }
}

void exa1(double *a, double *b,double *c,int n)
{
    int i,j,k;
    for(i = 0; i < n; i += 2)
       for(j = 0; j < n; j += 2)  {
            register int t = i*n+j; register int tt = t+n; 
            register double c00 = c[t]; register double c01 = c[t+1];  
            register double c10 = c[tt]; register double c11 = c[tt+1];

            for(k = 0; k < n; k += 2) {
                register int ta = i*n+k; register int tta = ta+n; register int tb = k*n+j; register int ttb = tb+n;
                register double a00 = a[ta]; register double a10 = a[tta]; 
                register double b00 = b[tb]; register double b01 = b[tb+1]; 

                c00 += a00*b00 ; c01 += a00*b01 ; c10 += a10*b00 ; c11 += a10*b01 ;

                a00 = a[ta+1]; a10 = a[tta+1]; b00 = b[ttb]; b01 = b[ttb+1];

                c00 += a00*b00 ; c01 += a00*b01 ; c10 += a10*b00 ; c11 += a10*b01 ;
             }
             c[t] = c00;
             c[t+1] = c01;
             c[tt] = c10;
             c[tt+1] = c11;
        }
}

void exa2(double *a, double *b,double *c,int n)
{
    int i,j,k;
    for(i = 0; i < n; i += 2)
    {
       for(j = 0; j < n; j += 2)  
       {
            register int t = i*n+j; register int tt = t+n; 
            register double c00 = c[t]; register double c01 = c[t+1];  
            register double c10 = c[tt]; register double c11 = c[tt+1];
            
            for(k = 0; k < n; k += 2) 
            {
                register int ta = i*n+k; register int tta = ta+n; register int tb = k*n+j; register int ttb = tb+n;
                register double a00 = a[ta]; register double a01 = a[ta+1]; 
                register double a10 = a[tta]; register double a11 = a[tta+1];
                register double b00 = b[tb]; register double b01 = b[tb+1]; 
                register double b10 = b[ttb]; register double b11 = b[ttb+1];
                c00 += a00*b00 + a01*b10;
                c01 += a00*b01 + a01*b11;
                c10 += a10*b00 + a11*b10;
                c11 += a10*b01 + a11*b11;
             }

             c[t] = c00;
             c[t+1] = c01;
             c[tt] = c10;
             c[tt+1] = c11;
             
        }
    }
}

void dge3(double *a, double *b,double *c,int n)
{
    int i,j,k;
    for(i = 0; i < n; i += 4)
    {
       for(j = 0; j < n; j += 4)  
       {
            
            for(k = 0; k < n; k += 4) 
            {
                register int t1 = i*n+j; register int t2 = (i+1)*n+j; register int t3 = (i+2)*n+j; register int t4 = (i+3)*n+j; 
                

                register double c00 = c[t1]; register double c01 = c[t2];  
                register double c02 = c[t3]; register double c03 = c[t4];
                register double c04 = c[t1+1]; register double c05 = c[t2+1];  
                register double c06 = c[t3+1]; register double c07 = c[t4+1];

                register int ta = i*n+k; register int tb = ta+n; register int tc = ta+2*n; register int td = ta+3*n;
                register int te = k*n+j; register int tf = te+n; register int tg = te+2*n; register int th = te+3*n;
                register double a00 = a[ta]; register double a01 = a[tb]; 
                register double a02 = a[tc]; register double a03 = a[td];
                register double b00 = b[te]; register double b01 = b[te+1]; 
                c00 += a00*b00; c01 += a01*b00; c02 += a02*b00; c03 += a03*b00;
                c04 += a00*b01; c05 += a01*b01; c06 += a02*b01; c07 += a03*b01;

                a00 = a[ta+1]; a01 = a[tb+1]; a02 = a[tc+1]; a03 = a[td+1];
                b00 = b[tf]; b01 = b[tf+1];
                c00 += a00*b00; c01 += a01*b00; c02 += a02*b00; c03 += a03*b00;
                c04 += a00*b01; c05 += a01*b01; c06 += a02*b01; c07 += a03*b01;

                a00 = a[ta+2]; a01 = a[tb+2]; a02 = a[tc+2]; a03 = a[td+2];
                b00 = b[tg]; b01 = b[tg+1];
                c00 += a00*b00; c01 += a01*b00; c02 += a02*b00; c03 += a03*b00;
                c04 += a00*b01; c05 += a01*b01; c06 += a02*b01; c07 += a03*b01;


                a00 = a[ta+3]; a01 = a[tb+3]; a02 = a[tc+3]; a03 = a[td+3];
                b00 = b[th]; b01 = b[th+1];
                c00 += a00*b00; c01 += a01*b00; c02 += a02*b00; c03 += a03*b00;
                c04 += a00*b01; c05 += a01*b01; c06 += a02*b01; c07 += a03*b01;

                c[t1] = c00;   c[t2] = c01;   c[t3] = c02;   c[t4] = c03;
                c[t1+1] = c04; c[t2+1] = c05; c[t3+1] = c06; c[t4+1] = c07;

                c00 = c[t1+2]; c01 = c[t2+2];  
                c02 = c[t3+2]; c03 = c[t4+2];
                c04 = c[t1+3]; c05 = c[t2+3];  
                c06 = c[t3+3]; c07 = c[t4+3];

                a00 = a[ta]; a01 = a[tb]; a02 = a[tc]; a03 = a[td];
                b00 = b[te+2]; b01 = b[te+3];
                c00 += a00*b00; c01 += a01*b00; c02 += a02*b00; c03 += a03*b00;
                c04 += a00*b01; c05 += a01*b01; c06 += a02*b01; c07 += a03*b01;

                a00 = a[ta+1]; a01 = a[tb+1]; a02 = a[tc+1]; a03 = a[td+1];
                b00 = b[tf+2]; b01 = b[tf+3];
                c00 += a00*b00; c01 += a01*b00; c02 += a02*b00; c03 += a03*b00;
                c04 += a00*b01; c05 += a01*b01; c06 += a02*b01; c07 += a03*b01;
                
                a00 = a[ta+2]; a01 = a[tb+2]; a02 = a[tc+2]; a03 = a[td+2];
                b00 = b[tg+2]; b01 = b[tg+3];
                c00 += a00*b00; c01 += a01*b00; c02 += a02*b00; c03 += a03*b00;
                c04 += a00*b01; c05 += a01*b01; c06 += a02*b01; c07 += a03*b01; 
                
                a00 = a[ta+3]; a01 = a[tb+3]; a02 = a[tc+3]; a03 = a[td+3];
                b00 = b[th+2]; b01 = b[th+3];
                c00 += a00*b00; c01 += a01*b00; c02 += a02*b00; c03 += a03*b00;
                c04 += a00*b01; c05 += a01*b01; c06 += a02*b01; c07 += a03*b01;
                
                c[t1+2] = c00; c[t2+2] = c01; c[t3+2] = c02; c[t4+2] = c03;
                c[t1+3] = c04; c[t2+3] = c05; c[t3+3] = c06; c[t4+3] = c07;
             }

                 
        }
    }
}

double calG(int n,float t)
{
    double g;
    g = 2 * pow(n,3) / t * pow(10,-9);
    return g;
}

double differ(double *c,double *d,int n)
{
    double max_dft=0;
    int i;
    for(i = 0; i < n * n; i++){
	if(fabs(c[i] - d[i]) > max_dft)
	    max_dft = fabs(c[i] - d[i]);

	}
    return max_dft;
}

int main()
{
    float t;
    struct timespec t1, t2;

    printf ("Calculating...\n");
    srand(time(NULL));
    int i,j,k,m=2048;
    double *A,*B,*C1,*C2,*D1,*D2,*D3,*D4,*D5;

    A=(double *)malloc(m*m*sizeof(double *));
    B=(double *)malloc(m*m*sizeof(double *));
    
    for(i=0;i<m;i++)
    {
        for(j=0;j<m;j++)
        {
            A[i*m+j]=rand()/(double)(RAND_MAX/100);
            B[i*m+j]=rand()/(double)(RAND_MAX/100);
         }
    }

    int n;
    for(n=64;n<=2048;n*=2)
    {
        printf("When n=%d:\n",n);

        C1=(double *)malloc(n*n*sizeof(double *));
        C2=(double *)malloc(n*n*sizeof(double *));

        clock_gettime(CLOCK_MONOTONIC, &t1);
        dge0(A,B,C1,n);
        clock_gettime(CLOCK_MONOTONIC, &t2);
        t = t2.tv_sec - t1.tv_sec + pow(10,-9)*(t2.tv_nsec - t1.tv_nsec);
        printf("The time of NO is %fs\n",t);
        printf("NO:%f\n",calG(n,t));

        clock_gettime(CLOCK_MONOTONIC, &t1);
        dge1(A,B,C2,n);
        clock_gettime(CLOCK_MONOTONIC, &t2);
        t = t2.tv_sec - t1.tv_sec + pow(10,-9)*(t2.tv_nsec - t1.tv_nsec);
        printf("The time of RE is %fs\n",t);
        printf("RE:%f\n",calG(n,t));
        printf("The difference from C1 is %lf\n",differ(C1,C2,n));

        D1=(double *)malloc(n*n*sizeof(double *));
        D2=(double *)malloc(n*n*sizeof(double *));
        D3=(double *)malloc(n*n*sizeof(double *));
        D4=(double *)malloc(n*n*sizeof(double *));
        D5=(double *)malloc(n*n*sizeof(double *));

        clock_gettime(CLOCK_MONOTONIC, &t1);
        dge2(A,B,D1,n);
        clock_gettime(CLOCK_MONOTONIC, &t2);
        t = t2.tv_sec - t1.tv_sec + pow(10,-9)*(t2.tv_nsec - t1.tv_nsec);
        printf("The time of 2*2NO is %fs\n",t);
        printf("N=%d(2*2 NO):%f\n",n,calG(n,t));
        printf("The difference from C1 is %lf\n",differ(C1,D1,n));

        clock_gettime(CLOCK_MONOTONIC, &t1);
        exa1(A,B,D2,n);
        clock_gettime(CLOCK_MONOTONIC, &t2);
        t = t2.tv_sec - t1.tv_sec + pow(10,-9)*(t2.tv_nsec - t1.tv_nsec);
        printf("The time of 2*2with8 is %fs\n",t);
        printf("N=%d(2*2 with8):%f\n",n,calG(n,t));
        printf("The difference from C1 is %lf\n",differ(C1,D2,n));

        clock_gettime(CLOCK_MONOTONIC, &t1);
        exa2(A,B,D3,n);
        clock_gettime(CLOCK_MONOTONIC, &t2);
        t = t2.tv_sec - t1.tv_sec + pow(10,-9)*(t2.tv_nsec - t1.tv_nsec);
        printf("The time of 2*2with12 is %fs\n",t);
        printf("N=%d(2*2 with12):%f\n",n,calG(n,t));
        printf("The difference from C1 is %lf\n",differ(C1,D3,n));

        clock_gettime(CLOCK_MONOTONIC, &t1);
        dge3(A,B,D4,n);
        clock_gettime(CLOCK_MONOTONIC, &t2);
        t = t2.tv_sec - t1.tv_sec + pow(10,-9)*(t2.tv_nsec - t1.tv_nsec);
        printf("The time of 4*4with14 is %fs\n",t);
        printf("N=%d(4*4 with12):%f\n",n,calG(n,t));
        printf("The difference from C1 is %lf\n",differ(C1,D4,n));



        printf("----------------------------------------------------\n");
    }

    return 0;

}
# CS211
