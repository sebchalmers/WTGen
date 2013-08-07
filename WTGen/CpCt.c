/*
 
 By Sebastien Gros
 Assistant Professor
 
 Department of Signals and Systems
 Chalmers University of Technology
 SE-412 96 Göteborg, SWEDEN
 
 Compute the evaluation of a surface B-spline (order 3), using the cox-deBoor/Böhm formula
 Provides 1st and second-order derivatives (function EvalSpline) or 1st and 2nd (function EvalSpline2)
 The blending functions are code-generated in Blender.h (for-looped version commented in this code, can be used with splines of different order than 3)
 
 */


#include <stdio.h>
  

int findspan(float xi, const float knots[], const int lenght_knots)
{
    // Binary search over the knots
    int index_up   =  lenght_knots;
    int index_low  =  0;
    int index_middle;

    while (index_up - index_low > 1)
    {
        index_middle = ((index_up+index_low)/2);
        if (xi < knots[index_middle])
        {
            index_up = index_middle;
        }
        else
        {
            index_low = index_middle;
        }
    }
    return index_low;
}

// Unrolled blending function
#include "Blender.h"

// Original blending function
/*
 float Pijeval(const float mat[],int i,int j,const int n)
 {
 // Line major "matrix" P -> take out i,j entry
 
 //printf("index = %d\n", i*n + j);
 return mat[i*n + j];
 }
 
 
 float Blend(float x_basis[], float y_basis[], int i_x, int i_y, const float P[], const int n, int length_x_basis, int length_y_basis)
 {
 //printf("Blending \n");
 float S = 0;
 float x_basis_k1;
 for (int k1 = 0; k1 < length_x_basis; k1++)
 {
 x_basis_k1 = x_basis[k1];
 for (int k2 = 0; k2 < length_y_basis; k2++)
 {
 S += Pijeval(P,i_x-k1,i_y-k2,n)*x_basis_k1*y_basis[k2];
 }
 }
 
 return S;
 }
*/

//#include "BasisFunc.h"




int basisFuncs(float basis[], float xi, const int order, const float U[], int i)
{
    //Minimal implementation of the Cox-deBoor formula

    
    // UNROLLING THESE LOOPS IMPROVES THE SPEED BY ONLY A FEW %
    
    int iplus = i+1;
    int pminus;
    int ipluspminusk;
    
    float Uiplus = U[iplus];
    float Ui     = U[i];
    float Den    = Uiplus - U[i];
    
    //Compute the first step (special branch, p=1, Ni,0 = 1)
    basis[0] = (xi     - Ui) / Den;
    basis[1] = (Uiplus -   xi) / Den;
    
    
    //Clear out the remaining values of basis
    for (int k = 2; k < order+1; k++)
    {
        basis[k] = 0;
    }
    
    //These loops could be unrolled into two different functions, depending on the order (2 or 3)
    for (int p = 2; p < order+1; p++)
    {
        pminus   = p-1;
        
        //Update of Ni-p,p is standalone (cross arrow):
        Den  = Uiplus - U[i-pminus];
        basis[p] = (Uiplus - xi)*basis[pminus] / Den;

        for (int k=p-1; k > 0; k--)
        {
            //Flat arrow
            basis[k]  = (xi - U[i-k])*basis[k] / Den;
            
            //Cross arrow
            ipluspminusk = iplus+p-k;
            Den   =  U[ipluspminusk] - U[iplus-k];
            basis[k] += (U[ipluspminusk] - xi)*basis[k-1] / Den;

        }
        // Update of Ni,p is standalone (flat arrow):
        basis[0] = (xi - Ui)*basis[0] / Den;
        
    }
    return 0;
}



int Cp0(float x, float y, float out[])
{
    
    #include "Cp.h"
    
    //interpolation points >= 0
    x -= x_shift;
    y -= y_shift;
    
    // Eval Spline
    int ix = findspan(x, knots_x, length_knots_x);
    int iy = findspan(y, knots_y, length_knots_y);
    
    float basis_x[p+1];
    float basis_y[q+1];
    
    basisFuncs(basis_x, x, p, knots_x, ix);
    basisFuncs(basis_y, y, q, knots_y, iy);
    
    out[0] = Blend44(basis_x, basis_y, ix, iy, P, n);
    
    return 0;
}

int Cp1(float x, float y, float out[])
{
    
    #include "Cp.h"
    
    //interpolation points >= 0
    x -= x_shift;
    y -= y_shift;
    
    // Eval Spline
    int ix = findspan(x, knots_x, length_knots_x);
    int iy = findspan(y, knots_y, length_knots_y);
     
    float basis_x[p+1];      
    float basis_y[q+1];
    
    basisFuncs(basis_x, x, p, knots_x, ix);
    basisFuncs(basis_y, y, q, knots_y, iy);
    
    out[0] = Blend44(basis_x, basis_y, ix, iy, P, n);

    // Eval derivatives
    int ix_tilde = findspan(x, Ux, length_Ux);
    int iy_tilde = findspan(y, Uy, length_Uy);
    
    float basis_x_tilde[p];
    float basis_y_tilde[q];
    
    basisFuncs(basis_x_tilde, x, p-1, Ux, ix_tilde);
    basisFuncs(basis_y_tilde, y, q-1, Uy, iy_tilde);

    out[1] = Blend34(basis_x_tilde, basis_y,       ix_tilde, iy,       Px, n   );
    out[2] = Blend43(basis_x      , basis_y_tilde, ix,       iy_tilde ,Py, n-1 );
    
    return 0;
}


int Cp2(float x, float y, float out[])
{
    
    #include "Cp.h"
    
    //interpolation points >= 0
    x -= x_shift;
    y -= y_shift;
    
    // Eval Spline
    int ix = findspan(x, knots_x, length_knots_x);
    int iy = findspan(y, knots_y, length_knots_y);
    
    float basis_x[p+1];
    float basis_y[q+1];
    
    basisFuncs(basis_x, x, p, knots_x, ix);
    basisFuncs(basis_y, y, q, knots_y, iy);
    
    out[0] = Blend44(basis_x, basis_y, ix, iy, P, n);
    
    // Eval derivatives
    int ix_tilde = findspan(x, Ux, length_Ux);
    int iy_tilde = findspan(y, Uy, length_Uy);
    
    float basis_x_tilde[p];
    float basis_y_tilde[q];
    
    basisFuncs(basis_x_tilde, x, p-1, Ux, ix_tilde);
    basisFuncs(basis_y_tilde, y, q-1, Uy, iy_tilde);
    
    out[1] = Blend34(basis_x_tilde, basis_y,       ix_tilde, iy,       Px, n   );
    out[2] = Blend43(basis_x      , basis_y_tilde, ix,       iy_tilde ,Py, n-1 );
    
    //Eval curvature
    int ixx = findspan(x, Uxx, length_Uxx);
    int iyy = findspan(y, Uyy, length_Uyy);
     
    float basis_xx[p-1];
    float basis_yy[q-1];
     
    basisFuncs(basis_xx,   x,   p-2,   Uxx,    ixx);
    basisFuncs(basis_yy,   y,   q-2,   Uyy,    iyy);
     
    out[3] = Blend24(basis_xx,      basis_y,            ixx,       iy,  Pxx,    n );
    out[4] = Blend42(basis_x,       basis_yy,            ix,      iyy,  Pyy,  n-2 );
    out[5] = Blend33(basis_x_tilde, basis_y_tilde, ix_tilde, iy_tilde,  Pxy,  n-1 );
    
    return 0;
}


int Ct0(float x, float y, float out[])
{
    
    #include "Ct.h"
    
    //interpolation points >= 0
    x -= x_shift;
    y -= y_shift;
    
    // Eval Spline
    int ix = findspan(x, knots_x, length_knots_x);
    int iy = findspan(y, knots_y, length_knots_y);
    
    float basis_x[p+1];
    float basis_y[q+1];
    
    basisFuncs(basis_x, x, p, knots_x, ix);
    basisFuncs(basis_y, y, q, knots_y, iy);
    
    out[0] = Blend44(basis_x, basis_y, ix, iy, P, n);
    
    return 0;
}

int Ct1(float x, float y, float out[])
{
    
    #include "Ct.h"
    
    //interpolation points >= 0
    x -= x_shift;
    y -= y_shift;
    
    // Eval Spline
    int ix = findspan(x, knots_x, length_knots_x);
    int iy = findspan(y, knots_y, length_knots_y);
    
    float basis_x[p+1];
    float basis_y[q+1];
    
    basisFuncs(basis_x, x, p, knots_x, ix);
    basisFuncs(basis_y, y, q, knots_y, iy);
    
    out[0] = Blend44(basis_x, basis_y, ix, iy, P, n);
    
    // Eval derivatives
    int ix_tilde = findspan(x, Ux, length_Ux);
    int iy_tilde = findspan(y, Uy, length_Uy);
    
    float basis_x_tilde[p];
    float basis_y_tilde[q];
    
    basisFuncs(basis_x_tilde, x, p-1, Ux, ix_tilde);
    basisFuncs(basis_y_tilde, y, q-1, Uy, iy_tilde);
    
    out[1] = Blend34(basis_x_tilde, basis_y,       ix_tilde, iy,       Px, n   );
    out[2] = Blend43(basis_x      , basis_y_tilde, ix,       iy_tilde ,Py, n-1 );
    
    return 0;
}


int Ct2(float x, float y, float out[])
{
    
    #include "Ct.h"
    
    //interpolation points >= 0
    x -= x_shift;
    y -= y_shift;
    
    // Eval Spline
    int ix = findspan(x, knots_x, length_knots_x);
    int iy = findspan(y, knots_y, length_knots_y);
    
    float basis_x[p+1];
    float basis_y[q+1];
    
    basisFuncs(basis_x, x, p, knots_x, ix);
    basisFuncs(basis_y, y, q, knots_y, iy);
    
    out[0] = Blend44(basis_x, basis_y, ix, iy, P, n);
    
    // Eval derivatives
    int ix_tilde = findspan(x, Ux, length_Ux);
    int iy_tilde = findspan(y, Uy, length_Uy);
    
    float basis_x_tilde[p];
    float basis_y_tilde[q];
    
    basisFuncs(basis_x_tilde, x, p-1, Ux, ix_tilde);
    basisFuncs(basis_y_tilde, y, q-1, Uy, iy_tilde);
    
    out[1] = Blend34(basis_x_tilde, basis_y,       ix_tilde, iy,       Px, n   );
    out[2] = Blend43(basis_x      , basis_y_tilde, ix,       iy_tilde ,Py, n-1 );
    
    //Eval curvature
    int ixx = findspan(x, Uxx, length_Uxx);
    int iyy = findspan(y, Uyy, length_Uyy);
    
    float basis_xx[p-1];
    float basis_yy[q-1];
    
    basisFuncs(basis_xx,   x,   p-2,   Uxx,    ixx);
    basisFuncs(basis_yy,   y,   q-2,   Uyy,    iyy);
    
    out[3] = Blend24(basis_xx,      basis_y,            ixx,       iy,  Pxx,    n );
    out[4] = Blend42(basis_x,       basis_yy,            ix,      iyy,  Pyy,  n-2 );
    out[5] = Blend33(basis_x_tilde, basis_y_tilde, ix_tilde, iy_tilde,  Pxy,  n-1 );
    
    return 0;
}

/*
int main()
{
    float x = 5.;
    float y = 3.;
    
    float out[6];
    
    // Cp table
    for (int k=0;k<1e6;k++)
    {
        //printf("k = %d",k);
        Cp0(x,y,out);
    }

    for (int k=0;k<1e6;k++)
    {
        //printf("k = %d",k);
        Cp1(x,y,out);
    }
    
    for (int k=0;k<1e6;k++)
    {
        //printf("k = %d",k);
        Cp2(x,y,out);
    }
    
    printf("Cp = %f\n",out[0]);
    printf("dCp/dbeta = %f \n",out[1]);
    printf("dCp/dlambda = %f \n",out[2]);
    printf("d2Cp/dbeta2 = %f \n",out[3]);
    printf("d2Cp/dlambda2 = %f \n",out[4]);
    printf("d2Cp/dbetadlambda = %f \n",out[5]);

    
    // Ct table
    for (int k=0;k<1e6;k++)
    {
        //printf("k = %d",k);
        Ct0(x,y,out);
    }
    
    for (int k=0;k<1e6;k++)
    {
        //printf("k = %d",k);
        Ct1(x,y,out);
    }
    
    for (int k=0;k<1e6;k++)
    {
        //printf("k = %d",k);
        Ct2(x,y,out);
    }
    
    printf("Ct = %f\n",out[0]);
    printf("dCt/dbeta = %f \n",out[1]);
    printf("dCt/dlambda = %f \n",out[2]);
    printf("d2Ct/dbeta2 = %f \n",out[3]);
    printf("d2Ct/dlambda2 = %f \n",out[4]);
    printf("d2Ct/dbetadlambda = %f \n",out[5]);


}
*/