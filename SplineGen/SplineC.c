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

float Pijeval(const float mat[],int i,int j,const int n)
{
// Line major "matrix" P -> take out i,j entry
    
    printf("index = %d\n", i*n + j);
    return mat[i*n + j];
}

float Blend(float x_basis[], float y_basis[], int i_x, int i_y, const float P[], const int n, int length_x_basis, int length_y_basis)
{
    printf("Blending \n");
    float S = 0;
    for (int k1 = 0; k1 < length_x_basis; k1++)
    {
        for (int k2 = 0; k2 < length_y_basis; k2++)
        {
            S += Pijeval(P,i_x-k1,i_y-k2,n)*x_basis[k1]*y_basis[k2];
        }
    }
    
    return S;
}



int basisFuncs(float basis[], float xi, const int order, const float U[], int i)
{
    //Minimal implementation of the Cox-deBoor formula
       
    //
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



int EvalSpline(float x, float y, float out[])
{
    
    #include "SplineData.h"
    
    int ix = findspan(x, knots_x, length_knots_x);
    int iy = findspan(y, knots_y, length_knots_y);
     
    
    float basis_x[p+1];
    float basis_y[q+1];
    
    basisFuncs(basis_x, x, p, knots_x, ix);
    basisFuncs(basis_y, y, q, knots_y, iy);
    
    
    //printf("Address basis_y[-1]: %p \n", &basis_y[q]);
    out[0] = Blend(basis_x, basis_y, ix, iy, P, n, p+1, q+1);

    int ix_tilde = findspan(x, Ux, length_Ux);
    int iy_tilde = findspan(y, Uy, length_Uy);
    
    float basis_x_tilde[p] = {0};
    float basis_y_tilde[q] = {0};
    
    basisFuncs(basis_x_tilde, x, p-1, Ux, ix_tilde);
    basisFuncs(basis_y_tilde, y, q-1, Uy, iy_tilde);
    
    
    out[1] = Blend(basis_x_tilde, basis_y,       ix_tilde, iy,       Px, n,   p,   q+1);
    out[2] = Blend(basis_x      , basis_y_tilde, ix,       iy_tilde ,Py, n-1, p+1, q  );
    
    return 0;
}

int main()
{
    float x = 5.;
    float y = 3.;
    
    float out[3];
    
    EvalSpline(x,y,out);
    
    printf("s = %f\n",out[0]);
    printf("dsdx = %f \n",out[1]);
    printf("dsdy = %f \n",out[2]);
}
