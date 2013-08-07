//
//  Blender.h
//  
//
//  Created by Sebastien Gros on 7/27/13.
//
//

float Blend44(float x_basis[], float y_basis[], int i_x, int i_y, const float P[], const int n)
{
    float S = 0;
    float x_basis_k1;
    
    // Access: P[i*n + j]
    
    //k1 = 0
    x_basis_k1 = x_basis[0];
    S += P[i_x*n + i_y  ]*x_basis_k1*y_basis[0];
    S += P[i_x*n + i_y-1]*x_basis_k1*y_basis[1];
    S += P[i_x*n + i_y-2]*x_basis_k1*y_basis[2];
    S += P[i_x*n + i_y-3]*x_basis_k1*y_basis[3];
    
    //k1 = 1
    x_basis_k1 = x_basis[1];
    S += P[(i_x-1)*n + i_y  ]*x_basis_k1*y_basis[0];
    S += P[(i_x-1)*n + i_y-1]*x_basis_k1*y_basis[1];
    S += P[(i_x-1)*n + i_y-2]*x_basis_k1*y_basis[2];
    S += P[(i_x-1)*n + i_y-3]*x_basis_k1*y_basis[3];
    
    //k1 = 2
    x_basis_k1 = x_basis[2];
    S += P[(i_x-2)*n + i_y  ]*x_basis_k1*y_basis[0];
    S += P[(i_x-2)*n + i_y-1]*x_basis_k1*y_basis[1];
    S += P[(i_x-2)*n + i_y-2]*x_basis_k1*y_basis[2];
    S += P[(i_x-2)*n + i_y-3]*x_basis_k1*y_basis[3];
    
    //k1 = 3
    x_basis_k1 = x_basis[3];
    S += P[(i_x-3)*n + i_y  ]*x_basis_k1*y_basis[0];
    S += P[(i_x-3)*n + i_y-1]*x_basis_k1*y_basis[1];
    S += P[(i_x-3)*n + i_y-2]*x_basis_k1*y_basis[2];
    S += P[(i_x-3)*n + i_y-3]*x_basis_k1*y_basis[3];
    
    return S;
}




float Blend43(float x_basis[], float y_basis[], int i_x, int i_y, const float P[], const int n)
{
    float S = 0;
    float x_basis_k1;
    
    // Access: P[i*n + j]
    
    //k1 = 0
    x_basis_k1 = x_basis[0];
    S += P[i_x*n + i_y  ]*x_basis_k1*y_basis[0];
    S += P[i_x*n + i_y-1]*x_basis_k1*y_basis[1];
    S += P[i_x*n + i_y-2]*x_basis_k1*y_basis[2];
    
    //k1 = 1
    x_basis_k1 = x_basis[1];
    S += P[(i_x-1)*n + i_y  ]*x_basis_k1*y_basis[0];
    S += P[(i_x-1)*n + i_y-1]*x_basis_k1*y_basis[1];
    S += P[(i_x-1)*n + i_y-2]*x_basis_k1*y_basis[2];
    
    //k1 = 2
    x_basis_k1 = x_basis[2];
    S += P[(i_x-2)*n + i_y  ]*x_basis_k1*y_basis[0];
    S += P[(i_x-2)*n + i_y-1]*x_basis_k1*y_basis[1];
    S += P[(i_x-2)*n + i_y-2]*x_basis_k1*y_basis[2];
    
    //k1 = 3
    x_basis_k1 = x_basis[3];
    S += P[(i_x-3)*n + i_y  ]*x_basis_k1*y_basis[0];
    S += P[(i_x-3)*n + i_y-1]*x_basis_k1*y_basis[1];
    S += P[(i_x-3)*n + i_y-2]*x_basis_k1*y_basis[2];
    
    return S;
}



float Blend34(float x_basis[], float y_basis[], int i_x, int i_y, const float P[], const int n)
{
    float S = 0;
    float x_basis_k1;
    
    // Access: P[i*n + j]
    
    //k1 = 0
    x_basis_k1 = x_basis[0];
    S += P[i_x*n + i_y  ]*x_basis_k1*y_basis[0];
    S += P[i_x*n + i_y-1]*x_basis_k1*y_basis[1];
    S += P[i_x*n + i_y-2]*x_basis_k1*y_basis[2];
    S += P[i_x*n + i_y-3]*x_basis_k1*y_basis[3];
    
    //k1 = 1
    x_basis_k1 = x_basis[1];
    S += P[(i_x-1)*n + i_y  ]*x_basis_k1*y_basis[0];
    S += P[(i_x-1)*n + i_y-1]*x_basis_k1*y_basis[1];
    S += P[(i_x-1)*n + i_y-2]*x_basis_k1*y_basis[2];
    S += P[(i_x-1)*n + i_y-3]*x_basis_k1*y_basis[3];
    
    //k1 = 2
    x_basis_k1 = x_basis[2];
    S += P[(i_x-2)*n + i_y  ]*x_basis_k1*y_basis[0];
    S += P[(i_x-2)*n + i_y-1]*x_basis_k1*y_basis[1];
    S += P[(i_x-2)*n + i_y-2]*x_basis_k1*y_basis[2];
    S += P[(i_x-2)*n + i_y-3]*x_basis_k1*y_basis[3];
    
    
    return S;
}

float Blend33(float x_basis[], float y_basis[], int i_x, int i_y, const float P[], const int n)
{
    float S = 0;
    float x_basis_k1;
    
    // Access: P[i*n + j]
    
    //k1 = 0
    x_basis_k1 = x_basis[0];
    S += P[i_x*n + i_y  ]*x_basis_k1*y_basis[0];
    S += P[i_x*n + i_y-1]*x_basis_k1*y_basis[1];
    S += P[i_x*n + i_y-2]*x_basis_k1*y_basis[2];
    
    //k1 = 1
    x_basis_k1 = x_basis[1];
    S += P[(i_x-1)*n + i_y  ]*x_basis_k1*y_basis[0];
    S += P[(i_x-1)*n + i_y-1]*x_basis_k1*y_basis[1];
    S += P[(i_x-1)*n + i_y-2]*x_basis_k1*y_basis[2];
    
    //k1 = 2
    x_basis_k1 = x_basis[2];
    S += P[(i_x-2)*n + i_y  ]*x_basis_k1*y_basis[0];
    S += P[(i_x-2)*n + i_y-1]*x_basis_k1*y_basis[1];
    S += P[(i_x-2)*n + i_y-2]*x_basis_k1*y_basis[2];
    
    
    return S;
}



float Blend24(float x_basis[], float y_basis[], int i_x, int i_y, const float P[], const int n)
{
    float S = 0;
    float x_basis_k1;
    
    // Access: P[i*n + j]
    
    //k1 = 0
    x_basis_k1 = x_basis[0];
    S += P[i_x*n + i_y  ]*x_basis_k1*y_basis[0];
    S += P[i_x*n + i_y-1]*x_basis_k1*y_basis[1];
    S += P[i_x*n + i_y-2]*x_basis_k1*y_basis[2];
    S += P[i_x*n + i_y-3]*x_basis_k1*y_basis[3];
    
    //k1 = 1
    x_basis_k1 = x_basis[1];
    S += P[(i_x-1)*n + i_y  ]*x_basis_k1*y_basis[0];
    S += P[(i_x-1)*n + i_y-1]*x_basis_k1*y_basis[1];
    S += P[(i_x-1)*n + i_y-2]*x_basis_k1*y_basis[2];
    S += P[(i_x-1)*n + i_y-3]*x_basis_k1*y_basis[3];
        
    
    return S;
}



float Blend42(float x_basis[], float y_basis[], int i_x, int i_y, const float P[], const int n)
{
    float S = 0;
    float x_basis_k1;
    
    // Access: P[i*n + j]
    
    //k1 = 0
    x_basis_k1 = x_basis[0];
    S += P[i_x*n + i_y  ]*x_basis_k1*y_basis[0];
    S += P[i_x*n + i_y-1]*x_basis_k1*y_basis[1];
    
    //k1 = 1
    x_basis_k1 = x_basis[1];
    S += P[(i_x-1)*n + i_y  ]*x_basis_k1*y_basis[0];
    S += P[(i_x-1)*n + i_y-1]*x_basis_k1*y_basis[1];
    
    //k1 = 2
    x_basis_k1 = x_basis[2];
    S += P[(i_x-2)*n + i_y  ]*x_basis_k1*y_basis[0];
    S += P[(i_x-2)*n + i_y-1]*x_basis_k1*y_basis[1];
    
    //k1 = 3
    x_basis_k1 = x_basis[3];
    S += P[(i_x-3)*n + i_y  ]*x_basis_k1*y_basis[0];
    S += P[(i_x-3)*n + i_y-1]*x_basis_k1*y_basis[1];
    
    return S;
}