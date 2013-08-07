//
//  Basisfunc.h
//  
//
//  Created by Sebastien Gros on 7/27/13.
//
//

int basisFuncs3(float basis[], float xi, const float U[], int i)
{
    //Minimal implementation of the Cox-deBoor formula
    
    int order = 3;
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
    basis[2] = 0;
    basis[3] = 0;

    
        // p = 2;
    
        
        //Update of Ni-p,p is standalone (cross arrow):
        Den  = Uiplus - U[i-1];
        basis[2] = (Uiplus - xi)*basis[1] / Den;
        

            int k = 1;
            //Flat arrow
            basis[1]  = (xi - U[i-1])*basis[1] / Den;
            
            //Cross arrow
            ipluspminusk = iplus+2-1;
            Den   =  U[ipluspminusk] - U[iplus-1];
            basis[1] += (U[ipluspminusk] - xi)*basis[0] / Den;
            
        
        // Update of Ni,p is standalone (flat arrow):
        basis[0] = (xi - Ui)*basis[0] / Den;
        
        
        //p = 3;
    
        //Update of Ni-p,p is standalone (cross arrow):
        Den  = Uiplus - U[i-2];
        basis[3] = (Uiplus - xi)*basis[2] / Den;
        
            // k = 2;
            //Flat arrow
            basis[2]  = (xi - U[i-2])*basis[2] / Den;
            
            //Cross arrow
            ipluspminusk = iplus+1;
            Den   =  U[ipluspminusk] - U[iplus-2];
            basis[2] += (U[ipluspminusk] - xi)*basis[1] / Den;

            // k = 1;
            //Flat arrow
            basis[1]  = (xi - U[i-1])*basis[1] / Den;
            
            //Cross arrow
            ipluspminusk = iplus+2;
            Den   =  U[ipluspminusk] - U[iplus-1];
            basis[1] += (U[ipluspminusk] - xi)*basis[0] / Den;

            
        
        // Update of Ni,p is standalone (flat arrow):
        basis[0] = (xi - Ui)*basis[0] / Den;

        
    
    return 0;
}