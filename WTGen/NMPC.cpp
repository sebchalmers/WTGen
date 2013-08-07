/**
 *    
 *    OCTOBER 2012 SEBASTIEN GROS & MARIO ZANON
 *    
 */


#include <acado_toolkit.hpp>
#include <ocp_export.hpp>



int main( int argc, char * const argv[] ){


	//======================================================================
	// PARSE THE INPUT ARGUMENTS:
	// ----------------------------------
	double IC[3];
	
	/* arguments are passed to the main function by string
	 *  there are 'int argc' arguments
	 *		the first one is the name of the program
	 *		the next ones are arguments passed in the call like
	 *			program number1 number2
	 *  in this stupid simple parser, we directly read doubles from the arguments
	 */
	 
	int k=1;
	while (k < argc) {
		// load the i-th string argument passed to the main function
		char* input = argv[k];
		// parse the string to a double value
		IC[k-1] = atof(input);
		k++;
	}
	
	cout << "Initial Conditions" << endl;
	cout << "------------" << endl;
	for (k = 0; k < argc-1; ++k) {
		cout << k+1 << ":\t" << IC[k] << endl;
	}
	//======================================================================	
	
	const int Ncvp  = (int)IC[0];       // = 50
	const double Tc = IC[1];       // = 10
    const int IntegratorSteps = (int)IC[2];// =  2 * 50
	
    const double n_XD = 8;
    const double n_U = 5;
     
    USING_NAMESPACE_ACADO

    DifferentialState Omega;    // Rotor speed
    DifferentialState xT;       // Nacelle fore-aft position
    DifferentialState dxT;      // Nacelle fore-aft speed

    DifferentialState theta;    // Pitch angle
    DifferentialState dtheta;   // Pitch rate
    DifferentialState Mg;       // Generator Torque
    
    DifferentialState Sigma_Omega;    // Dummy state to match the slack
    DifferentialState Sigma_Pow;    // Dummy state to match the slack
 
    // CONTROL :
    // -------------------------
	Control	    ddtheta;        // Pitch acceleration
    Control     dMg;            // Generator torqe rate  
    Control     S_Omega;              // Slack variable
    Control     S_Pow;              // Slack variable
    
    Control     v0;             // LIDAR wind profile
  
     
    double Mgfactor = 1e4; //Scale generator torque
    
    //double v0 = 10;
    double PI = 3.141592653589793;
 // PARAMETERS :
    double rho = 1.23;
	//double omega = 2*PI;
    double R = 63;
    double eta = 0.95;
    double hH = 90;

    double mT = 347460;
    double mN = 240000;
    double mB = 17740;
    double mH = 56780;
    double xi = 0.7; 
    double iratio = 0.010309278350515;
    double JH = 115926;
    double JB = 11776047;
    double JG = 534.116;
    double f0 = 0.32;
    double ds = 0.01;
    
    double mTe = mT/4 + mN + mH + 3*mB; 
    double cT = 4*PI*mTe*ds*f0;
    double J = JH + 3*JB + JG/iratio/iratio;
    double kT = mTe*(2*PI*f0)*(2*PI*f0);
    
    IntermediateState   vrel   = v0 - dxT;   // Nacelle relative velocity wrt the wind field
    IntermediateState lambda = Omega*R/vrel; // Tip-speed ratio

    
    
// Interpolation of the look-up tables    
//Cp interpolation
    
    float out[6];
    
    #include CpCt
    
    Cp1(theta,lambda,out);
    IntermediateState Cp         = out[0];
    IntermediateState dCpdtheta  = out[1];
    IntermediateState dCpdlambda = out[2];
    

    Ct1(theta,lambda,out);
    IntermediateState Ct         = out[0];
    IntermediateState dCtdtheta  = out[1];
    IntermediateState dCtdlambda = out[2];
     
    IntermediateState MA = 0.5*rho*PI*R*R*R*Cp*vrel*vrel/lambda;
    IntermediateState FA = 0.5*rho*PI*R*R*Ct*vrel*vrel;
   
//// THE EQUATIONS OF MOTION:
//// ---------------------------------------------------------------

   
    //IntermediateState Pel = 1e-6*eta*Mg*Omega/i;  
  
    IntermediateState dOmega = (MA - Mgfactor*Mg/iratio)/J;   // Rotor acceleration
    IntermediateState ddxT = (FA - cT*dxT - kT*xT)/mTe;  // Tower acceleration
	
    double a_Omega = 1e-1;
    double a_Pow = 5e-2;
 
    
  // THE "RIGHT-HAND-SIDE" OF THE ODE:
// ---------------------------------------------------------------
   DifferentialEquation f;

    f  << dot(Omega)   ==  dOmega                             ;
    f  << dot(xT)	   ==  dxT                             ;
    f  << dot(dxT)     ==  ddxT                             ;
    
    f  << dot(theta)   ==  dtheta                            ;
    f  << dot(dtheta)  ==  ddtheta                            ;
    f  << dot(Mg)	   ==  dMg                            ;	   
 	
  	f  << dot(Sigma_Omega) ==  S_Omega - a_Omega*Sigma_Omega                            ;
    f  << dot(Sigma_Pow)   ==  S_Pow - a_Pow*Sigma_Pow                            ;
     
    // DEFINE AN OPTIMAL CONTROL PROBLEM:
    // ----------------------------------
    OCP ocp( 0.0, Tc, Ncvp );
    
    ocp.subjectTo( f );

	
    // CONSTRAINTS:
    // ---------------------------------

    ocp.subjectTo(  0  <=  Mg <=  1.05*4.282138171017305   );
	ocp.subjectTo(  theta_min <= theta <=  theta_max   );

    ocp.subjectTo(  -7.0 <= dtheta <= 7.0    );
        
    ocp.subjectTo(   0   <= (Omega/0.837758040957278 - 1) + S_Omega );
    ocp.subjectTo(   (Omega/1.267109036947883 -1) - S_Omega <= 0  );
     
	//ocp.subjectTo(  0  <=  S_Omega  );
    //ocp.subjectTo(  0  <=  S_Pow  );
    
    IntermediateState  Pelec = (1e-6*eta*Mgfactor/iratio)*Mg*Omega;
    IntermediateState dPelec = (1e-6*eta*Mgfactor/iratio)*(dMg*Omega+Mg*dOmega);
  
    double Prated = 5;     // Rated power
    double Pmax = 1.01*Prated;  // Maximum steady-state power
    
    ocp.subjectTo(  (Pelec/Prated - 1) - S_Pow <= 0 );
    
    ocp.subjectTo( Sigma_Pow - (Pmax - Prated)/a_Pow <= 0 );
    
    Function h;
    Function hN;

    ExportVariable SS("SS", n_XD + n_U, n_XD + n_U, REAL);
    ExportVariable SSN("SSN", n_XD, n_XD, REAL);

    const double K = 1.823540394474115e+06;
   
    hN << 0.5*M*zeta*sqrt(v0*v0*v0); // Power optimization
    hN << dxT;                       // Nacelle fore-aft speed
    hN << 1e-3*v0*v0*v0*Sigma_Pow;   // Accumulated constraint violation
    
    h << 0.5*M*zeta*sqrt(v0*v0*v0);  // Power optimization
    h << dxT;                        // Nacelle fore-aft speed
    h << 1e-3*v0*v0*v0*Sigma_Pow;       // Accumulated constraint violation  
    
    // Inputs
    h << ddtheta;                    // Pitch acceleration
    h << dMg;                        // Generator torqe rate 
        
    h << dPelec;          // dPower/dt
    h << 1e-3*v0*v0*v0*S_Omega;         // Slack variable
    h << 1e-3*v0*v0*v0*S_Pow;           // Slack variable
    //h << v0;              // LIDAR wind profile

    ocp.minimizeLSQ(SS, h);
    ocp.minimizeLSQEndTerm(SSN, hN);
   
    // Other attempts, a bit artificial
    //h << dMg + a*Mg - Omega*(2*dOmega+a*Omega)*(iratio*K/Mgfactor);   // Generator torqe rate
    //h << dMg  - 2*Omega*dOmega*(iratio*K/Mgfactor);   // Generator torqe rate
    
    
	// EXPORT THE CODE:
	// ---------------------------------
	printf("EXPORTING LINUX/QPOASES CODE\n");
	OCPexport mpc( ocp );

    mpc.set( CG_HARDCODE_CONSTRAINT_VALUES, NO );
    mpc.set( MAX_NUM_QP_ITERATIONS,       1000              );
    
	mpc.set( HESSIAN_APPROXIMATION,       GAUSS_NEWTON    );
	mpc.set( DISCRETIZATION_TYPE,   MULTIPLE_SHOOTING );
	mpc.set( QP_SOLVER,             QP_QPOASES      );
	mpc.set( HOTSTART_QP,           YES              );
    
//	mpc.set( INTEGRATOR_TYPE,             INT_RK4         );
 	mpc.set( INTEGRATOR_TYPE,             INT_IRK_GL2     ); // Rien's shite
	mpc.set( NUM_INTEGRATOR_STEPS,        IntegratorSteps );
	
    //Options for implicit integrators
	mpc.set( IMPLICIT_INTEGRATOR_NUM_ITS, 3				);
	mpc.set( IMPLICIT_INTEGRATOR_NUM_ITS_INIT, 0		);
	//mpc.set( LINEAR_ALGEBRA_SOLVER,		  HOUSEHOLDER_QR );
	mpc.set( UNROLL_LINEAR_SOLVER,        NO	      );
	mpc.set( IMPLICIT_INTEGRATOR_MODE, IFTR );
	
	mpc.set(SPARSE_QP_SOLUTION, FULL_CONDENSING);
    
    mpc.set( GENERATE_MATLAB_INTERFACE, YES );
    
	mpc.set( CG_USE_C99,    YES              );  
    
    
    
	mpc.exportCode( "code_export_MPC" );
	mpc.printDimensionsQP();

    return 0;
}



