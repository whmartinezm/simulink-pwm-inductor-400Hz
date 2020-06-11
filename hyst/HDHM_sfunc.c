// Level 2 CMEX S-function gateway to the Fortran hysteresis model

#define S_FUNCTION_NAME  HDHM_sfunc
#define S_FUNCTION_LEVEL 2
#include "simstruc.h"

#ifdef _WIN32
#define HDHM_sfunc_init_ HDHM_SFUNC_INIT
#define HDHM_sfunc_ HDHM_SFUNC
#define HDHM_sfunc_terminate_ HDHM_SFUNC_TERMINATE
#endif

extern void HDHM_sfunc_terminate_();
extern void HDHM_sfunc_init_(int *nmodel, int *nrevmax);
extern void HDHM_sfunc_(float *bin, float *hout, float *dhout, int *iup);

typedef enum {IHY=0, NUM_SPARAMS } paramIndices;
#define IHY(S) (ssGetSFcnParam(S, IHY))

static void mdlInitializeSizes(SimStruct *S)
{
  ssSetNumSFcnParams(S,NUM_SPARAMS);
#if defined(MATLAB_MEX_FILE)
    if (ssGetNumSFcnParams(S) != ssGetSFcnParamsCount(S)) goto EXIT_POINT;
#endif

    {
      int iParam = 0;
      int nParam = ssGetNumSFcnParams(S);

      for ( iParam = 0; iParam < nParam; iParam++ )
      {
          ssSetSFcnParamTunable( S, iParam, SS_PRM_SIM_ONLY_TUNABLE );
      }
    }

  ssSetNumContStates( S, 0 );
  ssSetNumDiscStates( S, 0 );

  ssSetNumInputPorts(S, 1);
  ssSetInputPortWidth(S, 0, DYNAMICALLY_SIZED);
  ssSetInputPortDirectFeedThrough(S, 0, 1);
  ssSetInputPortRequiredContiguous(S, 0, 1);

  ssSetNumOutputPorts(S, 1);
  ssSetOutputPortWidth(S, 0, DYNAMICALLY_SIZED);

  ssSetSimStateCompliance(S, USE_DEFAULT_SIM_STATE);

  EXIT_POINT:
    return;
}


static void mdlInitializeSampleTimes(SimStruct *S)
{
    ssSetSampleTime(S, 0, INHERITED_SAMPLE_TIME);
    ssSetOffsetTime(S, 0, 0.0);
    ssSetModelReferenceSampleTimeDefaultInheritance(S);
}


// Function for running the HDHM
static void mdlOutputs(SimStruct *S, int_T tid)
{

  // Inputs and outputs 
  const double *bvec = (const double *) ssGetInputPortSignal(S,0);
  double *hvec  = (double *) ssGetOutputPortRealSignal(S,0);
  int     w  = ssGetInputPortWidth(S,0);
  int     k, yks, ind, iup;
  float   b, hout, dhout, t;
  
  // One
  yks = 1;

  // Which model to use
  //   1 -> Hysteresis model
  //   0 -> single-valued model
  int ihy = (int) mxGetScalar(IHY(S));

  // Time  
  t = ssGetTaskTime(S,0); 

  // Run through input vector
  for (k=0; k<w; k++) {

    // Input b and model index
    b   = (float) bvec[k];
    ind = k+1;

    // Run HDHM
    HDHM_sfunc_(&t, &yks, &ind, &b, &hout, &dhout, &ihy, &iup);
  
    // Assemble output vector
    hvec[k]   = (double) hout;
  }
}

static void mdlTerminate(SimStruct *S)
{ 
    HDHM_sfunc_terminate_();
}

#define MDL_INITIALIZE_CONDITIONS 
#if defined(MDL_INITIALIZE_CONDITIONS) && defined(MATLAB_MEX_FILE) 
static void mdlInitializeConditions(SimStruct *S)
{
 int nmodel = ssGetInputPortWidth(S,0);
 int ihy    = (int) mxGetScalar(IHY(S));
 int nrevmax = 30;

 HDHM_sfunc_init_(&nmodel, &nrevmax);
}
#endif


#ifdef  MATLAB_MEX_FILE 
#include "simulink.c" 
#else
#include "cg_sfun.h"
#endif


