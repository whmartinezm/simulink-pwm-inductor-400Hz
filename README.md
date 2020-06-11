# simulink-pwm-inductor-400Hz
The attached MATLAB / Simulink models demonstrate the implementation of the inductor iron-loss model discussed in manuscript “Simulink Model for PWM-Supplied Laminated Magnetic Cores Including Hysteresis, Eddy-Current and Excess Losses” and in "Iron Loss Evaluation of GaN PWM-Supplied Magnetic Cores for MHz Converters in More Electric Aircrafts - 400Hz case". They have been tested with MATLAB R2016b. The following compilers have been used to compile the hysteresis model:
-	Microsoft Visual C++ 2013 Professional
-	Intel Parallel Studio XE 2016 for Fortran with Microsoft Visual Studio 2013

The relevant files are briefly described below: 

- run_inductor.m
  
  Script file for initializing the model parameters and running the Simulink models.

- inductor.slx
  
  Model for PWM full-bridge converter supplied inductor, in which the skin effect in the core laminations is neglected.

- inductor_skineffect.slx
  
  Model for PWM full-bridge converter supplied inductor, in which the skin effect in the core laminations is taken into account.

- hyst/build_sfunc.m
  
  Script file for building the C MEX S-function gateway for the hysteresis model.

- hyst/HDHM.F
  
  FORTRAN implementation of the history-dependent hysteresis model described by Zirka et al. (2014).

- hyst/HDHM_sfunc_f.F
  
  FORTRAN implementation of the S-function interface to hysteresis model.

- hyst/HDHM_sfunc.c
  
  Level 2 S-function gateway to the FORTRAN hysteresis model.

- hyst/hydata_nippon_35h300.txt
  
  Text file containing magnetization curves for the studied 35H300 steel grade.
