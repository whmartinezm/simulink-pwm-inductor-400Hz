% Builds the FORTRAN and C codes
clear all; close all; clc

% Current directory needs to be here
d = fileparts(mfilename('fullpath'));
if ~strcmp(cd, d)
  cd(d)
end

% Build
mex -setup FORTRAN
mex -c HDHM.F HDHM_sfunc_f.F
mex -setup C
mex HDHM_sfunc.c HDHM_sfunc_f.obj HDHM.obj
