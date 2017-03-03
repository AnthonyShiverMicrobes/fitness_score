function y = splineFitFunc(x,p,k)
%SPLINEFITFUNC is the function used for splinefit
%
%written by Sean Collins (2006) as part of the EMAP toolbox

y=spline(k,p,x);
