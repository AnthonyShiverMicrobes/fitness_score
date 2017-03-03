function [f,pp] = splinefit(x,y,k,wt)
%SPLINEFIT fits cubic splines to the data (x,y) using knots k. The knots k
%are x coordinates. Initial guesses for y at the knots are the values at
%the closest points x.
%
%written by Sean Collins (2006) as part of the EMAP toolbox

k=k(:);
x=x(:);
y=y(:);
ind=~isnan(x) & ~isnan(y);
x=x(ind);
y=y(ind);
pin=ones(size(k))*mean(k);

stol=.0001;
niter=20;
if nargin<4
    wt=ones(size(y));
end

[f,p] = leasqr2(x,y,pin,'splineFitFunc',k,stol,niter,wt);
pp=spline(k,p);
