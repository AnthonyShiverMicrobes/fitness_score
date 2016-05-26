function sd=myNanSD(arr)
%returns the standard deviation of a list, ignoring NaN values
%
%written by Sean Collins (2006) as part of the EMAP toolbox

sd=sqrt(myNanVar(arr));
