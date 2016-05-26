function v=myNanVar(arr)
%computes the variance of a list, ignoring NaN values
%
%written by Sean Collins (2006) as part of the EMAP toolbox

arr=arr(~isnan(arr));
if length(arr)>1
    v=var(arr);
else
    v=NaN;
end
