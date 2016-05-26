function med=myNanMedian(arr)
%returns the median value of a list, ignoring NaN values
%
%written by Sean Collins (2006) as part of the EMAP toolbox

arr=arr(~isnan(arr));
if length(arr)>0
    med=median(arr);
else
    med=NaN;
end
