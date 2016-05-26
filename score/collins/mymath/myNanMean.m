function m=myNanMean(arr)
%returns the mean of a list, ignoring NaN values
%
%written by Sean Collins (2006) as part of the EMAP toolbox

arr=arr(~isnan(arr));
if length(arr)>0
    m=mean(arr);
else
    m=NaN;
end
