function res=nan75percentile(list)
%returns the value at the 75th percentile in the list
list=list(~isnan(list));
list=sort(list);
len=length(list);
ind50=round((len+1)*0.5);
res=median(list(ind50:len));
