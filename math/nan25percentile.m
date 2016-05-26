function res=nan25percentile(list)
%returns the value at the 25th percentile in the list
list=list(~isnan(list));
list=sort(list);
len=length(list);
ind50=round((len+0.99)*0.5);
res=median(list(1:ind50));
