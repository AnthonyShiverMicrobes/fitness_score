function plate=data2plate2(dat,row,col)
%converts data into 32x48 plate format
%made for interaction with sean collins computeErrorEstimates
plate=ones(32,48)*NaN;
for i=1:length(dat)
    plate(row(i),col(i))=dat(i);
end
