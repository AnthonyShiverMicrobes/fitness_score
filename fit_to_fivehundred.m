function [datamat]=fit_to_fivehundred(datamat,datameta,method,numDup)
%-------------------------------------------------------------------------!
% datamat=fit_to_fivehundred(datamat,datameta,method,numDup);             :
%-------------------------------------------------------------------------!
% fit_to_fivehundred.m: This is (hopefully) a temporary script that will  :
% assist in interacting with Sean Collins' scoring algorithms. I am       :
% worried that the error estimates require a mean of 500 to accurately    :
% guage the error of the distribution.                                    :
%-------------------------------------------------------------------------!
% "method"                                                                :
%  ------                                                                :
% 'median': median of the distribution                                    :
% 'middlemean': mean of the 25th-75th percentile of the distribution      : 
%-------------------------------------------------------------------------!
% Anthony Shiver (2014)                                                   :
%-------------------------------------------------------------------------!

fields=fieldnames(datamat);
for m=1:length(fields)
    matsize=size(datamat.(fields{m}));
    for i=1:matsize(1)
        if numDup==1
            for k=1:datameta.rep(i)
                av=[];
                switch method
                    case 'median'
                        av=nanmedian(datamat.(fields{m})(i,:,k));
                    case 'middlemean'
                        p25=nan25percentile(datamat.(fields{m})(i,:,k));
                        p75=nan75percentile(datamat.(fields{m})(i,:,k));
                        range=(datamat.(fields{m})(i,:,k)>p25 & datamat.(fields{m})(i,:,k)<p75);
                        av=nanmean(datamat.(fields{m})(i,range,k));
                end
                datamat.(fields{m})(i,:,k)=(500/av).*datamat.(fields{m})(i,:,k);
            end
        elseif numDup==2 %added on for chemical genomics screen
            for k=1:2:matsize(3)
                av=[];
                data=datamat.(fields{m})(i,:,k:k+1);
                switch method
                    case 'median'
                        av=nanmedian(data(:));
                    case 'middlemean'
                        p25=nan25percentile(data(:));
                        p75=nan75percentile(data(:));
                        middle=data(data>p25 & data<p75);
                        av=nanmean(middle(:));
                end
                datamat.(fields{m})(i,:,k:k+1)=(500/av)*data;
            end
        end
    end
end

