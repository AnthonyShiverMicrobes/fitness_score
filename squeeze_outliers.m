function datamat=squeeze_outliers(datamat,datameta)
%-------------------------------------------------------------------------!
% function datamat=squeeze_outliers(datamat,datameta)                     :
% ---------------------------------------------------                     :
% squeeze_outliers.mat searches for outliers that are > 3 times the size  :
% of their strain median (taken from cmd)                and sets their   :
% value to 3. This strategy was taken from the Collins workflow, but sets :
% the average relative to strain size, not plate average                  :
%-------------------------------------------------------------------------!
% Anthony Shiver 08-08-15                                                 :
%-------------------------------------------------------------------------!
fields=fieldnames(datamat);
for m=1:length(fields)
    cmd=nanmedian(nanmean(datamat.(fields{m}),3),1);
    relative.(fields{m})=datamat.(fields{m})*NaN;
    for i=1:length(datameta.rep)
        for k=1:datameta.rep(i)
            relative.(fields{m})(i,:,k)=datamat.(fields{m})(i,:,k)./cmd;
        end
    end
    outlierindex=(relative.(fields{m})>3);
    datamat.(fields{m})(outlierindex)=3; %set all outliers to a ceiling of 3
end