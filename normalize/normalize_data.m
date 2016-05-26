function [datamatN,datameta] = normalize_data(datamat,datameta,method,ignore,ignore_field)
%-------------------------------------------------------------------------!
% [datamatN,datameta]=normalize_data(datamat,datameta,method,...          :
%                                   ignore,ignore_field);                 :
%-------------------------------------------------------------------------!
% normalize_data.m takes the datamat passed to it, ignores the values     :
% passed through "ignore" and "ignore_field" (nulls the values whose      :
% ignore_field matches those in ignore), estimates the average size using :
% "method", and normalizes the plates so that they have the same average. :
%-------------------------------------------------------------------------!
% "method" values accepted                                                :
% ------------------------                                                :
% 'median': estimate median of distribution                               :
% 'middlemean': mean of 25th to 75th percentile                           :
%-------------------------------------------------------------------------!
% Anthony Shiver (2014)                                                   :
%-------------------------------------------------------------------------!

%build necessary datastructures
fields=fieldnames(datamat);
datamat_size=size(datamat.(fields{1})); %use first field as presumptive size
for m=1:length(fields)
    datamatN.(fields{m})=NaN.*datamat.(fields{m});
    datameta.([fields{m},'_nrm'])=NaN.*ones(datamat_size(1),1,datamat_size(3)); %datameta will store the average colony size parameter (datameta.field_nrm)
    norm.(fields{m})=datamat.(fields{m}); %norm will contain the clipped datasize to estimate average colony size
end

%remove fields that need to be ignored for normalization
if ~isempty(ignore) && isfield(datameta,ignore_field)
    for m=1:length(ignore)
        ind=strcmpi(ignore{m},datameta.(ignore_field));
        for n=1:length(fields)
            norm.(fields{n})(:,ind,:)=NaN;
        end
    end
end

%determine target for normalization
for m=1:length(fields)
    target.(fields{m})=nanmedian(datamat.(fields{m})(:));
end

%normalize entire plate
for m = 1 : length(fields)
    for i = 1 : datamat_size(1)
        for k = 1 : datameta.rep(i)
            switch method
                case 'median'
                    datameta.([fields{m},'_nrm'])(i,1,k)=nanmedian(norm.(fields{m})(i,:,k),2);
                case 'middlemean'
                    upqr=nan75percentile(norm.(fields{m})(i,:,k));
                    lwqr=nan25percentile(norm.(fields{m})(i,:,k));
                    platerange=(norm.(fields{m})(i,:,k)<upqr & norm.(fields{m})(i,:,k)>lwqr);
                    datameta.([fields{m},'_nrm'])(i,1,k)=nanmean(norm.(fields{m})(i,platerange,k),2);
                otherwise %median assumed
                    datameta.([fields{m},'_nrm'])(i,1,k)=nanmedian(norm.(fields{m})(i,:,k),2);
            end
            datamatN.(fields{m})(i,:,k)=datamat.(fields{m})(i,:,k).*(target.(fields{m})/datameta.([fields{m},'_nrm'])(i,1,k));
        end
    end
end

