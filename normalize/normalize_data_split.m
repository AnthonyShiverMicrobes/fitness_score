function [datamatN,datameta] = normalize_data_split(datamat,datameta,method,ignore,ignore_field,ind)
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
indsize=size(ind);
%build necessary datastructures
fields=fieldnames(datamat);
datamat_size=size(datamat.(fields{1})); %use first field as presumptive size
for m=1:length(fields)
    datamatN.(fields{m})=datamat.(fields{m});
    datameta.([fields{m},'_nrm'])=NaN.*ones(datamat_size(1),indsize(2),datamat_size(3)); %datameta will store the average colony size parameter (datameta.field_nrm)
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
%determine target
for m=1:length(fields)
    target.(fields{m})=nanmedian(norm.(fields{m})(:));
end
%iterate through the two indices and normalize separately.
for a=1:indsize(2)
    %estimate normalization factor for subpopulation
    for m = 1 : length(fields)
        for i = 1 : datamat_size(1)
            for k = 1 : datameta.rep(i)
                switch method
                    case 'median'
                        datameta.([fields{m},'_nrm'])(i,a,k)=nanmedian(norm.(fields{m})(i,ind(:,a),k),2);
                    case 'middlemean'
                        upqr=nan75percentile(norm.(fields{m})(i,ind(:,a),k));
                        lwqr=nan25percentile(norm.(fields{m})(i,ind(:,a),k));
                        platerange=(norm.(fields{m})(i,:,k)<upqr & norm.(fields{m})(i,:,k)>lwqr) & ind(:,a)';
                        datameta.([fields{m},'_nrm'])(i,a,k)=nanmean(norm.(fields{m})(i,platerange,k),2);
                    otherwise %median assumed
                        datameta.([fields{m},'_nrm'])(i,a,k)=nanmedian(norm.(fields{m})(i,ind(:,a),k),2);
                end
                %prevent the formation of inf values by dividing through by
                %zero
                if datameta.([fields{m},'_nrm'])(i,a,k)==0
                    datameta.([fields{m},'_nrm'])(i,a,k)=NaN;
                end
            %build normalization factor for subpopulation
            norm_vector=ones(size(datamat.(fields{m})(i,:,k)));
            norm_vector(ind(:,a))=target.(fields{m})/datameta.([fields{m},'_nrm'])(i,a,k);
            %normalize subpopulation, leaving other intact
            datamatN.(fields{m})(i,:,k)=datamatN.(fields{m})(i,:,k).*norm_vector;
            end
        end
    end
end

