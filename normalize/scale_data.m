function datamatV = scale_data(datamat,method)
%-------------------------------------------------------------------------!
% datamatV=scale_data(datamat,method);                                    :
%-------------------------------------------------------------------------!
% scale_data.m                                                            :
% ------------                                                            :
% scale_data takes the input datamat, normalizes the variance so that     :
% all plate distributions are equal ('method')                            :
% and finally shifts all distrubutions to ensure there are no negative    :
% values.                                                                 :
%-------------------------------------------------------------------------!
% method                                                              :
% ----------                                                              :
% 'iqr': distance between 25th and 75th percentile                        :
% 'mad': (m)edian (a)bsolute (d)eviation, median deviation from the median:
%-------------------------------------------------------------------------!
% Anthony Shiver (2014)                                                   :
%-------------------------------------------------------------------------!
fields=fieldnames(datamat);

for m=1:length(fields)
    datamatV.(fields{m})=datamat.(fields{m}).*NaN;
    S=size(datamat.(fields{m}));
    var=ones(S(1),S(3))*NaN;
    % shift data to center (median) at zero and calculate variance
    for k=1:S(3)
        datamatV.(fields{m})(:,:,k)=datamat.(fields{m})(:,:,k) - ...
                    (nanmedian(datamat.(fields{m})(:,:,k),2)*ones(1,S(2)));
        for i=1:S(1)
            switch method
                case 'iqr'
                    list=datamatV.(fields{m})(i,:,k)';
                    var(i,k)=nan75percentile(list)-nan25percentile(list);
                case 'mad'
                    list=abs(datamatV.(fields{m})(i,:,k)-nanmedian(datamatV.(fields{m})(i,:,k),2));
                    var(i,k)=nanmedian(list,2);
            end
        end
    end
    %calculate multiplicative factor to scale variance
    mean_var=nanmean(var(:));
    multi = var ./ mean_var;
    %scale variance
    for k=1:S(3)
        DIV=multi(:,k)*ones(1,S(2));
        datamatV.(fields{m})(:,:,k)=datamatV.(fields{m})(:,:,k)./DIV;
    end
    % shift scaled colony sizes so that all measurements are positive again
    min_sze=min(datamatV.(fields{m})(:));
    datamatV.(fields{m})=datamatV.(fields{m})+abs(min_sze)+1;
end

end
