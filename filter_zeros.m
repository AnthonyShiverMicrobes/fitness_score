function [datamat,datameta] = filter_zeros(datamat,datameta,frac)
fields=fieldnames(datamat);
for m=1:length(fields)
    S=size(datamat.(fields{m}));
    experimental_means=nanmean(datamat.(fields{m}),3);
    zero_threshold=false(S(2),1);
    for j=1:S(2)
        zero_threshold(j)=( sum(experimental_means(:,j)==0)/sum(~isnan(experimental_means(:,j))) > frac);
    end
    datamat.(fields{m})(:,zero_threshold,:)=NaN;
    datameta.([fields{m},'_zerothreshind'])=zero_threshold;
end
