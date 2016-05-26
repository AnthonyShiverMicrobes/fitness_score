function datamat=enforce_triplicates(datamat)
fields=fieldnames(datamat);
for m=1:length(fields)
    matsize=size(datamat.(fields{m}));
    notnan=~isnan(datamat.(fields{m}));
    numofreplicates=sum(notnan,3);
    removalindex=numofreplicates<3;
    for k=1:matsize(3)
        replacement=datamat.(fields{m})(:,:,k);
        replacement(removalindex)=NaN;
        datamat.(fields{m})(:,:,k)=replacement;
    end
end
end
        
        
    