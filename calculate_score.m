function scoremat=calculate_score(datamatF,datamatV,datameta,variance_limit,numDup)
fields=fieldnames(datamatF); %assume V and F have same fields

for m=1:length(fields)
    raw=mraw_error(datamatF,datameta,fields{m},'F',numDup);
    rawN=mraw_error(datamatV,datameta,fields{m},'V',numDup);
    rawS=mraw_score(datamatV,datameta,fields{m},'V');
    %score data
    err=computeErrorEstimates(raw,rawN);
    scores=computeScores(rawS,err,'median');
    %Remove scores resulting from unreliably noisy single measurements
    scores.data(rawS.sdsize>variance_limit)=NaN;
    scoremat.(fields{m})=scores.data;
end
    