function score=normalizeIQR(score)
fields=fieldnames(score);
for m=1:length(fields)
    S=size(score.(fields{m}));
    for i=1:S(1);
        iqr=nan75percentile(score.(fields{m})(i,:)')-nan25percentile(score.(fields{m})(i,:)');
        score.(fields{m})(i,:)=score.(fields{m})(i,:).*(1.35/iqr);
    end
end