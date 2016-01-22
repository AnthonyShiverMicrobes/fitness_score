function mat=remove_bad_strains(mat,bad)
fields=fieldnames(mat);
for m=1:length(fields)
    mat.(fields{m})(:,bad,:)=NaN;
end
