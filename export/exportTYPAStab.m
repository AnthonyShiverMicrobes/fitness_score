function exportTYPAStab(mat,meta,filename,field)
mat_size=size(mat.(field));
fid=fopen(filename,'w');
%print gene names
GENEname=strcat(meta.acc,'-',upper(meta.mut));
GENE=sprintf('\t%s', GENEname{1:length(GENEname)});
LINE=[ 'Condition' GENE ];
fprintf(fid,'%s\n',LINE);
%print row label then data
for i = 1 : mat_size(1)
    VECTOR=regexprep(mat2str(mat.(field)(i,:)),{'\s+','NaN','[',']'},{'\t','','',''});
    COND=[ char(meta.cnd(i)) ' [' char(meta.cnc(i)) '] {' char(meta.bch(i)) '}' ];
    LINE=[ COND sprintf('\t%s',VECTOR) ];
    fprintf(fid,'%s\n',LINE);
end
fclose(fid);
end

    