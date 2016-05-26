function datamat = filter_flags(datamat,datameta)
%-------------------------------------------------------------------------!
% datamat=filter_flags(datamat,datameta)                                  :
%-------------------------------------------------------------------------!
% filter_flags.m is an optional filter script that uses the FLAG field of :
% datamat (datamat.flg) to remove elements of the array that were manually:
% identified as bad before analysis. The flags come in the form of a      :
% format string, [r1-r2:c1-c2], which specify the edges of a rectangle    :
% chunk_index.m is used to convert this format string to a linear index   :
% that filter_flags.m then uses to remove elements from all appropriate   :
% data fields.                         :
%--------------------------------------!
% Anthony Shiver (2013)                :
%--------------------------------------!
fields=fieldnames(datamat);
flag_size=size(datameta.flg);
for i = 1 : flag_size(1)
    for k = 1 : flag_size(2)
        %ignore "-" and empty fields
        if(~strcmp('-',datameta.flg{i,k})&& ...
                ~strcmp('',datameta.flg{i,k})&&...
                ~strcmp('0',datameta.flg{i,k})&&...
                ~isempty(datameta.flg{i,k}))
            if(strcmp('discard',datameta.flg{i,k}))
                for m=1:length(fields)
                        datamat.(fields{m})(i,:,k)=NaN;
                end
            else
                ind=chunk_index([max(datameta.row),max(datameta.col)],datameta.flg{i,k});
                for m=1:length(fields)
                        datamat.(fields{m})(i,ind,k)=NaN;
                end
            end
        end
    end
end
end
           