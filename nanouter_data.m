function datamat=nanouter_data(datamat,datameta)
%-------------------------------------!
% datamat=nanouter_data(datamat)      :
%-------------------------------------!
% nanouter_data.m is an optional      :
% filtering step that replaces data   :
% in the outer two rows and columns   :
% with NaN.                           :
%-------------------------------------!
% Anthony Shiver (2013)               :
%-------------------------------------!
%extract necessary information
datafields=fieldnames(datamat);
max_row=max(datameta.row(:))-2;
max_col=max(datameta.col(:))-2;
outerrange=((datameta.row<3|datameta.row>max_row)|...
            (datameta.col<3|datameta.col>max_col));
%NaN the outer edges
for i=1:length(datafields)
        datamat.(datafields{i})(:,outerrange,:)=NaN;
end
end