function datamat=nomorezeros(datamat)
%-------------------------------------------------------------------------!
% Usage: datamat = nomorezeros(datamat);                                  :
%-------------------------------------------------------------------------!
% Before normalizing strains, replace all zeros with ones. Eliminates     :
% infinite values from dividing by zero.
%-------------------------------------------------------------------------!
f=fieldnames(datamat);
for i=1:numel(f)
    datamat.(f{i})(datamat.(f{i})==0)=1;
end
