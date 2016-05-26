function datamatP=transform_data(datamat,power)
%-------------------------------------------------------------------------!
% datamatP=transform_data(datamat,power);                                 :
%-------------------------------------------------------------------------!
% Transform Data transforms all fields of a given datamat by power.       :
% This was suggested in the workflow from Nichols and Sen, with a         :
% suggested power of 0.5. Supply a power of 1 to change nothing.          :
%-------------------------------------------------------------------------!
% Anthony Shiver (2013)                                                   :
%-------------------------------------------------------------------------!
fields=fieldnames(datamat);
for i = 1 : length(fields);
    datamatP.(fields{i})=(datamat.(fields{i})).^power;
end