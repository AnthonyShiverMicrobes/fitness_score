function key=read_array_key(filepath)
%-------------------------------------------------------------------------!
% array_meta=read_array_key(filepath)                                     :
%-------------------------------------------------------------------------!
% read_array_key reads in metadata associated with a plate format in the  :
% chemical genomics and genetic interaction screens (like RNAP1 or KEIO2 ):
% Expected File Format:                                                   :
% HeaderLines = 1                                                         :
% { Row | Column | 96-well plate | 96-well Well | Accession ID | ...      :
%   Clone ID | Mutation | MarkerName }                                    :
%-------------------------------------------------------------------------! 
% Row:           row number in array                                      :
% Column:        column number in array                                   :  
% 96-well plate: P1-P16, plate number of 96 well plate used to generate   :
%                the 1536 array. Only relevant for RNAP formats (so far). :
% 96-well well:  A01-G12, well ID in 96 well plate used to generate the   :
%                1536 array. Only relevant for RNAP formats (so far).     :
% Accession:     Accession name of mutated gene, e.g. ECK2983 for rpoB.   :
% Clone ID:      Unique identifier used to distinguish biological         :
%                replicates of particular strain. in scoreRNAP.m, this    :
%                field is used to group replicates of the same strain.    :
% Mutation:      Mnemonic identifying mutation, unique biological         :
%                replicates of the same mutation should have unique clone :
%                IDs but the same mutation.                               :
% MarkerName:    String identifying the marker used, e.g. ALS2*3 for rpoBC:
%-------------------------------------------------------------------------!
% Anthony Shiver (2013)                                                   :
%-------------------------------------------------------------------------!
fid=fopen(filepath,'r');
cell=textscan(fid,'%u%u%s%s%s%s%s%s','Delimiter','\t','HeaderLines',1);
[key.row,key.col,key.plt,key.wll,key.acc,key.uid,key.mut,key.mrk]=cell{:};
end
