function key=read_condition_key(filepath)
%-------------------------------------------------------------------------!
% condition_meta=read_condition_key(filepath)                             :
%-------------------------------------------------------------------------!
% read_condition_key reads in metadata associated with a chemical genetic :
% screen.                                                                 :
% Expected File Format:                                                   :
% HeaderLines = 1                                                         : 
% FileName | Condition | Concentration | Source | FLAG | Batch            :
%-------------------------------------------------------------------------!
% FileName:      Name of data file containing colony sizes                :
% Condition:     Name of condition, used for grouping experiments, very   :
%                important that this be consistent.                       :
% Concentration: Concentration of chemical, not used for comparison so its:
%                important that two identical chemical conditions in the  :
%                same batch be given unique names or they will be grouped :
% Source:        String describing the plate that was used to pin onto    :
%                the particular experiment.                               :
% FLAG:          String flagging elements of the array for removal.       :
% Batch:         Number describing the batch a particular experiment is   :
%                associated with. Used to match experiments from different:
%                arrays (K1&K2) without combining the same condition from :
%                unique batches.                                          :
%-------------------------------------------------------------------------!
% Anthony Shiver (2013)                                                   :
%-------------------------------------------------------------------------!
fid=fopen(filepath,'r');
cell=textscan(fid,'%s%s%s%s%s%s','Delimiter','\t','HeaderLines',1);
[key.nme,key.cnd,key.cnc,key.src,key.flg,key.bch]=cell{:};
end
