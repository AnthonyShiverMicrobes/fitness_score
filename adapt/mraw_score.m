function raw=mraw_score(datamat,datameta,field,format)
%-------------------------------------------------------------------------!
% raw=mraw_score(datamat,datameta,field,format);                          :
%-------------------------------------------------------------------------!
% mraw_score.m                                                            :
% ------------                                                            :
% this script builds the minimum structure required to score the data     :
% using scripts from Sean Collins.                                        :
%-------------------------------------------------------------------------!
% format                                                                  :
% ------                                                                  :
% 'F': the unnormalized data                                              :
% 'V': the normalized data                                                :
% this variable controls which parameters are extracted from datameta.    :
%-------------------------------------------------------------------------!
% Anthony Shiver (2014)                                                   :
%-------------------------------------------------------------------------!
raw.rowlabels=datameta.cnd;
raw.collabels=datameta.mut;
raw.sdsize=datameta.([field,'_std_',format]);
raw.size=datamat.(field);
raw.rep_num='dummy';
end