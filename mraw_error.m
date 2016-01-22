function raw = mraw_error(datamat,datameta,field,format,numDup)
%-------------------------------------------------------------------------!
% raw=mraw_error(datamat,datameta,field,format);                          :
%-------------------------------------------------------------------------!
% mraw_error.m                                                            :
% ------------                                                            :
% this script builds the minimum structure required to estimate error     :
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
raw.meansize=datameta.([field,'_men_',format]);
raw.size=datamat.(field);
raw.numDup=numDup;
%raw.row=datameta.row;
%raw.col=datameta.col;
end
