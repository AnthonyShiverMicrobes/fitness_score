function datameta = controlstats(datamat,datameta,matformat)
%-------------------------------------------------------------------------!
% datameta=controlstats(datamat,datameta,matformat)                       :
%-------------------------------------------------------------------------!
% controlstats.m                                                          :
% --------------                                                          :
% controlstats adds metadata about parameter distributions to datameta.   :
% the 'matformat' input allows multiple statistical measurements to be    :
% stored in datameta simultaneously.                                      :
%-------------------------------------------------------------------------!
% 'matformat'                                                             :
%  ---------                                                              :
%  'F': datamatF                                                          :
%  'V': datamatV                                                          :
%-------------------------------------------------------------------------!
% Anthony Shiver (2014)                                                   :
%-------------------------------------------------------------------------!
fields=fieldnames(datamat);
if ~strcmpi(matformat,'F')&&~strcmpi(matformat,'V')
    return
end

for m=1:length(fields)
    datameta.([fields{m},'_men_',matformat])=nanmean(datamat.(fields{m}),3);
    datameta.([fields{m},'_std_',matformat])=nanstd(datamat.(fields{m}),0,3);
    datameta.([fields{m},'_cmd_',matformat])=nanmedian(datameta.([fields{m},'_men_',matformat]));
    datameta.([fields{m},'_csd_',matformat])=nanmedian(datameta.([fields{m},'_std_',matformat]));
end