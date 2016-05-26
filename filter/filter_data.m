function [datamatF,datameta] = filter_data(datamat,datameta,plate,frac)
%-------------------------------------------------------------------------!
% datamatF=filter_data(datamat,datameta,plate);                           :
%-------------------------------------------------------------------------!
% filter_data.m is a driver script to perform the filtering steps         :
% (removal of specific array elements).                                   :
% datamat: data-matrix                                                    :
% datameta: meta-data                                                     :
% plate: type of plate being passed                                       :
%-------------------------------------------------------------------------!
% possible options for plate                                              :
% 'rnap': RNAP plate, remove outer two rows                               :
% 'keio': KEIO plate, keep outer two rows                                 :
% 'keio6':KEIO plate, distinction important later, keep outer two rows    :
%-------------------------------------------------------------------------!
% Anthony Shiver (2013)                                                   :
% (2015) added step to filter positions that are majority zero            :
%-------------------------------------------------------------------------!

%remove outer rows if RNAP plate
switch plate
    case 'rnap'
        datamatF=nanouter_data(datamat,datameta);
    case 'keio'
        datamatF=datamat;
    case 'keio6'
        datamatF=datamat;
    otherwise
        datamatF=datamat;
end
%remove positions that are flagged
datamatF=filter_flags(datamatF,datameta);
%remove positions that have majority zero
[datamatF,datameta]=filter_zeros(datamatF,datameta,frac);

