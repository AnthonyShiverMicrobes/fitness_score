function datamatN = normalize_strain(datamat,datameta,strain_method,norm_method,ignore,ignore_field)
%-------------------------------------------------------------------------!
% datamatN = normalize_strain(datamat,datameta,strain_method,...          :
%                             norm_method,ignore,ignore_field);           :
%-------------------------------------------------------------------------!
% normalize_strain.m                                                      :
% ------------------                                                      :
% This script was designed to normalize according to special needs of the :
% plate. 'datamat' is the data to be analyzed, 'datameta' the metadata,   :
% 'strain_method' is the switch for which data needs to be normalized     :
% 'norm_method' is the alg. for finding the center of the distribution    :
% 'ignore' is the cell array of mutations to disregard during norm.       :
% and 'ignore_field' is the field to search using 'ignore'.               :
%-------------------------------------------------------------------------!
% strain_method                                                           :
% -------------                                                           :
% RNAP_marker: use the 'control' strains from RNAP to normalize all       :
%              strains carrying the same marker i.e. rpoBC{cat} for 2*3   :
% mixed: use the 'control' strains for non-core mutations and the         :
%        entire set of core mutants for beta and beta prime.              :
% mixed_structure: use the 'control' strains for non-core, and split the  :
%                  core mutations among three groups: those that are      :
%                  similar to rpoBC{cat}, the group sensitized to         :
%                  pseudomonic acid (pA), and the group resistant to pA.  :
%-------------------------------------------------------------------------!
% norm_method                                                             :
% -----------                                                             :
% 'median': median of distribution                                        :
% 'middlemean': mean of 25th to 75th percentile.                          :
%-------------------------------------------------------------------------!
% Anthony Shiver (2014)                                                   :
%-------------------------------------------------------------------------!
% notes:                                                                  
%
%-------------------------------------------------------------------------!
%------> Main Body
fields=fieldnames(datamat);
%remove ignore before further manipulation
for m=1:length(fields)
    norm.(fields{m})=datamat.(fields{m});
    datamatN.(fields{m})=NaN.*datamat.(fields{m});
    for n=1:length(ignore)
        ind=strcmpi(ignore{n},datameta.(ignore_field));
        norm.(fields{m})(:,ind,:)=NaN;
    end
end
%normalize strains
switch strain_method
    case 'RNAP_marker'
        [controls,markers,control_field,marker_field]=populate_search('RNAP_marker');
        datamatN=normalize_driver(norm,datamat,datameta,norm_method,controls,markers,...
                                  control_field,marker_field);
    case 'mixed'
        [controls,markers,control_field,marker_field]=populate_search('marker');
        datamatN=normalize_driver(norm,datamat,datameta,norm_method,controls,markers,...
                                  control_field,marker_field);
        [controls,markers,control_field,marker_field]=populate_search('RNAPgroup');
        datamatN=normalize_driver(norm,datamatN,datameta,norm_method,controls,markers,...
                                  control_field,marker_field);
    case 'mixed_structure'
        [controls,markers,control_field,marker_field]=populate_search('marker');
        datamatN=normalize_driver(norm,datamat,datameta,norm_method,controls,markers,...
                                  control_field,marker_field);
        [controls,markers,control_field,marker_field]=populate_search('RNAPcontrol');
        datamatN=normalize_driver(norm,datamatN,datameta,norm_method,controls,markers,...
                                  control_field,marker_field);
        [controls,markers,control_field,marker_field]=populate_search('pAsensitive');
        datamatN=normalize_driver(norm,datamatN,datameta,norm_method,controls,markers,...
                                  control_field,marker_field);
        [controls,markers,control_field,marker_field]=populate_search('pAresistant');
        datamatN=normalize_driver(norm,datamatN,datameta,norm_method,controls,markers,...
                                  control_field,marker_field);
end
end
%-----> Subroutines
function av = plate_average(vector,norm_method)
%--->
% determine plate average
%--->
switch norm_method
    case 'median'
        av=nanmedian(vector);
    case 'middlemean'
        p25=nan25percentile(vector);
        p75=nan75percentile(vector);
        range=(vector>p25 & vector<p75);
        av=nanmean(vector(range));
end
end    
function datamatN=normalize_driver(norm,datamat,datameta,norm_method,controls,markers,control_field,marker_field)
            fields=fieldnames(norm);
            for m=1:length(fields)
                datamatN.(fields{m})=NaN*datamat.(fields{m});
                datamat_size=size(norm.(fields{m}));
                for i=1:datamat_size(1)
                    for k=1:datameta.rep(i)
                        av=plate_average(norm.(fields{m})(i,:,k)',norm_method);
                        datamatN.(fields{m})(i,:,k)=strain_norm(av,controls,markers,...
                                                                datameta,...
                                                                datamat.(fields{m})(i,:,k),...
                                                                control_field,...
                                                                marker_field);
                    end
                end
            end
end
function vector=strain_norm(av,controls,markers,datameta,input,control_field,marker_field)
%--->
% normalize by strain
%--->
cont_fields=fieldnames(controls); %assume cont and mark have same fieldnames
                                  %use controls to create list.
vector=input;
for a=1:length(cont_fields)
    %find the indices of control and marker
    cnt_ind=false(length(input),1);
    for c=1:length(controls.(cont_fields{a}))
        temp_ind=strcmpi(controls.(cont_fields{a}){c},datameta.(control_field));
        cnt_ind=cnt_ind|temp_ind;
    end
    mrk_ind=false(length(input),1);
    for b=1:length(markers.(cont_fields{a}))
        tmp_ind=strcmpi(markers.(cont_fields{a}){b},datameta.(marker_field));
        mrk_ind=mrk_ind|tmp_ind;
    end
    %normalize marker to average
    cntav=nanmean(vector(cnt_ind));
    vector(mrk_ind)=vector(mrk_ind).*(av/cntav);
end
end
    

%find mean size of marker strain
%normalize all marker to this value