function [scoremat,datametaC,datamatVC,datameta,datamatV,datamatN,datamatS,...
          datamatF,datamat,settings]...
          =scorePHEN(cond_file,array_file,data_path,file_format,ignore,bad,plate_type,keio6ind)
%[xxx]=read_in_stats('condition file','array file','data path','format','ignore','bad','type')
%--ALS 2015.
%--Store passed variables in settings strucure 
settings.cond_file=cond_file;
settings.array_file=array_file;
settings.data_path=data_path;
settings.file_format=file_format;
settings.ignore=ignore;
settings.bad=bad;
settings.plate_type=plate_type;
settings.keio6ind=keio6ind;
%--HARD CODED VALUES
settings.frac=0.25; %fraction of colony sizes that can be zero without before filtering colony
settings.variance_limit=1000; %limit of colony size variance before filtering S-score
settings.spatial_model_type='quartic'; %type of spatial model to use
settings.average_method='middlemean'; %how to estimate average
settings.variance_method='mad'; %how to estimate variance
settings.numDup=2; %in the end, replicates will be collapsed
settings.eliminate_disagreement_field='acc';
settings.normalize_field='mut';
settings.power_transform=true;
settings.variance_scale='RC4';
%--read data, populate datastructures
[datamat,datameta]=read_data(settings.data_path,...
                            read_array_key(settings.array_file),...
                            read_condition_key(settings.cond_file),...
                            settings.file_format);
settings.spatialmodel=generate_model(datameta.col,datameta.row,settings.spatial_model_type);
%-- filter
[datamatF,datameta]=filter_data(datamat,datameta,settings.plate_type,settings.frac);
%-- remove bad data
datamatF=remove_bad_strains(datamatF,settings.bad);
%-- smooth data
[datamatS,datameta]=smooth_data(settings.spatialmodel,datamatF,datameta);
%-- squeeze outliers and eliminate disagreements
datamatS=squeeze_outliers(datamatS,datameta);
datamatS=eliminate_disagreement(datamatS,datameta,settings.eliminate_disagreement_field);
datamatS=nomorezeros(datamatS);
%-- power transform data
if(settings.power_transform)
    datamatP=transform_data(datamatS,0.5);
else
    datamatP=datamatS;
end
%-- normalize size
if(strcmpi('keio6',settings.plate_type))
    [datamatN,datameta]=normalize_data_split(datamatP,datameta,settings.average_method,settings.ignore,settings.normalize_field,settings.keio6ind);
else
    [datamatN,datameta]=normalize_data(datamatP,datameta,settings.average_method,settings.ignore,settings.normalize_field);
end
%-- scale variance of data
switch settings.variance_scale
    case 'RC4'
        datamatV=RC4_scale(datamatN,datameta,settings.average_method,settings.variance_method);
    case 'old'
        datamatV=scale_data(datamatN,'mad');
    case 'none'
        datamatV=datamatN;
    otherwise
        datamatV=datamatN;
end
%-- collapse replicates for one score
[datamatVC,datametaC]=collapse_structure(datamatV,datameta);
[datamatFC,~]=collapse_structure(datamatF,datameta);
%-- enforce at least three measurements
datamatVC=enforce_triplicates(datamatVC);
datamatFC=enforce_triplicates(datamatFC);
%-- adapt for toolbox
datamatVscaled=fit_to_fivehundred(datamatVC,datametaC,settings.average_method,settings.numDup); %final data
datamatFscaled=fit_to_fivehundred(datamatFC,datametaC,settings.average_method,settings.numDup); %original data
datametaC=controlstats(datamatFscaled,datametaC,'F'); 
datametaC=controlstats(datamatVscaled,datametaC,'V');
%-- score the data
scoremat=calculate_score(datamatFscaled,datamatVscaled,datametaC,settings.variance_limit,settings.numDup);