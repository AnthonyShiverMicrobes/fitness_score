function [score,meta,scoremat,datameta,datamatV,datamatN,datamatS,...
          datamatF,datamat,spatialmodel,settings]...
          =scoreRNAP(cond_file,array_file,data_path,file_format,ignore,bad)
%-------------------------------------------------------------------!
%[xxx]=scoreRNAP('condition file','array file','data path','format','ignore','bad')
%ALS 2013.
%-------------------------------------------------------------------!
%store passed variables
settings.cond_file=cond_file;
settings.array_file=array_file;
settings.data_path=data_path;
settings.file_format=file_format;
settings.ignore=ignore;
settings.bad=bad;
%define variables
settings.frac=0.25; %fraction of colony sizes that can be zero without before filtering colony
settings.variance_limit=1000; %limit of colony size variance before filtering S-score
settings.spatial_model_type='quartic'; %type of spatial model to use
settings.average_method='middlemean'; %how to estimate average
settings.variance_method='mad'; %how to estimate variance
settings.numDup=1; %no collapsing
settings.normalization_rnap='RNAP_marker';
settings.platetype='rnap';
settings.normalize_field='mut';
settings.plate_average_field='uid';
settings.power_transform=false;
settings.variance_scale='RC4';
%read data, populate datastructures
[datamat,datameta]=read_data(settings.data_path,...
                            read_array_key(settings.array_file),...
                            read_condition_key(settings.cond_file),...
                            settings.file_format);
spatialmodel=generate_model(datameta.col,datameta.row,settings.spatial_model_type);
%-- filter
[datamatF,datameta]=filter_data(datamat,datameta,settings.platetype,settings.frac);
%-- remove bad data
datamatF=remove_bad_strains(datamatF,settings.bad);
%-- smooth data
[datamatS,datameta]=smooth_data(spatialmodel,datamatF,datameta);
%-- squeeze outliers
datamatS=squeeze_outliers(datamatS,datameta);
datamatS=nomorezeros(datamatS);
%-- power transform data
if(settings.power_transform==true)
    datamatP=transform_data(datamatS,0.5);
else
    datamatP=datamatS;
end
%-- normalize size
[datamatN,datameta]=normalize_data(datamatP,datameta,settings.average_method,settings.ignore,settings.normalize_field);
datamatN=normalize_strain(datamatN,datameta,settings.normalization_rnap,settings.average_method,settings.ignore,settings.normalize_field);
datamatN=squeeze_outliers(datamatN,datameta);
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
%-- remove low replicates
datamatV=enforce_triplicates(datamatV);
%-- adapt for toolbox
datamatVscaled=fit_to_fivehundred(datamatV,datameta,settings.average_method,settings.numDup);
datamatFscaled=fit_to_fivehundred(datamatF,datameta,settings.average_method,settings.numDup);
datameta=controlstats(datamatFscaled,datameta,'F'); 
datameta=controlstats(datamatVscaled,datameta,'V');
%score the data
scoremat=calculate_score(datamatFscaled,datamatVscaled,datameta,settings.variance_limit,settings.numDup);
%average degenerate sites
[score,meta]=average_ID(scoremat,datameta,settings.plate_average_field);
