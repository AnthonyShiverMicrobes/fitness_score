function [score,meta,scoremat,datameta,datamatV,datamatN,datamatS,...
          datamatF,datamat,spatialmodel]...
          =scoreKEIO(cond_file,array_file,data_path,format,ignore,bad,type,keio6ind)
%[xxx]=read_in_stats('condition file','array file','data path','format','ignore','bad','type')
%--ALS 2015.
%--HARD CODED VALUES
frac=0.25; %fraction of colony sizes that can be zero without before filtering colony
variance_limit=1000; %limit of colony size variance before filtering S-score
spatial_model_type='quartic'; %type of spatial model to use
average_method='middlemean'; %how to estimate average
variance_method='mad'; %how to estimate variance
uniq_field='acc'; %how to compare strains
numDup=1; %no collapsing
%--read data, populate datastructures
[datamat,datameta]=read_data(data_path,...
                            read_array_key(array_file),...
                            read_condition_key(cond_file),...
                            format);
spatialmodel=generate_model(datameta.col,datameta.row,spatial_model_type);
%-- filter
[datamatF,datameta]=filter_data(datamat,datameta,type,frac);
%-- remove bad data
datamatF=remove_bad_strains(datamatF,bad);
%-- smooth data
[datamatS,datameta]=smooth_data(spatialmodel,datamatF,datameta);
%-- squeeze outliers and eliminate disagreements
datamatS=squeeze_outliers(datamatS,datameta);
datamatS=eliminate_disagreement(datamatS,datameta,uniq_field);
%-- power transform data
%datamatP=transform_data(datamatS,0.5);
datamatP=datamatS;
%-- normalize size
if(strcmpi('keio6',type))
    [datamatN,datameta]=normalize_data_split(datamatP,datameta,average_method,ignore,'mut',keio6ind);
else
    [datamatN,datameta]=normalize_data(datamatP,datameta,average_method,ignore,'mut');
end
%-- scale variance of data
%datamatV=scale_data(datamatN,'mad');
%datamatV=kritikos_scale(datamatN,datameta,average_method,variance_method);
datamatV=datamatN;
%-- enforce at least three measurements
datamatV=enforce_triplicates(datamatV);
score=[];scoremat=[];meta=[];
%-- adapt for toolbox
%datamatVscaled=fit_to_fivehundred(datamatV,datameta,average_method,numDup);
%datamatFscaled=fit_to_fivehundred(datamatF,datameta,average_method,numDup);
%datameta=controlstats(datamatFscaled,datameta,'F'); 
%datameta=controlstats(datamatVscaled,datameta,'V');
%-- score the data
%scoremat=calculate_score(datamatFscaled,datamatVscaled,datameta,variance_limit,numDup);
%-- average degenerate sites
%[score,meta]=average_ID(scoremat,datameta,uniq_field);
