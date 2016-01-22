function [score,meta,scoremat,datameta,datamatV,datamatN,datamatS,...
          datamatF,datamat,model]...
          =scoreRNAP(cond_file,array_file,data_path,format,ignore,bad)
%[xxx]=scoreRNAP('condition file','array file','data path','format','ignore','bad')
%ALS 2013.
%read data, populate datastructures
[datamat,datameta]=read_data(data_path,...
                            read_array_key(array_file),...
                            read_condition_key(cond_file),...
                            format);
model=generate_model(datameta.col,datameta.row,'quartic');
%filter and transform
datamatF=filter_data(datamat,datameta,'rnap');
%datamatP=transform_data(datamatF,0.5);
datamatP=datamatF;
%remove bad strains
datamatP=remove_bad_strains(datamatP,bad);
%smooth data
[datamatS,datameta]=smooth_data(model,datamatP,datameta);
%remove low replicates
datamatS=enforce_triplicates(datamatS);
%normalize size and variance
[datamatN,datameta]=normalize_data(datamatS,datameta,'middlemean',ignore,'mut');
datamatN=normalize_strain(datamatN,datameta,'RNAP_marker','middlemean',ignore,'mut');
%datamatV=scale_data(datamatN,'mad');
datamatV=datamatN;
%adapt for toolbox
datamatVscaled=fit_to_fivehundred(datamatV,datameta,'middlemean');datamatPscaled=fit_to_fivehundred(datamatP,datameta,'middlemean');
datameta=controlstats(datamatPscaled,datameta,'F'); datameta=controlstats(datamatVscaled,datameta,'V');
%score the data
scoremat=calculate_score(datamatPscaled,datamatVscaled,datameta);
%average degenerate sites
[score,meta]=average_ID(scoremat,datameta,'uid');
