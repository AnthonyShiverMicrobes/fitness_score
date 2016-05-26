function [datamatS,datameta] = smooth_data(model,datamat,datameta)
%-------------------------------------------------------------------------!
% [datamatS,datameta]=smooth_data(model,datamat,datameta)                 :
%-------------------------------------------------------------------------!
% smooth_data.m fits a spatial model, generated earlier, to the fields    :
% contained in datamat. The fit is then used to multiplicatively smooth   :
% the data (spatially), centered on the median value for the parameter    :
% measured across the entire batch. Removed the step limiting values to be:
% above zero.                                                             :
%-------------------------------------------------------------------------!
% Anthony Shiver (2014)                                                   :
% 08022015: added outer row and column normalization                       :
%-------------------------------------------------------------------------!

% extract info for building data structures
fields=fieldnames(datamat);
datasize=size(datamat.(fields{1}));%take first field as representative
modelsize=size(model);

% extract info on outer rows and columns
outerrow=max(datameta.row);
outercol=max(datameta.col);
outer_index= (datameta.row <= 2 | datameta.row >= (outerrow-1)) | ...
             (datameta.col <= 2 | datameta.col >= (outercol-1)) ;

% create index for zero values in all fields
for m=1:length(fields)
    zeroindex.(fields{m})=datamat.(fields{m})==0;
end

% NaN the zeros to simplify model
for m=1:length(fields)
    datamat.(fields{m})(zeroindex.(fields{m}))=NaN;
end

% build skeletons
for m=1:length(fields)
    datamatS.(fields{m})=NaN*datamat.(fields{m});
    datameta.([fields{m},'_smt'])=NaN*ones(datasize(1),modelsize(2)+1,datasize(3));
    datameta.([fields{m},'_gene_model'])=NaN*ones(datasize(2));
    %for median of plate, ignore the outer rows and columns
    inner_for_average=datamat.(fields{m});
    inner_for_average(:,outer_index,:)=NaN;
    struct_middlemean.(fields{m})=nanmiddlemean(inner_for_average); %normalize data to middlemean across batch
    datamatN.(fields{m})=NaN*datamat.(fields{m});
end

% Build model for gene contribution by averaging all measurements of colony
for m=1:length(fields)
    %remove outer rows and columns for purpose of calculating normalization
    inner_for_average=datamat.(fields{m});
    inner_for_average(:,outer_index,:)=NaN;
    %normalize plate so that median is constant
    for i=1:datasize(1)
        for k=1:datasize(3)
            datamatN.(fields{m})(i,:,k)=datamat.(fields{m})(i,:,k)*(struct_middlemean.(fields{m})/nanmiddlemean(inner_for_average(i,:,k)));
        end
    end
    datameta.([fields{m},'_gene_model'])=nanmedian(nanmean(datamatN.(fields{m}),3),1)';
end


% smooth data  inner_for_average
for m = 1 : length(fields)
    for i = 1 : datasize(1)
        for k = 1 : datameta.rep(i)
            if(sum(~isnan(datamat.(fields{m})(i,:,k))))
                %remove outers normalize separately
                [outers]=extract_outers(datameta,datamatN.(fields{m})(i,:,k),outerrow,outercol);
                datamatN.(fields{m})(i,outer_index,k)=NaN;
                %remove gene contribution
                vector=datamatN.(fields{m})(i,:,k)'./datameta.([fields{m},'_gene_model']);
                %smooth interior
                datameta.([fields{m},'_smt'])(i,:,k)=robustfit(model,vector)';
                datamatS.(fields{m})(i,:,k)=(datamatN.(fields{m})(i,:,k)./(datameta.([fields{m},'_smt'])(i,:,k)*[ones(modelsize(1),1), model]'));            
                %reinsert outer rows and columns, but normalized
                datamatS.(fields{m})(i,:,k)=insert_outers(datamatS.(fields{m})(i,:,k),outers);
            end 
        end
    end
end
for m=1:length(fields)
%restore zeros and ensure no negative values
datamatS.(fields{m})(zeroindex.(fields{m}))=0;
datamatS.(fields{m})(datamatS.(fields{m})<0)=0;
end
end
%--subfunctions
function [outers]=extract_outers(datameta,vector,outerrow,outercol)
%--populate with index, then use index to extract data
outers.r1.ind=(datameta.row==1);
outers.r1.dat=vector(outers.r1.ind);
%
outers.r2.ind=(datameta.row==2 & ~(datameta.col==1) & ~(datameta.col==outercol));
outers.r2.dat=vector(outers.r2.ind);
%
outers.r3.ind=(datameta.row==(outerrow-1) & ~(datameta.col==1) & ~(datameta.col==outercol));
outers.r3.dat=vector(outers.r3.ind);
%
outers.r4.ind=(datameta.row==outerrow);
outers.r4.dat=vector(outers.r4.ind);
%
outers.c1.ind=(datameta.col==1 & ~(datameta.row==1) & ~(datameta.row==outerrow));
outers.c1.dat=vector(outers.c1.ind);
%
outers.c2.ind=(datameta.col==2 & ~(datameta.row<=2) & ~(datameta.row>=(outerrow-1)));
outers.c2.dat=vector(outers.c2.ind);
%
outers.c3.ind=(datameta.col==(outercol-1) & ~(datameta.row<=2) & ~(datameta.row>=(outerrow-1)));
outers.c3.dat=vector(outers.c3.ind);
%
outers.c4.ind=(datameta.col==outercol & ~(datameta.row==1) & ~(datameta.row==outerrow));
outers.c4.dat=vector(outers.c4.ind);
end
function [vector]=insert_outers(vector,outers)
norm=nanmiddlemean(vector(:));
chunk={'r1','r2','r3','r4','c1','c2','c3','c4'};
for i=1:length(chunk)
    vector(outers.(chunk{i}).ind)=outers.(chunk{i}).dat*(norm/nanmiddlemean(outers.(chunk{i}).dat(:)));
end
end
function average=nanmiddlemean(vector)
    up=nan75percentile(vector(:));
    lw=nan25percentile(vector(:));
    middle=vector( vector < up & vector > lw);
    average=nanmean(middle(:));
end


