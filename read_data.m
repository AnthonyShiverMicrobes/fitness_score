function [datamat,datameta]=read_data(path,array_meta,cond_meta,format)
%-------------------------------------------------------------------------!
% datamat=read_data(path,array_meta,cond_meta)                            :
%-------------------------------------------------------------------------!
% read_data builds a datamat structure using information from both the    :
% array (array_meta) and conditions (cond_meta). Attempts to take care of :
% all known possible formats of colony size files, using format as the    :
% control string.                                                         :
% CURRENTLY ASSUMES THAT DATA FILE FORMAT IS RELAYED as                   :
% row col ...                                                             :
% 1   1   ...                                                             :
% 1   2   ...                                                             :
% ...                                                                     :
% 1   48  ...                                                             :
% 2   1   ...                                                             :
% ...                                                                     :
% ...                                                                     :
%-------------------------------------------------------------------------!
% See slurp_dat.m for supported data formats                              :
%-------------------------------------------------------------------------!
% Anthony Shiver (2013)                                                   :
% 09-12-13: Split read data between datamat and datameta                  :
% 08-05-13: Added code to take cnd-cnc as unique for datamat population   :
%-------------------------------------------------------------------------!

%---
%estimate size of dataset, initialize the datamat
uni_cnd=strcat(cond_meta.cnd,cond_meta.cnc);
conditions=unique(uni_cnd);
depth=max(cellfun(@(x) sum(ismember(uni_cnd,x)),conditions));
[datameta,datamat]=initialize_datastructures([length(conditions),length(array_meta.row),depth],format); %fn @ EOF

%---
%tack on array information
datameta.row=array_meta.row;
datameta.col=array_meta.col;
datameta.plt=array_meta.plt;
datameta.wll=array_meta.wll;
datameta.acc=array_meta.acc;
datameta.uid=array_meta.uid;
datameta.mut=array_meta.mut;
datameta.mrk=array_meta.mrk;

%---
%fill datamat with condition-dependent information
for i = 1: length(cond_meta.nme)
        data_file=slurp_dat([path,filesep,char(cond_meta.nme(i))],format);
        [x1,datameta]=condition_update(cond_meta.cnd{i},cond_meta.cnc{i},cond_meta.bch{i},datameta); %fn @EOF, assigns values to .cnd,.cnc,.bch if new
        datamat=insert_new_datafile(datamat,data_file,x1,datameta.rep(x1)); %fn @ EOF
        datameta.src(x1,datameta.rep(x1))=cond_meta.src(i);
        datameta.flg(x1,datameta.rep(x1))=cond_meta.flg(i);
        file_elements=regexp(char(cond_meta.nme(i)),'\.','split','once');
        datameta.file{x1,datameta.rep(x1)}=file_elements{1};
        datameta.app{x1,datameta.rep(x1),:}=file_elements{2:length(file_elements)}; %uninitialized, no assumption of numel(file_elements)
end
end

%--------------------------------------------------------------------------
function [x1,datameta]=condition_update(condition,concentration,batch,datameta)
%for the passed condition-concentration, finds within datamat, updates replicate number
%and returns position of this condition. If not found, creates the
%condition and associated information
uni_cnd=strcat(datameta.cnd,datameta.cnc);
cndcnc=strcat(condition,concentration);
if(max(strcmp(cndcnc,uni_cnd))==1)
    x1=find(strcmp(cndcnc,uni_cnd));
    datameta.rep(x1)=datameta.rep(x1)+1;
else
    x1=length(datameta.cnd)+1;
    datameta.cnd{x1,1}=condition;
    datameta.cnc{x1,1}=concentration;
    datameta.bch{x1,1}=batch;
end
end

%--------------------------------------------------------------------------
function [datameta,datamat]=initialize_datastructures(size,format)
%creates datamat fields, taking into account the format of the
%datafile containing colony sizes. Leaves datamat.app uninitialized.
%initializes fields for datameta.
%---
%datamat, format-specific
switch format
    case 'collins_v1'
        datamat.sze=ones(size)*NaN;
        datamat.crc=ones(size)*NaN;
    case 'collins_v2'
        datamat.sze=ones(size)*NaN;
        datamat.crc=ones(size)*NaN;
        datamat.int=ones(size)*NaN;
    case 'iris_v0_ecogrowth'
        datamat.sze=ones(size)*NaN;
    case 'iris_v0_ecoopacity'
        datamat.sze=ones(size)*NaN;
        datamat.crc=ones(size)*NaN;
        datamat.opc=ones(size)*NaN;
    case 'iris_kritikos'
        datamat.opc=ones(size)*NaN;
    otherwise
        disp([filepath ': Proper format not specified!']);
end
%---
%datamat, format-independent
datameta.rep=ones(size(1),1);
datameta.cnd={}; %.cnd is uninitialized, because it will be added piece by piece
datameta.cnc={}; %.cnc is uninitialized, and will be added piece by piece
datameta.bch=cell(size(1),1);
datameta.src=cell(size(1),size(3));
datameta.flg=cell(size(1),size(3));
datameta.file=cell(size(1),size(3));
end

%--------------------------------------------------------------------------
function datamat=insert_new_datafile(datamat,datafile,row,replicate)
% adds fields in data file to appropriate location in datamat
% independently of the specific fields passed in datafile
fields=fieldnames(datafile);
for i=1:length(fields)
    if isfield(datamat,fields{i})
        datamat.(fields{i})(row,:,replicate)=datafile.(fields{i})';
    else
        disp('Error: Field found in file but not in datamat incorrect format assumed!');
    end
end
end

    