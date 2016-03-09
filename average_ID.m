function [uniqmat,uniqmeta] = average_ID(mat,meta,ID)
%-------------------------------------------------------------------------!
% [uniqmat,uniqmeta]=average_ID(mat,meta,fieldID)                         :       :
%-------------------------------------------------------------------------!
% average_ID.m                                                            :
% -------------                                                           :
% average_uid uses the specified field in metadata to average measurements: 
%  in                                                                     :
% the scoremat provided. Returned is a unique datamatrix and metadata     :
% structure.                                                              :
%-------------------------------------------------------------------------!
% Anthony Shiver (2014)                                                   :
%-------------------------------------------------------------------------!
[uid,un2deg,deg2un]=unique(meta.(ID));
[uniqmat,uniqmeta]=prepstructure(mat,meta,uid);
fields=fieldnames(mat);
for i=1:length(uid)
    deg_elem=(deg2un==i);
    uniq_elem=un2deg(i);
    for m=1:length(fields)
        uniqmat.(fields{m})(:,i,:)=nanmean(mat.(fields{m})(:,deg_elem,:),2);
        uniqmeta.([fields{m},'_pltstd'])(:,i,:)=nanstd(mat.(fields{m})(:,deg_elem,:),0,2);
    end
    uniqmeta.acc(i)=meta.acc(uniq_elem);
    uniqmeta.uid(i)=meta.uid(uniq_elem);
    uniqmeta.mut(i)=meta.mut(uniq_elem);
    uniqmeta.mrk(i)=meta.mrk(uniq_elem);
end
end
%--------------------------------------------------------------------------
function [uniqmat,uniqmeta]=prepstructure(mat,meta,uid)
%unchanged fields
uniqmeta.cnd=meta.cnd;
uniqmeta.cnc=meta.cnc;
uniqmeta.bch=meta.bch;
%reduced/dependent fields
fields=fieldnames(mat);
for m=1:length(fields)
    originalsize=size(mat.(fields{m}));
    %behave accordingly as to whether this is a raw or S-score matrix
    if(length(originalsize)>=3)
        uniqmat.(fields{m})=NaN*ones(originalsize(1),length(uid),originalsize(3));
    elseif(length(originalsize)<3)
        uniqmat.(fields{m})=NaN*ones(originalsize(1),length(uid));
    end
end
uniqmeta.acc=cell(length(uid),1);
uniqmeta.uid=cell(length(uid),1);
uniqmeta.mut=cell(length(uid),1);
uniqmeta.mrk=cell(length(uid),1);
end

