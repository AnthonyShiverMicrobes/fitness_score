function [collmat,collmeta] = collapse_structure(mat,meta)
%-------------------------------------------------------------------------!
% [uniqmat,uniqmeta]=collapse_structure(mat,meta,fieldID)                 :
%-------------------------------------------------------------------------!
% collapse_structure.mat                                                  :
% -------------                                                           :
% collapses replicates according to a side by side arrangement, specifical:
% -ly written for 1536 format KEIO plates, please never use for anything  :
% else.                                                                   : 
%                                                                         :
%-------------------------------------------------------------------------!
% Anthony Shiver (2015)                                                   :
%-------------------------------------------------------------------------!
[collmat,collmeta]=prepstructure(mat,meta);
fields=fieldnames(mat);
for m=1:length(fields)
    S=size(mat.(fields{m}));
    SPLIT1=mat.(fields{m})(:,1:2:1536,:);
    SPLIT2=mat.(fields{m})(:,2:2:1536,:);
    for k=1:S(3)
        collmat.(fields{m})(:,:,(2*k-1))=SPLIT1(:,:,k);
        collmat.(fields{m})(:,:,2*k)=SPLIT2(:,:,k);
    end
end
end
%--------------------------------------------------------------------------
function [collmat,collmeta]=prepstructure(mat,meta)
%unchanged fields
collmeta.cnd=meta.cnd;
collmeta.cnc=meta.cnc;
collmeta.bch=meta.bch;
%reduced/dependent fields
fields=fieldnames(mat);
for m=1:length(fields)
    originalsize=size(mat.(fields{m}));
    collmat.(fields{m})=NaN*ones(originalsize(1),originalsize(2)/2,originalsize(3)*2); %assumes only doublets
end
collmeta.acc=meta.acc(1:2:1536);
collmeta.uid=meta.uid(1:2:1536);
collmeta.mut=meta.mut(1:2:1536);
collmeta.mrk=meta.mrk(1:2:1536);
end