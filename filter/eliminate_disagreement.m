function datamat=eliminate_disagreement(datamat,datameta,ID)
%-------------------------------------------------------------------------!
% datamat=eliminate_disagreement(datamat,datameta)                        :
% ------------------------------------------------                        :
% eliminate_disagreement.mat finds replicates where one strain has a zero :
% value and the second is normal, and NaNs both values to eliminate       :
% spurious zeros.Assumes two unique element (specific for KEIO layout     :
%-------------------------------------------------------------------------!
% Anthony Shiver 08-08-15                                                 :
%-------------------------------------------------------------------------!
fields=fieldnames(datamat);
[uid,~,deg2un]=unique(datameta.(ID));
for m=1:length(fields)
    cmd=nanmedian(nanmean(datamat.(fields{m}),3),1);
    for i=1:length(datameta.rep)
        for k=1:datameta.rep(i)
            %now operating plate by plate, look through unique IDs
            for a=1:length(uid)
                deg_elem=(deg2un==a);%index of the elements corresponding to uid(a)
                rel_values=(datamat.(fields{m})(i,deg_elem,k))./cmd(deg_elem);
                if( ((rel_values(1)==0) && (rel_values(2)>0.6)) || ...
                     ((rel_values(1)>0.6) && (rel_values(2)==0)) ) %one or other is zero
                    datamat.(fields{m})(i,deg_elem,k)=NaN;
                end
            end
        end
    end
end