function datamatV = RC4_scale(datamat,datameta,methodcenter,methodscale)
%--
%
%
%--
fields=fieldnames(datamat);

for m=1:length(fields)
    datamatV.(fields{m})=datamat.(fields{m}).*NaN;
    S=size(datamat.(fields{m}));
    for i = 1:S(1)
        for k = 1:S(3)
            %step through the five groups, center and scale them
            vectorik=datamat.(fields{m})(i,:,k);
            scaledvector=scalesubsets(vectorik,methodcenter,methodscale,datameta);
            vik=scaleplate(scaledvector,methodcenter,methodscale);
            datamatV.(fields{m})(i,:,k)=vik;
        end
    end
    minvar=min(datamatV.(fields{m})(:));
    datamatV.(fields{m})=datamatV.(fields{m})+abs(minvar)+1;
end
end

function vscaled=scalesubsets(v,methodcenter,methodscale,meta)
%
%
vscaled=v*NaN;
%build indices
row=meta.row;
maxrow=max(row);
col=meta.col;
maxcol=max(col);
rc1i=(row==1)|(col==1)|(row==maxrow)|(col==maxcol);
rc2i=(row==2&(col>1&col<maxcol))|(row==maxrow-1&(col>1&col<maxcol))|...
     ((row>1&row<maxrow)&col==2)|((row>1&row<maxrow)&col==maxcol-1);
rc3i=(row==3&(col>2&col<maxcol-1))|(row==maxrow-2&(col>2&col<maxcol-1))|...
     ((row>2&row<maxrow-1)&col==3)|((row>2&row<maxrow-1)&col==maxcol-2);
rc4i=(row==4&(col>3&col<maxcol-2))|(row==maxrow-3&(col>3&col<maxcol-2))|...
     ((row>3&row<maxrow-2)&col==4)|((row>3&row<maxrow-2)&col==maxcol-3);
inneri=(row>4&row<maxrow-3)&(col>4&col<maxcol-3);
%grab subsets
set.r1=v(rc1i);
set.r2=v(rc2i);
set.r3=v(rc3i);
set.r4=v(rc4i);
set.in=v(inneri);
%--- center and scale
fields=fieldnames(set);
%center
for a=1:length(fields)
    switch methodcenter
        case 'middlemean'
            center.(fields{a})=middlemean(set.(fields{a})(:));
        case 'median'
            center.(fields{a})=nanmedian(set.(fields{a})(:));
        otherwise
            center.(fields{a})=nanmedian(set.(fields{a})(:));
    end
    set.(fields{a})=set.(fields{a})-center.(fields{a});
end
%scale
for a=1:length(fields)
    switch methodscale
        case 'mad'
            scale.(fields{a})=medianabsolutedeviation(set.(fields{a})(:));
        case 'iqr'
            scale.(fields{a})=nan75percentile(set.(fields{a})(:))-...
                              nan25percentile(set.(fields{a})(:));
        otherwise
            scale.(fields{a})=nan75percentile(set.(fields{a})(:))-...
                              nan25percentile(set.(fields{a})(:));
    end
    %prevent creation of Inf values if majority zero
    if scale.(fields{a})~=0
        set.(fields{a})=set.(fields{a})./scale.(fields{a});
    else
        set.(fields{a})=set.(fields{a}).*NaN;
    end
    %reconstruct
    vscaled(rc1i)=set.r1;
    vscaled(rc2i)=set.r2;
    vscaled(rc3i)=set.r3;
    vscaled(rc4i)=set.r4;
    vscaled(inneri)=set.in;
end
end
function vscaled=scaleplate(v,methodcenter,methodscale)
%
%
switch methodcenter
    case 'median'
        center=nanmedian(v(:));
    case 'middlemean'
        center=middlemean(v(:));
    otherwise
        center=nanmedian(v(:));
end
switch methodscale
    case 'iqr'
        scale=nan75percentile(v(:))-nan25percentile(v(:));
    case 'mad'
        scale=medianabsolutedeviation(v(:));
    otherwise 
        scale=nan75percentile(v(:))-nan25percentile(v(:));
end
vscaled=(v-center)./scale;
end
function scale=medianabsolutedeviation(list)
list=abs(list-nanmedian(list));
scale=nanmedian(list(:));
end
function center=middlemean(list)
    ind=list<nan75percentile(list(:))&list>nan25percentile(list(:));
    list=list(ind);
    center=nanmean(list(:));
end
