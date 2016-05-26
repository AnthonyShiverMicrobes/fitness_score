function plate = data2plate(dat)
%DATA2PLATE takes a 1x384x2 data matrix and converts it to a 32x48 matrix
%that represents physical positions on the plate. Half of the values will
%be empty (NaN) because of the staggered alignment of the colonies.
%
%written by Sean Collins (2006) as part of the EMAP toolbox

[r c n]=size(dat);
if r==1 && c==384 && n==2
    plate = zeros(32,48) * NaN;
    for j = 1:384
        for k=1:2
            row = index2row(j,k);
            col = index2col(j,k);
            plate(row,col)=dat(1,j,k);
        end
    end
else
    if n==2 || r~=1
        %error('Inconsistent colony matrix dimension in the normalization procedure');
        %fprintf('Not properly reconstructing plate organization\n');
        r1=round(sqrt(2*r*c*n/3)); c1=round(r*c*n/r1);
        plate=reshape(dat,[c1 r1])';
        return;
    elseif n==4
        numCol=c*n;
        switch numCol
            case 96
                r1=8; c1=12;
            case 384
                r1=16; c1=24;
            case 1536
                r1=32; c1=48;
        end
        v1=ones(1,c*n/2)*NaN;
        v2=ones(1,c*n/2)*NaN;
        v1(1:2:length(v1))=dat(1,:,1);
        v1(2:2:length(v1))=dat(1,:,2);
        v2(1:2:length(v1))=dat(1,:,3);
        v2(2:2:length(v1))=dat(1,:,4);
        p1=reshape(v1,[c1 r1/2])';
        p2=reshape(v2,[c1 r1/2])';
        plate=zeros(r1,c1)*NaN;
        for i=1:(r1/2)
            plate(i*2-1,1:c1)=p1(i,:);
            plate(i*2,1:c1)=p2(i,:);
        end
    else
        numCol=c;
        switch numCol
            case 96
                r=8; c=12;
            case 384
                r=16; c=24;
            case 1536
                r=32; c=48;
        end
        plate=reshape(dat,[c r])';
    end
end
