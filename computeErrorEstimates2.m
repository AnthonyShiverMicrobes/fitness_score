function err = computeErrorEstimates2(raw,rawn)
%computes minimum bounds for the standard error as a function of the
%unnormalized and normalized colony sizes
%written by Sean Collins (2006) as part of the EMAP toolbox
%modified version to adapt to possibly split datasets

[r c n] = size(raw.size);

binsize=50;

natMax=max(1300,max(round(max(raw.meansize))));

kanMax=max(1300,max(round(max(rawn.meansize))));

nats = zeros(r,c,n)*NaN;
for k=1:n
    means(:,:,k)=rawn.meansize;
end
error=rawn.size-means;
relerr=error./means;

fn=fieldnames(rawn);
if ismember('numDup',fn)
    numDup=max(rawn.numDup(:));
    if numDup>2 || numDup<1
        fprintf('Unusual value for number of duplicates: %i\n',numDup);
    end
else
    numDup=2;
end

for i=1:r
    for k=1:numDup:n
        plate=data2plate2(raw.size(i,:,k:(k+numDup-1)),raw.row,raw.col);
        [r1 c1]=size(plate);
        nats(i,k)=myNanMedian(plate(2:(r1-1),2:(c1-1)));  %the typical unnormalized colony size on the plate
        nats(i,k+numDup-1)=nats(i,k);   %trivial if numDup==1
        nats2(i,1:c,(k:(k+numDup-1)))=nats(i,k);
    end
end
for j=1:c
    temp = rawn.size(:,j,:);
    kans(j)=myNanMedian(temp(:));   %the typical normalized size for a given KAN marked strain
    kans2(1:r,j,1:n)=kans(j);
end

%computing error estimates as a function of the phenotype of the NAT marked
%strain (the typical unnormalized size)
ind=find(~isnan(nats2));
natSizes=nats2(ind); natSizes=natSizes(:);
relErrors=relerr(ind);relErrors=relErrors(:);
m2=means(ind);m2=m2(:);
mat=[natSizes relErrors m2];
sortedNats=sortrows(mat,1);
b1=1; b2=1;
for i=10:natMax
    b1=findmin(sortedNats(:,1),i,binsize,b1);
    b2=findmax(sortedNats(:,1),i,binsize,b2);
    list=sortedNats(b1:b2,2);
    bynat(i)=myNanSD(list);
    natm(i)=myNanMean(sortedNats(b1:b2,3));
end
err.natm=natm;
sn2=round(sortedNats); sn2=sn2(:,1);
counts=hist(sn2,0:natMax);
ind1=~isnan(bynat);
tot = sum(counts(ind1));
w = counts/tot;
wmean = sum(w(ind1).*bynat(ind1));  %This is a weighted mean of the relative error measurements
bynat = bynat/wmean; %The nat data will be used as a correction to the kan-based estimate, so we compute the ratio of the observed relative error to the weighted mean to see if we are noisier or less noisy than average
bn1=bynat;
bynat = -1*(pav(-1*bynat(:)));     %This is to smooth the nat estimate a bit (making it monotonic)

%computing relative error as a function of the phenotype of the KAN marked
%strain
ind=find(~isnan(kans2));
k2=kans2(ind);k2=k2(:);
e2=relerr(ind);e2=e2(:);
e3=error(ind);e3=e3(:);
mat=[k2 e2 k2 e3];
sortedKans=sortrows(mat,1);

b1=1; b2=1;
for i=10:kanMax
    b1=findmin(sortedKans(:,1),i,binsize,b1);
    b2=findmax(sortedKans(:,1),i,binsize,b2);
    list=sortedKans(b1:b2,2);
    kanm(i)=myNanMean(sortedKans(b1:b2,3));
    bykan(i)=myNanSD(list);                 %based on the relative error
    bykan2(i)=myNanSD(sortedKans(b1:b2,4)); %based on the absolute error
end
err.kanm=kanm;
err.bykan2=bykan2;

%computing expected errors
err.exp=zeros(r,c,n)*NaN;
for i=1:r
    for j=1:c
        for k=1:n
            nat=round(nats2(i,j,k));
            kan=round(kans2(i,j,k));
            if isnan(nat) || isnan(kan) || nat==0 || kan==0
                err.exp(i,j,k)=NaN;
            else
                err.exp(i,j,k)=bynat(nat)*bykan(kan)*means(i,j,k);
            end
        end
    end
end

err.binsize=binsize;
err.bynat=bynat;
err.bykan=bykan;
end

%*************************************************************************

function ind=findmin(list,middle,binsize,prev)
target=middle-binsize/2;
ind=prev;
while (list(ind)<target)&(ind<length(list))
    ind=ind+1;
end
end

%*************************************************************************

function ind=findmax(list,middle,binsize,prev)
target=middle+binsize/2;
ind=prev;
while (list(ind)<target)&(ind<length(list))
    ind=ind+1;
end
if (ind>1)
    ind=ind-1;
end
end
