function y=linearInterp(x1,y1,guide)
%Linearly interpolates from (x1,y1) onto the new set of x-values guide. x1
%and guide should be sorted lists with max(guide)<=max(x1) and
%min(guide)>=min(x1)

if min(guide)<min(x1) | max(guide)>max(x1)
    fprintf('Interpolation limit out of bounds -- Interpolation not performed\n');
    x=x1;
    y=y1;
else
    guide=sort(guide);
    temp=[x1(:) y1(:)];
    temp=sortrows(temp,1);
    x1=temp(:,1); y1=temp(:,2);
    for i=1:length(guide)
        c=x1-guide(i);
        d=abs(c);
        temp1=min(d(c>=0)); temp1=temp1(1);
        temp2=min(d(c<=0)); temp2=temp2(1);
        ind1=find((c>=0) & d==temp1); ind1=ind1(1);
        ind2=find((c<=0) & d==temp2); ind2=ind2(1);
        diff=x1(ind1)-x1(ind2);
        x(i)=guide(i);
        if diff>0
            frac=(x(i)-x1(ind2))/diff;
            y(i)=y1(ind1)*(frac) + y1(ind2)*(1-frac);
        else
            y(i)=y1(ind1);
        end
    end
end
        
