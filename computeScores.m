function scoremat = computeScores(raw,err,method)
%this function converts the raw colony size data into unaveraged S scores.
%It can use one of two methods for computing the expected colony size.
%
%written by Sean Collins (2006) as part of the EMAP toolbox

[r c n] = size(raw.size);
scoremat.rowlabels = raw.rowlabels;
scoremat.collabels = raw.collabels;
scoremat.rep_num = raw.rep_num ;
fn=fieldnames(raw);
if ismember('row_mut',fn)
    scoremat.row_mut = raw.row_mut ;
end
if ismember('col_mut',fn)
    scoremat.col_mut = raw.col_mut ;
end

% Compute the expected colony size and standard error values
for j=1:c
    SD(j)=myNanMedian(raw.sdsize(:,j));
end
if strcmp(method,'parzen')
    
    for j=1:c
        
        %Estimate required window size from median of median of plates
        %Should be what you normalized to earlier. Take the median to be 
        %one-tenth of this value.
        
        window = round( nanmedian( nanmedian( raw.size( :, :, 1 ), 2 ), 1 ) / 10 ) ;
        
        m( j ) = estimateCenterParzenWindow( raw.size( :, j, : ), window ) ;
    
    end

else
    %method = 'median' is the default
    for j=1:c
        m(j)=myNanMedian(raw.size(:,j,:));
    end
end

% Compute the epistasis scores
for i=1:r
    for j=1:c
        if ~strcmp(scoremat.collabels(j),'N/A')
            scoremat.data(i,j)=scoreS(raw,i,j,err,m(j),SD(j));
        else
            scoremat.data(i,j)=NaN;
        end
    end
end

% Associate gene names
fn=fieldnames(raw);
if ismember('geneToOrf',fn)
    scoremat.geneToOrf=raw.geneToOrf;
end

