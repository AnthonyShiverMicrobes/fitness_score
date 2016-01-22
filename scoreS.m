function t = scoreS(raw,i,j,err,m1,SDwt)
%computes an unaveraged S score for the pair of genes corresponding to row
%i and column j in the colony size matrix. It also takes the input variable
%err which is used to place a minimum bound on experimental variance, and
%it takes parameters m1 (the expected colony size) and SDwt (the estimated
%standard deviation for the control sample).
%
%written by Sean Collins (2006) as part of the EMAP toolbox

ind = find( ~isnan( raw.size( i, j, : ) ) ) ;

exper = raw.size( i, j, ind ) ;

exper = exper( : ) ;

experl = length( exper ) ;

thiserr = err.exp( i, j, : ) ;

wtl = 6 ;

%Minimum bound on experimental SD

if ( experl > 1 )
    
    varVect = thiserr( ind ) .* thiserr( ind ) ;
    
    var2 = sum( varVect ) / ( experl - 1 ) ;
    
    var2 = max( var2, var( exper ) ) ;

else

    var2 = NaN ;

end

%Minimum bound on SDwt

if ( SDwt < 0.064 * m1 )
    
    SDwt = 0.064 * m1 ;

end

var1 = SDwt * SDwt ;

if ( ( wtl > 1 ) && ( experl > 1 ) )
    
    svar = ( var1 * ( wtl - 1 ) + var2 * ( experl - 1 ) ) / ( wtl + experl - 2 ) ;
    
    m2 = mean( exper ) ;
    
    t = ( m2 - m1 ) / sqrt( svar * ( 1 / wtl + 1 / experl ) ) ;
    
else
    
    t = NaN;

end
