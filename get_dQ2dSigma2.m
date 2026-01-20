function [d2QdSigma2] = get_dQ2dSigma2 ( TSigma  , J2 ) 
txx = TSigma(1) ; 
tyy = TSigma(2) ; 
txy = TSigma(3) ; 
tzz = TSigma(4) ; 
d2Qdsxxdsxx = 1 .* (1 .* J2) .^ (-1 ./ 2) ./ 3 - 1 .^ 2 .* txx .^ 2 .* (1 .* J2) .^ (-3 ./ 2) ./ 4;
d2Qdsxxdsyy = -1 .* (1 .* J2) .^ (-1 ./ 2) ./ 6 - 1 .^ 2 .* txx .* tyy .* (1 .* J2) .^ (-3 ./ 2) ./ 4;
d2Qdsxxdsxy = -1 .^ 2 .* txx .* txy .* (1 .* J2) .^ (-3 ./ 2) ./ 2;
d2Qdsxxdszz = -1 .* (1 .* J2) .^ (-1 ./ 2) ./ 6 - 1 .^ 2 .* txx .* tzz .* (1 .* J2) .^ (-3 ./ 2) ./ 4;
d2Qdsyydsxx = -1 .* (1 .* J2) .^ (-1 ./ 2) ./ 6 - 1 .^ 2 .* txx .* tyy .* (1 .* J2) .^ (-3 ./ 2) ./ 4;
d2Qdsyydsyy = 1 .* (1 .* J2) .^ (-1 ./ 2) ./ 3 - 1 .^ 2 .* tyy .^ 2 .* (1 .* J2) .^ (-3 ./ 2) ./ 4;
d2Qdsyydsxy = -1 .^ 2 .* tyy .* txy .* (1 .* J2) .^ (-3 ./ 2) ./ 2;
d2Qdsyydszz = -1 .* (1 .* J2) .^ (-1 ./ 2) ./ 6 - 1 .^ 2 .* tyy .* tzz .* (1 .* J2) .^ (-3 ./ 2) ./ 4;
d2Qdsxydsxx = -1 .^ 2 .* txx .* txy .* (1 .* J2) .^ (-3 ./ 2) ./ 2;
d2Qdsxydsyy = -1 .^ 2 .* tyy .* txy .* (1 .* J2) .^ (-3 ./ 2) ./ 2;
d2Qdsxydsxy = 1 .* (1 .* J2) .^ (-1 ./ 2) - 1 .^ 2 .* txy .^ 2 .* (1 .* J2) .^ (-3 ./ 2);
d2Qdsxydszz = -1 .^ 2 .* txy .* tzz .* (1 .* J2) .^ (-3 ./ 2) ./ 2;
d2Qdszzdsxx = -1 .* (1 .* J2) .^ (-1 ./ 2) ./ 6 - 1 .^ 2 .* txx .* tzz .* (1 .* J2) .^ (-3 ./ 2) ./ 4;
d2Qdszzdsyy = -1 .* (1 .* J2) .^ (-1 ./ 2) ./ 6 - 1 .^ 2 .* tyy .* tzz .* (1 .* J2) .^ (-3 ./ 2) ./ 4;
d2Qdszzdsxy = -1 .^ 2 .* txy .* tzz .* (1 .* J2) .^ (-3 ./ 2) ./ 2;
d2Qdszzdszz = 1 .* (1 .* J2) .^ (-1 ./ 2) ./ 3 - 1 .^ 2 .* tzz .^ 2 .* (1 .* J2) .^ (-3 ./ 2) ./ 4;
    
    d2QdSigma2 =  [
                    d2Qdsxxdsxx d2Qdsxxdsyy d2Qdsxxdsxy d2Qdsxxdszz  
                    d2Qdsyydsxx d2Qdsyydsyy d2Qdsyydsxy d2Qdsyydszz  
                    d2Qdsxydsxx d2Qdsxydsyy d2Qdsxydsxy d2Qdsxydszz  
                    d2Qdszzdsxx d2Qdszzdsyy d2Qdszzdsxy d2Qdszzdszz  
                    ]; 

end