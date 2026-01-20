function [local] = element_coordinate (POINT , sctr , node )
% this function gives element's local coordinate
% inputs: 
% POINT : global coordinate of points
% sctr : element's scatter vector 
% node : coordinate of all nodes
local = zeros(size(POINT,1),2) ; 
for K = 1 : size(POINT,1)
point = POINT(K,:) ; 
X = node(sctr,1) ; 
Y = node(sctr,2) ; 

if size(sctr,2) == 4
    elemType = 'Q4';
else
    elemType = 'T3' ;
end

xi = 0 ;    eta = 0 ; 
% xp = 9.77 ;
% yp = 9.8 ;
xp = point(1); 
yp = point(2); 
normmm = 19 ; 
itr = 0;  
while(normmm)>1e-8
    itr = itr + 1;  
    pt = [xi , eta ] ;    
    [N,dNdx]  = lagrange_basis(elemType,pt);

    f(1) = xp - N' * X ;
    f(2) = yp - N' * Y ;

    r = f' ;
    df(1,1) = - dNdx(:,1)' * X ; 
    df(1,2) = - dNdx(:,2)' * X ;
    df(2,1) = - dNdx(:,1)' * Y ; 
    df(2,2) = - dNdx(:,2)' * Y ;    
%     df
% pause
    J = - df ; 
    dd = J\r ;
    dxi = dd(1) ; 
    deta = dd(2) ; 
    xi = xi + dxi ;
    eta = eta + deta ; 
    normmm = norm(dd);
    if ( itr >  15 ) 
        itr
        error ( 'something went wrong') 
    end
end

local(K,:) = [xi eta] ; 
% local
% pause
end