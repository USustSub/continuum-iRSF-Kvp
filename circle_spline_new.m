% *************************************************************************
%                 TWO DIMENSIONAL ELEMENT FREE GALERKIN CODE
%                            Nguyen Vinh Phu
%                        LTDS, ENISE, Juillet 2006
% *************************************************************************
function [w,dwdx,dwdy,dwdxx,dwdyy,dwdxy] = circle_spline_new(X,xI,d,form)
% Compute cubic and quartic spline function
% Inputs:
% x (1x2)  : coordinate of point at which w is to be evaluated
% xI (1x2) : coord of node I
% d        : size of the support

r = sqrt( (X(1,1) - xI(1,1)) .* (X(1,1) - xI(1,1)) + (X(1,2) - xI(1,2)).*(X(1,2) - xI(1,2)) )/d ;

switch form
  case 'cubic_spline' 
     [w,dwdr,dwdrr] = cubic_spline_new(r);
  case 'quartic_spline'
     [w,dwdr] = quartic_spline(r);
  otherwise 
     error('Grr. Unknown functional form');
end

if (r ~= 0)
    drdx = (X(1,1) - xI(1,1))/(r*d*d) ;
    drdy = (X(1,2) - xI(1,2))/(r*d*d) ;
    
    
    x = X(1,1) ; 
    y = X(1,2) ; 
    xi = xI(1,1) ;
    yi = xI(1,2) ;
    drdx = (x - xi)/(d*((x - xi)^2 + (y - yi)^2)^(1/2));
    drdy = (y - yi)/(d*((x - xi)^2 + (y - yi)^2)^(1/2));
    drdxx =(y - yi)^2/(d*((x - xi)^2 + (y - yi)^2)^(3/2));
    drdyy = (x - xi)^2/(d*((x - xi)^2 + (y - yi)^2)^(3/2));
    drdxy = -((2*y - 2*yi)*(x - xi))/(2*d*((x - xi)^2 + (y - yi)^2)^(3/2));

else
    drdx = 0 ;
    drdy = 0 ;
    drdxx = 0 ; 
    drdyy = 0 ; 
    drdxy = 0 ; 
end

% dwdx = dwdr * drdx ;
% dwdy = dwdr * drdy ;

%% 
% % syms r x y xI yI d
% % r = sqrt( (x - xI) * (x - xI) + (y - yI)*(y - yI) )/d ;
% % drdx = simplify(diff(r,x));
% % drdy = simplify(diff(r,y));
% % drdxx = simplify(diff(drdx,x));
% % drdyy = simplify(diff(drdy,y));
% % drdxy = simplify(diff(drdx,y));


dwdx = dwdr * drdx ;
dwdy = dwdr * drdy ;


dwdxx = dwdrr * (drdx)^2 + drdxx * dwdr ; 
dwdyy = dwdrr * (drdy)^2 + drdyy * dwdr ; 
dwdxy = dwdrr * (drdx)*(drdy) + drdxy * dwdr ; 

