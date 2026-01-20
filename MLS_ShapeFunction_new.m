% *************************************************************************
%                 TWO DIMENSIONAL ELEMENT FREE GALERKIN CODE
%                            Nguyen Vinh Phu
%                        LTDS, ENISE, Juillet 2006
% *************************************************************************

function [phi,dphidx,dphidy,dphidxx,dphidyy,dphidxy] = MLS_ShapeFunction_new(pt,index,node,di,form)

% Compute the MLS shape function at point pt for all nodes within the
% support of this point pt.
% Basis used is linear basis pT = [1 x y]

% --------------------------------------
%      compute the moment matrix A
% --------------------------------------

A    = zeros(3,3) ;
dAdx = zeros(3,3) ;
dAdy = zeros(3,3) ;
dAdxx = zeros(3,3) ;
dAdyy = zeros(3,3) ;
dAdxy = zeros(3,3) ;
w = zeros(1,length(index)); 
dwdx = zeros(1,length(index)); 
dwdy = zeros(1,length(index)); 
dwdxx = zeros(1,length(index)); 
dwdyy = zeros(1,length(index)); 
dwdxy = zeros(1,length(index)); 
for m = 1 : size(index,2)
    xi = [node(index(m),1) node(index(m),2)] ;
    [wi,dwidx,dwidy,dwidxx,dwidyy,dwidxy] = circle_spline_new(pt,xi,di(index(m)),form);
%     pTp = [1 xi(1,1) xi(1,2)]'*[1 xi(1,1) xi(1,2)] ;
    x = xi(1,1) ; y = xi(1,2) ;
    pTp = [1 x y ; x x*x x*y ; y x*y y*y] ;
    A    = A    + wi*pTp ;
    dAdx = dAdx + dwidx*pTp ;
    dAdy = dAdy + dwidy*pTp ;
    dAdxx = dAdxx + dwidxx*pTp ;
    dAdyy = dAdyy + dwidyy*pTp ;
    dAdxy = dAdxy + dwidxy*pTp ;
    % store weight function and its derivative at node I for later use
    w(m)    = wi ;
    dwdx(m) = dwidx ;
    dwdy(m) = dwidy ;
    dwdxx(m) = dwidxx ;
    dwdyy(m) = dwidyy ;
    dwdxy(m) = dwidxy ;
end

% clear wi; clear dwidx; clear dwidy ; clear xi;

p  = [1; pt(1,1); pt(1,2)];

% --------------------------------------
%         compute  matrix c(x)
% --------------------------------------

% A(x)c(x)   = p(x)
% A(x)c,k(x) = b,k(x)
% Backward substitutions, two times for c(x), two times for c,k(x) k
% =1,2

% Using LU factorization for A
% [L,U,PERM] = lu(A) ;
% c = zeros(3,3) ; 
% for i = 1 : 3
%     if i == 1         % backward substitution for c(x)
%         C = PERM*p;
%     elseif i == 2     % backward substitution for c,x(x)
%         C = PERM*([0 1 0]' - dAdx*c(1:3,1));
%     elseif i == 3     % backward substitution for c,y(x)
%         C = PERM*([0 0 1]' - dAdy*c(1:3,1));
%     end
% 
%     D1 = C(1);
%     D2 = C(2) - L(2,1)*D1;
%     D3 = C(3) - L(3,1)*D1 - L(3,2)*D2 ;
% 
%     c(3,i) = D3/U(3,3) ;
%     c(2,i) = (D2 - U(2,3)*c(3,i))/(U(2,2));
%     c(1,i) = (D1 - U(1,2)*c(2,i) - U(1,3)*c(3,i))/(U(1,1));
% end

c = zeros(3,3) ; 
cc = inv(A)*p;
c(:,1) = inv(A)*p ;
px = [0 1 0] ;
py = [0 0 1] ;
pxx = [ 0 0 0 ]; 
pxy = [ 0 0 0 ]; 
pyy = [ 0 0 0 ]; 

c(:,2) = inv(A)*(-dAdx*cc+px');
c(:,3) = inv(A)*(-dAdy*cc+py');

cx = c(:,2) ; 
cy = c(:,3) ; 

cxx = inv(A)*(pxx'-dAdx*cx-dAdx*cx-dAdxx*cc);
cyy = inv(A)*(pyy'-dAdy*cy-dAdy*cy-dAdyy*cc);
cxy = inv(A)*(pxy'-dAdy*cx-dAdx*cy-dAdxy*cc);


phi = zeros(1,length(index)) ; 
dphidx = zeros(1,length(index)) ; 
dphidy = zeros(1,length(index)) ; 
dphidxx = zeros(1,length(index)) ; 
dphidyy = zeros(1,length(index)) ; 
dphidxy = zeros(1,length(index)) ; 

for m = 1 : size(index,2)
    xi = [node(index(m),1) node(index(m),2)] ;
    piT = [1 xi(1,1) xi(1,2)]';
    phi(m) = c(:,1)'* piT*w(m) ;
    dphidx(m) = c(:,2)'*piT*w(m) + c(:,1)'*piT*dwdx(m) ;
    dphidy(m) = c(:,3)'*piT*w(m) + c(:,1)'*piT*dwdy(m);

    
    dphidx (m) = cx' * piT * w(m) + c(:,1)'*piT*dwdx(m) ; 
    dphidy (m) = cy' * piT * w(m) + c(:,1)'*piT*dwdy(m) ; 

    dphidxx (m) = cxx' * piT * w(m) + cx'*piT*dwdx(m) + cx'*piT*dwdx(m) + c(:,1)'*piT*dwdxx(m) ; 
    dphidyy (m) = cyy' * piT * w(m) + cy'*piT*dwdy(m) + cy'*piT*dwdy(m) + c(:,1)'*piT*dwdyy(m) ; 
    dphidxy (m) = cxy' * piT * w(m) + cx'*piT*dwdy(m) + cy'*piT*dwdx(m) + c(:,1)'*piT*dwdxy(m) ; 
end





% % % %% symbolic 
% % % syms x y xI yI d 
% % % syms a0 a1 a2 
% % % assume(x>0)
% % % assume(y>0)
% % % assume(a0>0)
% % % assume(a1>0)
% % % assume(a2>0)
% % % 
% % % a = [a0 ; a1 ; a2 ] ;
% % % p = [ 1 ; x ; y ] ; 
% % % uh = p'* a ; 
% % % 
% % % A = w * p * p' ; 
% % % 
% % % 
% % % 
% % % r = sqrt( (x - xI) * (x - xI) + (y - yI)*(y - yI) )/d ;
% % % 
% % %      [w,dwdr] = cubic_spline(r);
% % % %%
% % % if (r ~= 0)
% % %     drdx = (x(1,1) - xI(1,1))/(r*d*d) ;
% % %     drdy = (x(1,2) - xI(1,2))/(r*d*d) ;
% % % else
% % %     drdx = 0 ;
% % %     drdy = 0 ;
% % % end
% % % 
% % % dwdx = dwdr * drdx ;
% % % dwdy = dwdr * drdy ;




