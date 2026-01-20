
function [globalK,globalF] = Essential_FEM...
    (globalK,globalF,dispNodes,presDofs,anode,il,ap)
if nargin == 3 
    bcwt = mean(diag(globalK)); % a measure of the average size of an element in K
    % used to keep the conditioning of the K matrix

    udofs = dispNodes.*2-1;
    vdofs = dispNodes.*2;

    globalF(udofs) = 0 ;
    globalF(vdofs) = 0 ;

    globalK(udofs,:) = 0;   % zero out the rows and columns of the K matrix
    globalK(vdofs,:) = 0;
    globalK(:,udofs) = 0;
    globalK(:,vdofs) = 0;

    globalK(udofs,udofs) = bcwt*speye(length(udofs)); % put ones*bcwt on the diagonal
    globalK(vdofs,vdofs) = bcwt*speye(length(vdofs));

else
    A = globalK ; 
    rhs = globalF ; 
    appp = setdiff(anode,presDofs) ; 
    Kff = A(appp,appp);
    Kfp = A(appp,presDofs);
    r  =   Kfp*ap ;
    if il == 1 
        globalF = rhs(appp) + r;
    else
        globalF = rhs(appp) + 0;
    end
    globalK = Kff ; 
end