% *************************************************************************
%                 TWO DIMENSIONAL ELEMENT FREE GALERKIN CODE
%                              S.Sh.Ghorashi
%                                  2009
% *************************************************************************
function [node] = square_node_array_ir(xA,xB,xC,xD,densx,densy)
% Compute irregular coordinates of node in the domain
% densx(:,1) indicate the number of division in x direction and densx(:,2)
% indicate portion of that case
% densy(:,1) indicate the number of division in y direction and densy(:,2)
% indicate portion of that case


nnx = 1 + sum(densx(:,1)) ;
nny = 1 + sum(densy(:,1)) ;
node = zeros(nnx*nny,2) ;
nnode = 0 ;
for j = 1 : nny
    if ( j-1 <= densy(1,1) )
        x1(1)=xA(1)+(j-1)*(xD(1)-xA(1))*densy(1,2)/densy(1,1);
        x1(2)=xA(2)+(j-1)*(xD(2)-xA(2))*densy(1,2)/densy(1,1);
        x2(1)=xB(1)+(j-1)*(xC(1)-xB(1))*densy(1,2)/densy(1,1);
        x2(2)=xB(2)+(j-1)*(xC(2)-xB(2))*densy(1,2)/densy(1,1);
    else
        for jj = 2 : size(densy,1) % row number in densy matrix
            if ( j-1-sum(densy(1:jj-1,1)) <= densy(jj,1) )
                x1(1)=xA(1)+(j-1-sum(densy(1:jj-1,1)))*(xD(1)-xA(1))...
                    *densy(jj,2)/densy(jj,1)+(xD(1)-xA(1))*sum(densy(1:jj-1,2));
                x1(2)=xA(2)+(j-1-sum(densy(1:jj-1,1)))*(xD(2)-xA(2))...
                    *densy(jj,2)/densy(jj,1)+(xD(2)-xA(2))*sum(densy(1:jj-1,2));
                x2(1)=xB(1)+(j-1-sum(densy(1:jj-1,1)))*(xC(1)-xB(1))...
                    *densy(jj,2)/densy(jj,1)+(xC(1)-xB(1))*sum(densy(1:jj-1,2));
                x2(2)=xB(2)+(j-1-sum(densy(1:jj-1,1)))*(xC(2)-xB(2))...
                    *densy(jj,2)/densy(jj,1)+(xC(2)-xB(2))*sum(densy(1:jj-1,2));
                break
            end
        end
    end
    for i = 1 : nnx
        nnode = nnode + 1 ;
        if ( i-1 <= densx(1,1) )
            node(nnode,1)=x1(1)+(i-1)*(x2(1)-x1(1))*densx(1,2)/densx(1,1);
            node(nnode,2)=x1(2)+(i-1)*(x2(2)-x1(2))*densx(1,2)/densx(1,1);
        else
            for ii = 2 : size(densx,1) % row number in densx matrix
                if ( i-1-sum(densx(1:ii-1,1)) <= densx(ii,1) )
                    node(nnode,1)=x1(1)+(i-1-sum(densx(1:ii-1,1)))*...
                        (x2(1)-x1(1))*densx(ii,2)/densx(ii,1)+...
                        (x2(1)-x1(1))*sum(densx(1:ii-1,2));
                    node(nnode,2)=x1(2)+(i-1-sum(densx(1:ii-1,1)))*...
                        (x2(2)-x1(2))*densx(ii,2)/densx(ii,1)+...
                        (x2(2)-x1(2))*sum(densx(1:ii-1,2));
                    break
                end
            end
        end
    end
end
clear nnode;