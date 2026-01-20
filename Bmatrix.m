function [B , N ] = Bmatrix(pt,index,node,di,form,...
                     snode,tnode,xCr,xTip,alpha)

[phi,dphidx,dphidy] = MLS_ShapeFunction(pt,index,node,di,form);

le   = length(index);
StdB = zeros(3,2*le);
StdB(1,1:2:2*le)  = dphidx ;
StdB(2,2:2:2*le)  = dphidy ;
StdB(3,1:2:2*le)  = dphidy ;
StdB(3,2:2:2*le)  = dphidx ;

% number of H and tip enriched nodes
% in index
num_split_node = size(find(snode ~= 0),2);
num_tip_node   = size(find(tnode ~= 0),2);

QT = [cos(alpha) sin(alpha); -sin(alpha) cos(alpha)];

% no enriched nodes, B is simply standard B matrix
if (num_split_node == 0) && (num_tip_node == 0)
    B = StdB;
    N = phi ; 
    % there are enriched nodes, compute enriched part, Benr
else
    Benr = []; Nenr = [];
    for m = 1 : le
        % a split enriched node
        if (snode(m) ~= 0)
            % compute Heaviside and derivatives
            dist = signed_distance(xCr,pt);
            [H,dHdx,dHdy] = heaviside_(dist);

            BI_enr = [dphidx(m)*H 0 ;
                0 dphidy(m)*H;
                dphidy(m)*H dphidx(m)*H];
            NI_enr = [ phi(m)*H ] ; 
            
            Benr   = [Benr BI_enr];
            Nenr   = [Nenr NI_enr];
            clear BI_enr ; clear NI_enr ;
            % a tip enriched node
        elseif (tnode(m) ~= 0)
            % compute branch functions
            xp = QT*(pt-xTip)';                % local coordinates
            [theta,r] = cart2pol(xp(1),xp(2)); % local polar coordinates
            [Br,dBdx,dBdy] = branch(r,theta,alpha);

            aa = dphidx(m)*Br(1) + phi(m)*dBdx(1) ;
            bb = dphidy(m)*Br(1) + phi(m)*dBdy(1);
            B1_enr = [aa 0 ; 0 bb ; bb aa];

            aa = dphidx(m)*Br(2) + phi(m)*dBdx(2) ;
            bb = dphidy(m)*Br(2) + phi(m)*dBdy(2);
            B2_enr = [aa 0 ; 0 bb ; bb aa];

            aa = dphidx(m)*Br(3) + phi(m)*dBdx(3) ;
            bb = dphidy(m)*Br(3) + phi(m)*dBdy(3);
            B3_enr = [aa 0 ; 0 bb ; bb aa];

            aa = dphidx(m)*Br(4) + phi(m)*dBdx(4) ;
            bb = dphidy(m)*Br(4) + phi(m)*dBdy(4);
            B4_enr = [aa 0 ; 0 bb ; bb aa];

            BI_enr = [B1_enr B2_enr B3_enr B4_enr];
            clear B1_enr; clear B2_enr; clear B3_enr; clear B4_enr;

            Benr = [Benr BI_enr];
            clear BI_enr ;
        end
    end   % end of loop on nodes in neighbour of Gp
    % Total B matrix
    B = [StdB Benr];
    N = [phi Nenr] ; 
    clear StdB; clear Benr; clear Nenr ; 
end      % end of check enriched nodes