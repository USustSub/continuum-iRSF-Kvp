    if plast_it == 1 && (inc == 1 & plast_it==1) 
        [phi,B,en] = get_data ( node_c , node , di , form ) ;
        Shape2{cc_ij2}{1} = B ;  
        Shape2{cc_ij2}{2} = en ;  
        Shape2{cc_ij2}{3} = phi ; 
    else
        B = Shape2{cc_ij2}{1};
        en = Shape2{cc_ij2}{2};
        phi = Shape2{cc_ij2}{3};
    end
    cc_ij2 = cc_ij2 + 1 ;

    
if strcmp(Ex,'Duretz')    
    if node_c(1) == 0 % left edge
% set ux = 0 
% if node_c(2)<0.05
%             K(2*ij-1,en(1:2:end)) = K(2*ij-1,en(1:2:end)) + phi ;
            if (inc == 1 & plast_it==1) C_BC = [C_BC ; [en(1:2:end)'*0+2*ij-1 en(1:2:end)' phi']  ] ;  end
            F(2*ij-1,1) = F(2*ij-1,1) + 0 ;
            Residual(2*ij-1) = Residual(2*ij-1) + phi * dU(en(1:2:end)) ; 
% end
% set uy = y*incr0 
% % %             K(2*ij-0,en(2:2:end)) = K(2*ij-0,en(2:2:end)) + phi ;
% %             if (inc == 1 & plast_it==1) C_BC = [C_BC ; [en(2:2:end)'*0+2*ij-0 en(2:2:end)' phi']  ] ; end
% %             F(2*ij-0,1) = F(2*ij-0,1) + (node_c(2)-D)*inc*incr0 ;
% %             Residual(2*ij-0) = Residual(2*ij-0) + phi * dU(en(2:2:end)) ; 
% set duydx = 0 
            dphidx = B(1,1:2:end) ;
%             K(2*ij-0,en(2:2:end)) = K(2*ij-0,en(2:2:end)) + dphidx ;
            if (inc == 1 & plast_it==1) C_BC = [C_BC ; [en(2:2:end)'*0+2*ij-0 en(2:2:end)' dphidx']  ] ; end
            F(2*ij-0,1) = F(2*ij-0,1) + 0 ;
            Residual(2*ij-0) = Residual(2*ij-0) + dphidx * dU(en(2:2:end)) ; 
    elseif node_c(1) == L % right edge
% set ux = -x*incr0 
%             K(2*ij-1,en(1:2:end)) = K(2*ij-1,en(1:2:end)) + phi ;
            if (inc == 1 & plast_it==1) C_BC = [C_BC ; [en(1:2:end)'*0+2*ij-1 en(1:2:end)' phi']  ] ; end
            F(2*ij-1,1) = F(2*ij-1,1) - node_c(1)*inc*incr0 ;
            Residual(2*ij-1) = Residual(2*ij-1) + phi * dU(en(1:2:end)) ; 
% set uy = y*incr0 
% % %             K(2*ij-0,en(2:2:end)) = K(2*ij-0,en(2:2:end)) + phi ;
% %             if (inc == 1 & plast_it==1) C_BC = [C_BC ; [en(2:2:end)'*0+2*ij-0 en(2:2:end)' phi']  ] ; end
% %             F(2*ij-0,1) = F(2*ij-0,1) + (node_c(2)-D)*inc*incr0 ;
% %             Residual(2*ij-0) = Residual(2*ij-0) + phi * dU(en(2:2:end)) ; 
% % set duydx = 0 
            dphidx = B(1,1:2:end) ;
%             K(2*ij-0,en(2:2:end)) = K(2*ij-0,en(2:2:end)) + dphidx ;
            if (inc == 1 & plast_it==1) C_BC = [C_BC ; [en(2:2:end)'*0+2*ij-0 en(2:2:end)' dphidx']  ] ; end
            F(2*ij-0,1) = F(2*ij-0,1) + 0 ;
            Residual(2*ij-0) = Residual(2*ij-0) + dphidx * dU(en(2:2:end)) ; 

    end

    if node_c(2) == 0% & ~( node_c(1) ==0 | node_c(1) == L ) % bot edge
% set ux = -x*incr0 
% % %             K(2*ij-1,en(1:2:end)) = K(2*ij-1,en(1:2:end)) + phi ;
% %             if (inc == 1 & plast_it==1) C_BC = [C_BC ; [en(1:2:end)'*0+2*ij-1 en(1:2:end)' phi']  ] ; end
% %             F(2*ij-1,1) = F(2*ij-1,1) - node_c(1)*inc*incr0 ;
% %             Residual(2*ij-1) = Residual(2*ij-1) + phi * dU(en(1:2:end)) ; 

% set uy = y*incr0
%             K(2*ij-0,en(2:2:end)) = K(2*ij-0,en(2:2:end)) + phi ;
            if (inc == 1 & plast_it==1) C_BC = [C_BC ; [en(2:2:end)'*0+2*ij-0 en(2:2:end)' phi']  ] ; end
            F(2*ij-0,1) = F(2*ij-0,1) + (node_c(2)-0.7)*inc*incr0 ;
            Residual(2*ij-0) = Residual(2*ij-0) + phi * dU(en(2:2:end)) ; 
% set duxdy = 0 
            dphidy = B(2,2:2:end) ;
%             K(2*ij-1,en(1:2:end)) = K(2*ij-1,en(1:2:end)) + dphidy ;
            if (inc == 1 & plast_it==1) C_BC = [C_BC ; [en(1:2:end)'*0+2*ij-1 en(1:2:end)' dphidy']  ] ; end
            F(2*ij-1,1) = F(2*ij-1,1) + 0 ;
            Residual(2*ij-1) = Residual(2*ij-1) + dphidy * dU(en(1:2:end)) ; 
    elseif node_c(2) == D% & ~( node_c(1) ==0 | node_c(1) == L ) % top edge
% % % set ux = -x*incr0 
% % %             K(2*ij-1,en(1:2:end)) = K(2*ij-1,en(1:2:end)) + phi ;
% %             if (inc == 1 & plast_it==1) C_BC = [C_BC ; [en(1:2:end)'*0+2*ij-1 en(1:2:end)' phi']  ] ; end
% %             F(2*ij-1,1) = F(2*ij-1,1) - node_c(1)*inc*incr0 ;
% %             Residual(2*ij-1) = Residual(2*ij-1) + phi * dU(en(1:2:end)) ; 
% set uy = y*incr0 
%             K(2*ij-0,en(2:2:end)) = K(2*ij-0,en(2:2:end)) + phi ;
            if (inc == 1 & plast_it==1) C_BC = [C_BC ; [en(2:2:end)'*0+2*ij-0 en(2:2:end)' phi']  ] ; end
            F(2*ij-0,1) = F(2*ij-0,1) + (node_c(2)-0.7)*inc*incr0 ;
            Residual(2*ij-0) = Residual(2*ij-0) + phi * dU(en(2:2:end)) ;% set duxdy = 0 
% set duxdy = 0 
            dphidy = B(2,2:2:end) ;
%             K(2*ij-1,en(1:2:end)) = K(2*ij-1,en(1:2:end)) + dphidy ;
            if (inc == 1 & plast_it==1) C_BC = [C_BC ; [en(1:2:end)'*0+2*ij-1 en(1:2:end)' dphidy']  ] ; end
            F(2*ij-1,1) = F(2*ij-1,1) + 0 ;
            Residual(2*ij-1) = Residual(2*ij-1) + dphidy * dU(en(1:2:end)) ;  
    end
elseif strcmp(Ex,'Slip') || strcmp(Ex,'Herr')
     if norm(node_c(1)-0)<1e-5 % left edge
            if (inc == 1 & plast_it==1)  DD_BC = [DD_BC ; node_c ] ; end
% set uy = 0 
            K(2*ij-0,en(2:2:end)) = K(2*ij-0,en(2:2:end)) + phi ;
            F(2*ij-0,1) = F(2*ij-0,1) + 0 ;
            Residual(2*ij-0) = Residual(2*ij-0) + phi * dU(en(2:2:end)) ; 
% set duydx = 0 
            dphidx = B(1,1:2:end) ;
            K(2*ij-1,en(1:2:end)) = K(2*ij-1,en(1:2:end)) + dphidx ;
            F(2*ij-1,1) = F(2*ij-1,1) + 0 ;
            Residual(2*ij-1) = Residual(2*ij-1) + dphidx * dU(en(1:2:end)) ; 
    elseif norm(node_c(1)-L)<1e-5 % right edge
            if (inc == 1 & plast_it==1)  DD_BC = [DD_BC ; node_c ] ; end
% set uy = 0 
            K(2*ij-0,en(2:2:end)) = K(2*ij-0,en(2:2:end)) + phi ;
            F(2*ij-0,1) = F(2*ij-0,1) + 0 ;
            Residual(2*ij-0) = Residual(2*ij-0) + phi * dU(en(2:2:end)) ; 
% set duydx = 0 
            dphidx = B(1,1:2:end) ;
            K(2*ij-1,en(1:2:end)) = K(2*ij-1,en(1:2:end)) + dphidx ;
            F(2*ij-1,1) = F(2*ij-1,1) + 0 ;
            Residual(2*ij-1) = Residual(2*ij-1) + dphidx * dU(en(1:2:end)) ; 
    end

    if norm(node_c(2)-0)<1e-5% & ~( node_c(1) ==0 | node_c(1) == L ) % bot edge
            if (inc == 1 & plast_it==1)  DD_BC = [DD_BC ; node_c ] ; end
% set uy = 0 
            K(2*ij-0,en(2:2:end)) = K(2*ij-0,en(2:2:end)) + phi ;
            F(2*ij-0,1) = F(2*ij-0,1) + 0 ;
            Residual(2*ij-0) = Residual(2*ij-0) + phi * dU(en(2:2:end)) ; 
% set ux = -incr0 
            K(2*ij-1,en(1:2:end)) = K(2*ij-1,en(1:2:end)) + phi ;
            F(2*ij-1,1) = F(2*ij-1,1) - inc*incr0 ;
            Residual(2*ij-1) = Residual(2*ij-1) + phi * dU(en(1:2:end)) ; 
    elseif norm(node_c(2)-D)<1e-5% & ~( node_c(1) ==0 | node_c(1) == L ) % top edge
            if (inc == 1 & plast_it==1)  DD_BC = [DD_BC ; node_c ] ; end
% set ux = incr0 
            K(2*ij-1,en(1:2:end)) = K(2*ij-1,en(1:2:end)) + phi ;
            F(2*ij-1,1) = F(2*ij-1,1) + inc*incr0 ;
            Residual(2*ij-1) = Residual(2*ij-1) + phi * dU(en(1:2:end)) ; 
% set uy = 0 
            K(2*ij-0,en(2:2:end)) = K(2*ij-0,en(2:2:end)) + phi ;
            F(2*ij-0,1) = F(2*ij-0,1) + 0 ;
            Residual(2*ij-0) = Residual(2*ij-0) + phi * dU(en(2:2:end)) ; 
    end
end