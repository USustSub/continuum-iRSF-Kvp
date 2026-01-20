    if plast_it == 1 && inc == 1 
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

    if norm(node_c(1)-0)<1e-5 % left edge
        DD = [DD ; node_c ];  
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
        DD = [DD ; node_c ];  
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
        DD = [DD ; node_c ];  
% set uy = 0 
            K(2*ij-0,en(2:2:end)) = K(2*ij-0,en(2:2:end)) + phi ;
            F(2*ij-0,1) = F(2*ij-0,1) + 0 ;
            Residual(2*ij-0) = Residual(2*ij-0) + phi * dU(en(2:2:end)) ; 
% set ux = -incr0 
            K(2*ij-1,en(1:2:end)) = K(2*ij-1,en(1:2:end)) + phi ;
            F(2*ij-1,1) = F(2*ij-1,1) - inc*incr0 ;
            Residual(2*ij-1) = Residual(2*ij-1) + phi * dU(en(1:2:end)) ; 
    elseif norm(node_c(2)-D)<1e-5% & ~( node_c(1) ==0 | node_c(1) == L ) % top edge
        DD = [DD ; node_c ];  
% set ux = incr0 
            K(2*ij-1,en(1:2:end)) = K(2*ij-1,en(1:2:end)) + phi ;
            F(2*ij-1,1) = F(2*ij-1,1) + inc*incr0 ;
            Residual(2*ij-1) = Residual(2*ij-1) + phi * dU(en(1:2:end)) ; 
% set uy = 0 
            K(2*ij-0,en(2:2:end)) = K(2*ij-0,en(2:2:end)) + phi ;
            F(2*ij-0,1) = F(2*ij-0,1) + 0 ;
            Residual(2*ij-0) = Residual(2*ij-0) + phi * dU(en(2:2:end)) ; 
    end