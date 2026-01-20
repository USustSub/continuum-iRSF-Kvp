    


for ij = 1 : numnode % loop over nodes 
            
    node_c = node(ij,:) ; 




if strcmp(Ex,'Slip') || strcmp(Ex,'Herr') || strcmp(Ex,'Junction') || strcmp(Ex,'Preuss') || strcmp(Ex,'Preuss2')
     if norm(node_c(1)-0)<1e-5 % left edge
         [phi,B,en] = get_data ( node_c , node , di , form ) ;
% set duxdx = 0 
                dphidx = B(1,1:2:end) ;
                app__ = a3/a1*(dphidx*vu0(en(1:2:end))) + a5/a1*(dphidx*au0(en(1:2:end))) ;
                F(2*ij-1,1) = F(2*ij-1,1) + app__ ;
    elseif norm(node_c(1)-L)<1e-5 % right edge        
        [phi,B,en] = get_data ( node_c , node , di , form ) ;
% set duxdx = 0 
                dphidx = B(1,1:2:end) ;
                app__ = a3/a1*(dphidx*vu0(en(1:2:end))) + a5/a1*(dphidx*au0(en(1:2:end))) ;
                F(2*ij-1,1) = F(2*ij-1,1) + app__ ;
     end
end

end