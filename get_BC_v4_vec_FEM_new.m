    if plast_it == 1 && (inc == 1 & plast_it==1) 
% % % % %         con_cand = conn_nodes (ij,2:conn_nodes(ij,1)+1) ;
% % % % %         [phi,B,en] = get_data_Lag ( node_c  , con_cand ,  node , element , 'L' ) ;
% % % % %         Shape2{cc_ij2}{1} = B ;  
% % % % %         Shape2{cc_ij2}{2} = en ;  
% % % % %         Shape2{cc_ij2}{3} = phi ; 
        [phi,B,en] = get_data ( node_c , node , di , form ) ;
        Shape2{cc_ij2}{1} = B ;  
        Shape2{cc_ij2}{2} = en ;  
        Shape2{cc_ij2}{3} = phi ; 
%     else
%         B = Shape2{cc_ij2}{1};
%         en = Shape2{cc_ij2}{2};
%         phi = Shape2{cc_ij2}{3};
    end
    cc_ij2 = cc_ij2 + 1 ;


if strcmp(Ex,'Duretz')    
    if node_c(1) == 0 % left edge
% set ux = 0 
% if node_c(2)<0.05
            if (inc == 1 & plast_it==1)   
                K_BC_(2*ij-1,en(1:2:end)) = K_BC_(2*ij-1,en(1:2:end)) + phi ;
            end
            F(2*ij-1,1) = F(2*ij-1,1) + 0 ;
% set duydx = 0 
            if (inc == 1 & plast_it==1)  
                dphidx = B(1,1:2:end) ; 
                K_BC_(2*ij-0,en(2:2:end)) = K_BC_(2*ij-0,en(2:2:end)) + dphidx ;
            end
            F(2*ij-0,1) = F(2*ij-0,1) + 0 ;
    elseif node_c(1) == L % right edge
% set ux = -x*incr0 
            if (inc == 1 & plast_it==1) 
                K_BC_(2*ij-1,en(1:2:end)) = K_BC_(2*ij-1,en(1:2:end)) + phi ;
            end
            F(2*ij-1,1) = F(2*ij-1,1) - node_c(1)*app ;
% set uy = y*incr0 
            if (inc == 1 & plast_it==1) 
                K_BC_(2*ij-0,en(2:2:end)) = K_BC_(2*ij-0,en(2:2:end)) + phi ;
            end
            F(2*ij-0,1) = F(2*ij-0,1) + node_c(2)*app ;
    end

    if node_c(2) == 0% & ~( node_c(1) ==0 | node_c(1) == L ) % bot edge
% set uy = 0 
            if (inc == 1 & plast_it==1) 
                K_BC_(2*ij-0,en(2:2:end)) = K_BC_(2*ij-0,en(2:2:end)) + phi ;
            end
            F(2*ij-0,1) = F(2*ij-0,1) + 0 ;
% set duxdy = 0 
            if (inc == 1 & plast_it==1) 
                dphidy = B(2,2:2:end) ;
                K_BC_(2*ij-1,en(1:2:end)) = K_BC_(2*ij-1,en(1:2:end)) + dphidy ;
            end
            F(2*ij-1,1) = F(2*ij-1,1) + 0 ;
    elseif node_c(2) == D% & ~( node_c(1) ==0 | node_c(1) == L ) % top edge
% set ux = -x*incr0 
            if (inc == 1 & plast_it==1) 
                K_BC_(2*ij-1,en(1:2:end)) = K_BC_(2*ij-1,en(1:2:end)) + phi ;
            end
            F(2*ij-1,1) = F(2*ij-1,1) - node_c(1)*app ;
% set uy = y*incr0 
            if (inc == 1 & plast_it==1) 
                K_BC_(2*ij-0,en(2:2:end)) = K_BC_(2*ij-0,en(2:2:end)) + phi ;
            end
            F(2*ij-0,1) = F(2*ij-0,1) + node_c(2)*app ;
    end
elseif strcmp(Ex,'Slip') || strcmp(Ex,'Herr') || strcmp(Ex,'Junction') || strcmp(Ex,'Preuss') || strcmp(Ex,'Preuss2')
     if norm(node_c(1)-0)<1e-5 % left edge
            if (inc == 1 & plast_it==1)  
                DD_BC = [DD_BC ; node_c ] ; 
% set uy = 0 
                K_BC_(2*ij-0,en(2:2:end)) = K_BC_(2*ij-0,en(2:2:end)) + phi ;
            end
                F(2*ij-0,1) = F(2*ij-0,1) + 0 ;
% set duxdx = 0 
%             if (inc == 1 & plast_it==1)  
%                 dphidx = B(1,1:2:end) ;
%                 K_BC_(2*ij-1,en(1:2:end)) = K_BC_(2*ij-1,en(1:2:end)) + dphidx ;
%             end
            
%                 app__ = a3/a1*(dphidx*vu0(en(1:2:end))) + a5/a1*(dphidx*au0(en(1:2:end))) ;
%                 F(2*ij-1,1) = F(2*ij-1,1) + app__ ;
    elseif norm(node_c(1)-L)<1e-5 % right edge
            if (inc == 1 & plast_it==1)
                DD_BC = [DD_BC ; node_c ] ; 
% set uy = 0 
                K_BC_(2*ij-0,en(2:2:end)) = K_BC_(2*ij-0,en(2:2:end)) + phi ;
            end
                F(2*ij-0,1) = F(2*ij-0,1) + 0 ;
% set duxdx = 0 
%             if (inc == 1 & plast_it==1)
%                 dphidx = B(1,1:2:end) ;
%                 K_BC_(2*ij-1,en(1:2:end)) = K_BC_(2*ij-1,en(1:2:end)) + dphidx ;
%             end
%                 app__ = a3/a1*(dphidx*vu0(en(1:2:end))) + a5/a1*(dphidx*au0(en(1:2:end))) ;
%                 F(2*ij-1,1) = F(2*ij-1,1) + app__ ;
    end

    if norm(node_c(2)-0)<1e-5% & ~( node_c(1) ==0 | node_c(1) == L ) % bot edge
% set uy = 0 
            if (inc == 1 & plast_it==1)  
                DD_BC = [DD_BC ; node_c ] ; 
                K_BC_(2*ij-0,en(2:2:end)) = K_BC_(2*ij-0,en(2:2:end)) + phi ;
            end
            F(2*ij-0,1) = F(2*ij-0,1) + 0 ;
% set ux = -incr0 
            if (inc == 1 & plast_it==1)  
                K_BC_(2*ij-1,en(1:2:end)) = K_BC_(2*ij-1,en(1:2:end)) + phi ;
            end
            app =  1*( -uy + a3*vu0(ij*2-1) + a5*au0(ij*2-1) ) / a1 ; 
            F(2*ij-1,1) = F(2*ij-1,1) + app ;
    elseif norm(node_c(2)-D)<1e-5% & ~( node_c(1) ==0 | node_c(1) == L ) % top edge
            if (inc == 1 & plast_it==1)  
                DD_BC = [DD_BC ; node_c ] ; 
% set ux = incr0 
                K_BC_(2*ij-1,en(1:2:end)) = K_BC_(2*ij-1,en(1:2:end)) + phi ;
            end
            app =  1*( uy + a3*vu0(ij*2-1) + a5*au0(ij*2-1) ) / a1 ;
            F(2*ij-1,1) = F(2*ij-1,1) + app ;
% set uy = 0 
            if (inc == 1 & plast_it==1)  
                K_BC_(2*ij-0,en(2:2:end)) = K_BC_(2*ij-0,en(2:2:end)) + phi ;
            end
            F(2*ij-0,1) = F(2*ij-0,1) + 0 ;
    end
end



%%
% if ij == 1 
%     R1 = find(abs(node(:,1)-0)<0.00001);
%     R2 = find(abs(node(:,1)-dY)<0.00001);
%     R3 = find(abs(node(:,1)-D)<0.00001);
%     R4 = find(abs(node(:,1)-(D-dY))<0.00001);
% 
% %     figure
% %     hold on
% %     plot(node(:,1),node(:,2),'r.')
% %     axis equal
% %     plot(node(R1,1),node(R1,2),'bsq')
% %     plot(node(R2,1),node(R2,2),'gsq')
% %     plot(node(R3,1),node(R3,2),'ksq')
% %     plot(node(R4,1),node(R4,2),'csq')
% 
%     for yhh =1 : size(R1,1)
%         ixx = R1(yhh);
%         K_BC_(2*ixx-1,2*ixx-1) = 1 ;
%         ixx = R2(yhh);
%         K_BC_(2*ixx-1,2*ixx-1) = -1 ;
% 
%         ixx = R3(yhh);
%         K_BC_(2*ixx-1,2*ixx-1) = 1 ;
%         ixx = R4(yhh);
%         K_BC_(2*ixx-1,2*ixx-1) = -1 ;
%     end
% end