if strcmp(side,'L')
    idx  = 1 ; 
    xx = x_l ; 
elseif strcmp(side,'R')
    idx  = 2 ; 
    xx = x_r ; 
elseif strcmp(side,'B')
    idx  = 3 ; 
    xx = x_b ; 
elseif strcmp(side,'T')
    idx  = 4 ; 
    xx = x_t ; 
end

    % get strain
    id1 = ij+idx*numnode-numnode ; 
% %     en2 = En(id1,:); en2 = en2(en2~=0); en2 = en2(en2>0);
% % 
% %     en = en2 ; 
    Eii = Eii_(id1,:); 
% %     Sigma_t = Sigma_t_(id1,:) ; 
% % 
% %     F_ = F__(id1);
% %     
% % % Plastic corrections
% %     if F_>0 % if trial stresses exceed yield 
% %         plastic = 1;
% %         Sigma_t = Sigma_t_2(id1,:) ; 
% %     end % F_ > 0
% % 
% %     maxF = max(abs(F_),maxF);
% %     kx = 0 ; ky = 0 ;
% % if idx == 1
% %     % x momentum equation - left
% %         Rx = FX(id1,1:length(en));
% %         kx = kx - Rx ;
% %         
% %     % y momentum equation - left
% %         Ry = FY(id1,1:length(en)) ; 
% %         ky = ky - Ry ;
% % end
% % 
% % 
% % if idx == 2
% %     % x momentum equation - right
% %         Rx = FX(id1,1:length(en));
% %         kx = kx + Rx ;
% % 
% %     % y momentum equation - right
% %         Ry = FY(id1,1:length(en)) ; 
% %         ky = ky + Ry ;
% % end
% % 
% % if idx == 3
% %     % x momentum equation - bot
% %         Rx = FX2(id1,1:length(en));
% %         kx = kx - Rx ;
% % 
% %     % y momentum equation - bot
% %         Ry = FY2(id1,1:length(en)) ; 
% %         ky = ky - Ry ;
% % end
% % 
% % 
% % if idx == 4
% %     % x momentum equation - top
% %         Rx = FX2(id1,1:length(en));
% %         kx = kx + Rx ;
% % 
% %     % y momentum equation - top
% %         Ry = FY2(id1,1:length(en)) ; 
% %         ky = ky + Ry ; 
% % end
% % 
% %  
%     K(2*ij-1,en) = K(2*ij-1,en) + kx ; 
%     K(2*ij  ,en) = K(2*ij  ,en) + ky ;
    
% %     [i j s] = find(kx); 
% %     Cell{c_} = [en(i)'*0+2*ij-1 ,en(j)', s'];
% %     [i j s] = find(ky); 
% %     Cell2{c_} = [en(i)'*0+2*ij ,en(j)', s'];
% %     c_ = c_ + 1 ; 

% add inertia terms 
if dynamic == 1 
    K(2*ij-1,2*ij-1) = K(2*ij-1,2*ij-1) - a0*rho;
    Converged_last_step = rho*a0*u0(2*ij-1) + rho*a2*vu0(2*ij-1) + rho*a4*au0(2*ij-1) ;
    Residual(2*ij-1) = Residual(2*ij-1) + Converged_last_step - a0*rho*dU(2*ij-1);

    K(2*ij-0,2*ij-0) = K(2*ij-0,2*ij-0) - a0*rho;
    Converged_last_step = rho*a0*u0(2*ij-0) + rho*a2*vu0(2*ij-0) + rho*a4*au0(2*ij-0) ;
    Residual(2*ij-0) = Residual(2*ij-0) + Converged_last_step - a0*rho*dU(2*ij-0);
end
    
    SS(cc_ij,:) = [ xx log10(abs(Eii))] ;
    cc_ij = cc_ij + 1 ; 
