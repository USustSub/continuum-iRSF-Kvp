%     idx  = 1 ; 
%     xx = node_c ; 
% 
%     TSigma0 = HH(ij,idx*6-5:idx*6-2)';
%     P0 = HH(ij,idx*6-1);
%     C0 = HH(ij,idx*6-0); 
% 
%     % get strain
%     [~,B,en] = get_data ( xx , node , di , form ) ;
%     dEps_t = B*(dU(en)-dU0(en)); dEps_t = [dEps_t;0]; % add z strain
%     dEps_0 = dEps_t; % save total strain  
% 
%     % get (visco)elastic tangent 
%     [Dmat_t,Dmat,Gve,K_] = identify_tangent_v2 ( xx , mat) ;
%     dSigma = Dmat_t * dEps_t ;% total stress increment 
%     dP    =-1/3*(sum(dSigma([1 2 4]))); % pressure increment
% 
%     Sigma_t = VE1* TSigma0 - [P0;P0;0;P0] + dSigma; % compute total trial stress (eq. 1)
% 
%     P   = P0 + dP; % get pressure
%     TSigma = [P;P;0;P] + Sigma_t; % get deviatoric components
%     J2 = 1/2*(TSigma(1)^2+TSigma(2)^2+TSigma(4)^2)+TSigma(3)^2;
%     C     = C0; 
%     F_     =  sqrt(J2) - C.*cos_phi - P.*sin_phi; % trial yield function
% %     F_
%     % Plastic corrections
% %     %{
%     if F_>0 % if trial stresses exceed yield 
% %         disp('plasticity')
%         plastic = 1;
%         dQdsigma = get_dQdSigma ( TSigma , sin_psi , J2 );
%         h1 = cos_phi.*h.*sqrt(2/3).*sqrt((dQdsigma(1).^2 + 2*dQdsigma(3).^2 + dQdsigma(2).^2 + dQdsigma(4).^2));
%         dg = (F_./(Gve + K_.*sin_phi*sin_psi + eta_vp/dt + h1 ));% increment of plastic multiplier
%         dEps_t = [dEps_0(1); dEps_0(2) ; dEps_0(3); dEps_0(4)]  - ...
%             dg*[dQdsigma(1);dQdsigma(2);dQdsigma(3)/1;dQdsigma(4)]; % compute corrected strain
%         dSigma = Dmat_t * dEps_t ;% total stress increment 
%         dP    =-1/3*(sum(dSigma([1 2 4]))); % pressure increment
%         Sigma_t = VE1 .* TSigma0 - [P0;P0;0;P0] + dSigma; % compute total stresses
%         P   = P0 + dP; % get pressure
%         TSigma = [P;P;0;P] + Sigma_t; % get deviatoric components
% 
%         % Check yield function (F_vp=0)
%         J2 = 1/2*(TSigma(1).^2+TSigma(2).^2+TSigma(4).^2)+TSigma(3).^2;
%         dep=sqrt(2/3).*sqrt((dg.*dQdsigma(1)).^2+(dg.*dQdsigma(2)).^2+(dg.*dQdsigma(4)).^2+2*(dg.*dQdsigma(3)).^2); % eq. 
% 
%         C     = C0 + h.*dep; 
%         F_     =  sqrt(J2) - (C.*cos_phi + P.*sin_phi) ;                       % Backbone plastic yield function
%         F_vp   =  sqrt(J2) - (C.*cos_phi + P.*sin_phi) - (dg./dt).*eta_vp;  % Consistency viscoplasticity model
% % %         fprintf('    Backbone EP max. Fc = %2.2e \n', max(F_(:)))
% % %         fprintf('    Viscoplast. max. Fc = %2.2e \n', max(F_vp(:)))
% %         pause
% 
% 
%         % now get the tangent operator:
%         dQdsigma = get_dQdSigma ( TSigma , sin_psi , J2 )';
%         dFdsigma = get_dQdSigma ( TSigma , sin_phi , J2 )';
%         d2Qdsigma2 = get_dQ2dSigma2 ( TSigma  , J2 );
%         [Dmat_t,~,Gve,K_] = identify_tangent_v2 ( xx , mat) ;
% 
%         I = [ 1 0 0 0 ; 0 1 0 0 ; 0 0 1 0 ; 0 0 0 1] ;    
%         M = I + dg * Dmat_t * d2Qdsigma2 ;
%         Mi = inv(M) ; 
%         Dvp = Mi * Dmat_t -  ( Mi * Dmat_t * dQdsigma * dFdsigma' * Mi * Dmat_t ) / (  eta_vp/dt + h1 + dFdsigma' * Mi * Dmat_t * dQdsigma ) ; 
% %                     
%         Dmat = Dvp(1:3,1:3); % use viscoplastic tangent
% %         maxFvp = max(abs(F_vp),maxFvp);
%     end % F_ > 0
    
    
% Drichlet boundary conditions
if node_c(1) == 0 % left edge
%        plot(node_c(1),node_c(2),'bsq')
        [phi,B,en] = get_data ( node_c , node , di , form ) ;
        K(2*ij-1,en(1:2:end)) = K(2*ij-1,en(1:2:end)) + phi ;
        K(2*ij-0,en(2:2:end)) = K(2*ij-0,en(2:2:end)) + phi ;
        F(2*ij-1,1) = F(2*ij-1,1) + 0 ;
        F(2*ij-0,1) = F(2*ij-0,1) + 0 ;
        Residual(2*ij-1) = Residual(2*ij-1) + phi * dU(en(1:2:end)) ; 
        Residual(2*ij-0) = Residual(2*ij-0) + phi * dU(en(2:2:end)) ; 
elseif norm(node_c(1)-L)<1e-10 % right edge
%        plot(node_c(1),node_c(2),'bsq')
        [phi,B,en] = get_data ( node_c , node , di , form ) ;
        K(2*ij-1,en(1:2:end)) = K(2*ij-1,en(1:2:end)) + phi ;
        F(2*ij-1,1) = F(2*ij-1,1) + inc*1*incr0 ; %node_c(1)/L*0.05 ;

        K(2*ij-0,en) = K(2*ij-0,en) + Dmat(3,:)*B ;
        F(2*ij-0,1) = F(2*ij-0,1) + 0 ;
        Residual(2*ij-1) = Residual(2*ij-1) + phi * dU(en(1:2:end)) ; 
        Residual(2*ij-0) = Residual(2*ij-0) + Dmat(3,:)*B * dU(en) ; 
end

% Neuman boundary condition
if node_c(2) == 0 & ~( node_c(1) ==0 | norm(node_c(1)-L)<1e-10 ) % bot edge
%        plot(node_c(1),node_c(2),'gsq')
        [phi,B,en] = get_data ( node_c , node , di , form ) ;
        K(2*ij-1,en) = K(2*ij-1,en) + Dmat(3,:)*B; %+ Dmat(2,:)*B;
        K(2*ij-0,en) = K(2*ij-0,en) + Dmat(2,:)*B;
        F(2*ij-1,1) = F(2*ij-1,1) + 0 ;
        F(2*ij-0,1) = F(2*ij-0,1) + 0 ;
        Residual(2*ij-1) = Residual(2*ij-1) + Dmat(3,:)*B * dU(en) ; 
        Residual(2*ij-0) = Residual(2*ij-0) + Dmat(2,:)*B * dU(en) ; 
elseif norm(node_c(2)-D)<1e-10 & ~( node_c(1) ==0 | norm(node_c(1)-L)<1e-10 ) % top edge
%        plot(node_c(1),node_c(2),'gsq')
        [phi,B,en] = get_data ( node_c , node , di , form ) ;
        K(2*ij-1,en) = K(2*ij-1,en) + Dmat(3,:)*B; % ; + Dmat(2,:)*B;
        K(2*ij-0,en) = K(2*ij-0,en) + Dmat(2,:)*B;
        F(2*ij-1,1) = F(2*ij-1,1) + 0 ;
        F(2*ij-0,1) = F(2*ij-0,1) + 0 ;
        Residual(2*ij-1) = Residual(2*ij-1) + Dmat(3,:)*B * dU(en) ; 
        Residual(2*ij-0) = Residual(2*ij-0) + Dmat(2,:)*B * dU(en) ; 
end


%     HH_(ij,idx*6-5:idx*6-2) = TSigma' ;
%     HH_(ij,idx*6-1) = P ;
%     HH_(ij,idx*6-0) = C ; 
%     SS = [SS ; xx Sigma_t' J2];

