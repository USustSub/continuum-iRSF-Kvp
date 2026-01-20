
        if any(F__>0)
% plasticity occured!
            plastic = 1;
            dQdsigma_(:,1) = TSigma_(:,1) .* (J2_).^ (-0.5) / 2 + sin_psi / 3;
            dQdsigma_(:,2) = TSigma_(:,2) .* (J2_).^ (-0.5) / 2 + sin_psi / 3;
            dQdsigma_(:,3) = TSigma_(:,3) .* (J2_).^ (-0.5);
            dQdsigma_(:,4) = TSigma_(:,4) .* (J2_).^ (-0.5) / 2 + sin_psi / 3;

            h1 = 0;%cos_phi.*h.*sqrt(2/3).*sqrt((dQdsigma(1).^2 + 2*dQdsigma(3).^2 + dQdsigma(2).^2 + dQdsigma(4).^2));
            dg_ = (F__./(Gve_ + K__.*sin_phi_.*sin_psi + eta_vp/dt + h1 ));% increment of plastic multiplier
        
            dEps_0_ = [ dEps_0_l ; dEps_0_r; dEps_0_b; dEps_0_t ] ; 
            dEps_t_2 = [dEps_0_(:,1)  dEps_0_(:,2)  dEps_0_(:,3)  dEps_0_(:,4)]  - ...
                        dg_.*[dQdsigma_(:,1) dQdsigma_(:,2) dQdsigma_(:,3)/1 dQdsigma_(:,4)]; % compute corrected strain
                        dEps_t_l2 = dEps_t_2(1:numnode,:);
                        dEps_t_r2 = dEps_t_2(numnode+1:2*numnode,:);
                        dEps_t_b2 = dEps_t_2(2*numnode+1:3*numnode,:);
                        dEps_t_t2 = dEps_t_2(3*numnode+1:4*numnode,:);

            get_dSigma_vec2 ;
%             dSigma_2 = [dSigma_1 dSigma_2 dSigma_3 dSigma_4];
        
            dP_2   =-1/3*(dSigma_2(:,1)+dSigma_2(:,2)+dSigma_2(:,4)); 
            P_2   = P0_ + dP_2; % get pressure
            Sigma_t_2 = VE1* TSigma0_ - [P0_ P0_ P0_*0 P0_] + dSigma_2; % compute total trial stress (eq. 1)

            TSigma_2 = [P_2 P_2 P_2*0 P_2] + Sigma_t_2; % get deviatoric components

            J2_2 = 1/2*(TSigma_2(:,1).^2+TSigma_2(:,2).^2+TSigma_2(:,4).^2)+TSigma_2(:,3).^2;
            F__2_     =  sqrt(J2_2) - C.*cos_phi - P_2.*sin_phi_; % trial yield function
            F_vp_   =  sqrt(J2_2) - (C.*cos_phi + P_2.*sin_phi_) - (dg_./dt).*eta_vp;  % Consistency viscoplasticity model
%             max(abs(F_vp_(F__>0)))
            dQdsigma_2(:,1) = TSigma_2(:,1) .* (J2_2).^ (-0.5) / 2 + sin_psi / 3;
            dQdsigma_2(:,2) = TSigma_2(:,2) .* (J2_2).^ (-0.5) / 2 + sin_psi / 3;
            dQdsigma_2(:,3) = TSigma_2(:,3) .* (J2_2).^ (-0.5);
            dQdsigma_2(:,4) = TSigma_2(:,4) .* (J2_2).^ (-0.5) / 2 + sin_psi / 3;
            
            dFdsigma_(:,1) = TSigma_2(:,1) .* (J2_2).^ (-0.5) / 2 + sin_phi_./ 3;
            dFdsigma_(:,2) = TSigma_2(:,2) .* (J2_2).^ (-0.5) / 2 + sin_phi_./ 3;
            dFdsigma_(:,3) = TSigma_2(:,3) .* (J2_2).^ (-0.5);
            dFdsigma_(:,4) = TSigma_2(:,4) .* (J2_2).^ (-0.5) / 2 + sin_phi_./ 3;
            
            txx = TSigma_2(:,1) ; tyy = TSigma_2(:,2) ; txy = TSigma_2(:,3) ; tzz = TSigma_2(:,4) ; 
            d2Qdsxxdsxx_ = 1 .* (1 .* J2_2) .^ (-1 ./ 2) ./ 3 - 1 .^ 2 .* txx .^ 2 .* (1 .* J2_2) .^ (-3 ./ 2) ./ 4;
            d2Qdsxxdsyy_ = -1 .* (1 .* J2_2) .^ (-1 ./ 2) ./ 6 - 1 .^ 2 .* txx .* tyy .* (1 .* J2_2) .^ (-3 ./ 2) ./ 4;
            d2Qdsxxdsxy_ = -1 .^ 2 .* txx .* txy .* (1 .* J2_2) .^ (-3 ./ 2) ./ 2;
            d2Qdsxxdszz_ = -1 .* (1 .* J2_2) .^ (-1 ./ 2) ./ 6 - 1 .^ 2 .* txx .* tzz .* (1 .* J2_2) .^ (-3 ./ 2) ./ 4;
            d2Qdsyydsxx_ = -1 .* (1 .* J2_2) .^ (-1 ./ 2) ./ 6 - 1 .^ 2 .* txx .* tyy .* (1 .* J2_2) .^ (-3 ./ 2) ./ 4;
            d2Qdsyydsyy_ = 1 .* (1 .* J2_2) .^ (-1 ./ 2) ./ 3 - 1 .^ 2 .* tyy .^ 2 .* (1 .* J2_2) .^ (-3 ./ 2) ./ 4;
            d2Qdsyydsxy_ = -1 .^ 2 .* tyy .* txy .* (1 .* J2_2) .^ (-3 ./ 2) ./ 2;
            d2Qdsyydszz_ = -1 .* (1 .* J2_2) .^ (-1 ./ 2) ./ 6 - 1 .^ 2 .* tyy .* tzz .* (1 .* J2_2) .^ (-3 ./ 2) ./ 4;
            d2Qdsxydsxx_ = -1 .^ 2 .* txx .* txy .* (1 .* J2_2) .^ (-3 ./ 2) ./ 2;
            d2Qdsxydsyy_ = -1 .^ 2 .* tyy .* txy .* (1 .* J2_2) .^ (-3 ./ 2) ./ 2;
            d2Qdsxydsxy_ = 1 .* (1 .* J2_2) .^ (-1 ./ 2) - 1 .^ 2 .* txy .^ 2 .* (1 .* J2_2) .^ (-3 ./ 2);
            d2Qdsxydszz_ = -1 .^ 2 .* txy .* tzz .* (1 .* J2_2) .^ (-3 ./ 2) ./ 2;
            d2Qdszzdsxx_ = -1 .* (1 .* J2_2) .^ (-1 ./ 2) ./ 6 - 1 .^ 2 .* txx .* tzz .* (1 .* J2_2) .^ (-3 ./ 2) ./ 4;
            d2Qdszzdsyy_ = -1 .* (1 .* J2_2) .^ (-1 ./ 2) ./ 6 - 1 .^ 2 .* tyy .* tzz .* (1 .* J2_2) .^ (-3 ./ 2) ./ 4;
            d2Qdszzdsxy_ = -1 .^ 2 .* txy .* tzz .* (1 .* J2_2) .^ (-3 ./ 2) ./ 2;
            d2Qdszzdszz_ = 1 .* (1 .* J2_2) .^ (-1 ./ 2) ./ 3 - 1 .^ 2 .* tzz .^ 2 .* (1 .* J2_2) .^ (-3 ./ 2) ./ 4;
            
            get_M ;
            
            D_ = [D11 D12 D13 D14 D21 D22 D23 D24 D31 D32 D33 D34 D41 D42 D43 D44 ];
%             Dvp = Mi * Dmat_t -  ( Mi * Dmat_t * dQdsigma * dFdsigma' * Mi * Dmat_t ) / (  eta_vp/dt + 0 + dFdsigma' * Mi * Dmat_t * dQdsigma ) ;
            FF = Multiply_ ( M_inv_ , D_) ;
%             F2 = dQdsigma * dFdsigma'
            F2 = [ dQdsigma_(:,1).*dFdsigma_(:,1) dQdsigma_(:,1).*dFdsigma_(:,2) dQdsigma_(:,1).*dFdsigma_(:,3)  dQdsigma_(:,1).*dFdsigma_(:,4) ...
                dQdsigma_(:,2).*dFdsigma_(:,1) dQdsigma_(:,2).*dFdsigma_(:,2) dQdsigma_(:,2).*dFdsigma_(:,3)  dQdsigma_(:,2).*dFdsigma_(:,4)...
                dQdsigma_(:,3).*dFdsigma_(:,1) dQdsigma_(:,3).*dFdsigma_(:,2) dQdsigma_(:,3).*dFdsigma_(:,3)  dQdsigma_(:,3).*dFdsigma_(:,4)...
                dQdsigma_(:,4).*dFdsigma_(:,1) dQdsigma_(:,4).*dFdsigma_(:,2) dQdsigma_(:,4).*dFdsigma_(:,3)  dQdsigma_(:,4).*dFdsigma_(:,4)];
            F3 = Multiply_ ( FF , F2 ) ;
            F4 = Multiply_ ( F3 , FF ) ;
            
%             dFdsigma' * Mi * Dmat_t * dQdsigma
            G = [dFdsigma_(:,1).*FF(:,1) + dFdsigma_(:,2).*FF(:,5) + dFdsigma_(:,3).*FF(:,9) + dFdsigma_(:,4).*FF(:,13) ... 
                 dFdsigma_(:,1).*FF(:,2) + dFdsigma_(:,2).*FF(:,6) + dFdsigma_(:,3).*FF(:,10) + dFdsigma_(:,4).*FF(:,14) ...
                 dFdsigma_(:,1).*FF(:,3) + dFdsigma_(:,2).*FF(:,7) + dFdsigma_(:,3).*FF(:,11) + dFdsigma_(:,4).*FF(:,15) ...
                 dFdsigma_(:,1).*FF(:,4) + dFdsigma_(:,2).*FF(:,8) + dFdsigma_(:,3).*FF(:,12) + dFdsigma_(:,4).*FF(:,16)] ;
            G2=G(:,1).*dQdsigma_(:,1)+G(:,2).*dQdsigma_(:,2)+G(:,3).*dQdsigma_(:,3)+G(:,4).*dQdsigma_(:,4);
            Dvp_ = FF - F4./(eta_vp/dt + 0 + G2);

            inx = find(F__>0) ; 
            TSigma__(inx,:) = TSigma_2(inx,:) ;
            Sigma_t__(inx,:) = Sigma_t_2(inx,:) ;
            P__(inx,:) = P_2(inx,:) ;
            J2__(inx) = J2_2(inx) ;
            dg__(inx) = dg_(inx) ; 
            dQdsigma__(inx,:) = dQdsigma_2(inx,:) ; 
        end
        