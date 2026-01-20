
        dEps_t_l = [B_x_l*(dU-dU0) B_y_l*(dU-dU0) B_xy_l*(dU-dU0) 0*B_xy_l*(dU-dU0)];
        dEps_t_r = [B_x_r*(dU-dU0) B_y_r*(dU-dU0) B_xy_r*(dU-dU0) 0*B_xy_r*(dU-dU0)];
        dEps_t_t = [B_x_t*(dU-dU0) B_y_t*(dU-dU0) B_xy_t*(dU-dU0) 0*B_xy_t*(dU-dU0)];
        dEps_t_b = [B_x_b*(dU-dU0) B_y_b*(dU-dU0) B_xy_b*(dU-dU0) 0*B_xy_b*(dU-dU0)];
        dEps_0_l = dEps_t_l; % save total strain  
        dEps_0_r = dEps_t_r; % save total strain  
        dEps_0_b = dEps_t_b; % save total strain  
        dEps_0_t = dEps_t_t; % save total strain  
        dd = [dEps_t_l; dEps_t_r; dEps_t_b; dEps_t_t];
        dd(:,3) = dd(:,3)/2 ;
        Eps_t_ = [HH_2(1:numnode,1*4-3:1*4);
                  HH_2(1:numnode,2*4-3:2*4);
                  HH_2(1:numnode,3*4-3:3*4);
                  HH_2(1:numnode,4*4-3:4*4)]+dd;
        Ekkc_ = Eps_t_(:,1) + Eps_t_(:,2) + Eps_t_(:,4) ;
        Epsdc_ = Eps_t_ ;
        Epsdc_(:,[1 2 4]) = Eps_t_(:,[1 2 4]) - 1/3*Ekkc_ ; 
        Eii_  = sqrt(1/2*(Epsdc_(:,1)).^2 + 1/2*(Epsdc_(:,2)).^2 + 1/2*(Epsdc_(:,4)).^2 + Epsdc_(:,3).^2);

        P0_ = [HH(:,1*6-1);HH(:,2*6-1);HH(:,3*6-1);HH(:,4*6-1)];

        dd__=[dEps_t_l;dEps_t_r;dEps_t_b;dEps_t_t];
        dSigma_ = [D11.*dd__(:,1)+D12.*dd__(:,2)+D13.*dd__(:,3)+D14.*dd__(:,4) ...
        D21.*dd__(:,1)+D22.*dd__(:,2)+D23.*dd__(:,3)+D24.*dd__(:,4)...
        D31.*dd__(:,1)+D32.*dd__(:,2)+D33.*dd__(:,3)+D34.*dd__(:,4)...
        D41.*dd__(:,1)+D42.*dd__(:,2)+D43.*dd__(:,3)+D44.*dd__(:,4)];
%         dSigma_(:,4) = 0*dSigma_(:,4) ; 
%         dSigma_ = [dSigma_1 dSigma_2 dSigma_3 dSigma_4];

        dP_   =-1/3*(dSigma_(:,1)+dSigma_(:,2)+dSigma_(:,4)); 

        P_   = P0_ + dP_; % get pressure
        TSigma0_ = [ HH(:,1*6-5:1*6-2);  HH(:,2*6-5:2*6-2);  HH(:,3*6-5:3*6-2);  HH(:,4*6-5:4*6-2)];
        Sigma_t_ = VE1* TSigma0_ - [P0_ P0_ P0_*0 P0_] + dSigma_; % compute total trial stress (eq. 1)

        TSigma_ = [P_ P_ P_*0 P_] + Sigma_t_; % get deviatoric components

        J2_ = 1/2*(TSigma_(:,1).^2+TSigma_(:,2).^2+TSigma_(:,4).^2)+TSigma_(:,3).^2;
        C     = 1*C0; 
        F__     =  sqrt(J2_) - C.*cos_phi - P_.*sin_phi_; % trial yield function
%         pause
        TSigma__ = TSigma_ ;
        P__ = P_ ; 
        J2__ = J2_ ; 
        Sigma_t__ = Sigma_t_ ;
        dg__ = 0*F__ ; 
        dQdsigma__ = [0*F__ 0*F__ 0*F__ 0*F__] ; 

        inx = find(F__>0) ; 
        inx2 = find(C~=0);
        inx = find(C==0) ;

%         if any(F__>0) 
        if plast_it > 1
% plasticity occured!
            plastic = 1;
            dQdsigma_(:,1) = TSigma_(:,1) .* (J2_).^ (-0.5) / 2 + sin_psi / 3;
            dQdsigma_(:,2) = TSigma_(:,2) .* (J2_).^ (-0.5) / 2 + sin_psi / 3;
            dQdsigma_(:,3) = TSigma_(:,3) .* (J2_).^ (-0.5);
            dQdsigma_(:,4) = TSigma_(:,4) .* (J2_).^ (-0.5) / 2 + sin_psi / 3;

            h1 = 0;
            dg_ = (F__./(Gve_ + K__.*sin_phi_.*sin_psi + eta_vp/dt + h1 ));% increment of plastic multiplier

            normmm = 100 ; 
            if inc == 1 
                dg_2 = 0*dg_+0*1e-14 ;
            end
            
                
%%
            hhhh = 0.0000000000000000000000000000000001 ;
            ittt = 0 ; 
%             dg_2 = dg_2/100 ; 
            while (normmm)>1e-16
                ittt = ittt + 1;
                [sin_phi_0,Vp,Eii_p,theta_new] = get_mu_new ...
                    ( dg_2 , dQdsigma_ , dt , inx , fault_width , b_ , theta_ , Mat_);
                sin_phi_(inx) = sin_phi_0(inx) ; 
                theta_n = theta_new ; 
                res2 = 0 - sqrt(J2_) + P_.*sin_phi_*(1-lambda__) + dg_2.*Gve_ + dg_2.*K__.*sin_phi_.*sin_psi + dg_2*eta_vp/dt ; %+ min(0,1e14*(dg_2));
                
                if inc == 1
                    [sin_phi_hhh,~,~,~] = get_mu_new ( (dg_2+hhhh) , dQdsigma_ , dt , inx , fault_width , b_ , theta_ , Mat_);
                    [sin_phi_hhh2,~,~,~] = get_mu_new ( (dg_2+i*hhhh) , dQdsigma_ , dt , inx , fault_width , b_ , theta_ , Mat_);
                    res3 = 0 - sqrt(J2_) + P_.*sin_phi_hhh2*(1-lambda__) + (dg_2+i*hhhh).*Gve_ + (dg_2+i*hhhh).*K__.*sin_phi_hhh.*sin_psi + (dg_2+i*hhhh)*eta_vp/dt ;%+ min(0,1e14*(dg_2+hhhh));
                    J_ = imag(res3)/hhhh;
                    J = J_ ; 
                end
                
                get_der_F_2;
                if inc > 1
                    J = Der ;
                end
                ddg_ = J.\res2 ;
                dg_2 = dg_2 - ddg_ ;
                normmm_ = norm(ddg_(inx));
                if all (abs(ddg_(inx)) < 1e-19)
                    normmm = 1e-28;
                else
                    normmm = 1 ;
                end
%                 if ittt > 10
%                     warning('ittt>10')
%                 end
            end
            if ittt > 10
                warning('ittt>10')
                normmm_
%                 break ;
            end
dg_c = dg_2 ; 
resc = res2 ; 

%% Do bisection for failed nodes:
% %{
if size(find(dg_c(inx(:))<0),1) > 0
    warning (['Bisection needed for ' num2str(size(find(dg_c(inx(:))<0),1)) ' nodes'])
end
for iii_  = 1 : length(inx)
% dg_c(inx(iii_))
    if dg_c(inx(iii_))<0
%     iii_/length(inx)
%         disp('Bisection correction due to negative root') 
        dg_a = 0 ;
        dg_b = 0.5 ; 
        get_Fa_Fb
        dg_p = (dg_a + dg_b)/2;
        get_Fp ;

        err = abs(Fp);
       while err > 1e-9
%            err
       if Fa*Fp<0 
           dg_b = dg_p;
       else
           dg_a = dg_p;          
       end
        dg_p = (dg_a + dg_b)/2; 
        get_Fp ;
       err = abs(Fp);
       end
        dg_c(inx(iii_)) = dg_p;
        sin_phi_(inx(iii_)) = sin_phi_p;
        Vp(inx(iii_)) = Vp_p ; 
        theta_n(inx(iii_)) = theta_new_p;
    end
end
%}
% dg_c(dg_c<1e-28) = 1e-28 ; 
dg_2 = dg_c ; 
%% 
            dg_ = dg_2 ; 

            dEps_0_ = [ dEps_0_l ; dEps_0_r; dEps_0_b; dEps_0_t ] ; 
            dEps_t_2 = [dEps_0_(:,1)  dEps_0_(:,2)  dEps_0_(:,3)  dEps_0_(:,4)]  - ...
                        dg_.*[dQdsigma_(:,1) dQdsigma_(:,2) dQdsigma_(:,3)/1 dQdsigma_(:,4)]; % compute corrected strain
                        dEps_t_l2 = dEps_t_2(1:numnode,:);
                        dEps_t_r2 = dEps_t_2(numnode+1:2*numnode,:);
                        dEps_t_b2 = dEps_t_2(2*numnode+1:3*numnode,:);
                        dEps_t_t2 = dEps_t_2(3*numnode+1:4*numnode,:);

            dd__=[dEps_t_l2;dEps_t_r2;dEps_t_b2;dEps_t_t2];
            dSigma_2 = [D11.*dd__(:,1)+D12.*dd__(:,2)+D13.*dd__(:,3)+D14.*dd__(:,4) ...
                D21.*dd__(:,1)+D22.*dd__(:,2)+D23.*dd__(:,3)+D24.*dd__(:,4)...
                D31.*dd__(:,1)+D32.*dd__(:,2)+D33.*dd__(:,3)+D34.*dd__(:,4)...
                D41.*dd__(:,1)+D42.*dd__(:,2)+D43.*dd__(:,3)+D44.*dd__(:,4)];
%             dSigma_2 = [dSigma_1 dSigma_2 dSigma_3 dSigma_4];
%             dSigma_2(:,4) = 0*dSigma_2(:,4) ; 

            dP_2   =-1/3*(dSigma_2(:,1)+dSigma_2(:,2)+dSigma_2(:,4)); 
            P_2   = P0_ + dP_2; % get pressure
            Sigma_t_2 = VE1* TSigma0_ - [P0_ P0_ P0_*0 P0_] + dSigma_2; % compute total trial stress (eq. 1)

            TSigma_2 = [P_2 P_2 P_2*0 P_2] + Sigma_t_2; % get deviatoric components

            J2_2 = 1/2*(TSigma_2(:,1).^2+TSigma_2(:,2).^2+TSigma_2(:,4).^2)+TSigma_2(:,3).^2;

            F__2_     =  sqrt(J2_2) - C.*cos_phi - P_2.*sin_phi_; % trial yield function
            F_vp_   =  sqrt(J2_2) - (C.*cos_phi + P_2.*sin_phi_) - (dg_./dt).*eta_vp;  % Consistency viscoplasticity model
%             disp(['Max F_vp_ --> ' num2str(max(abs(F_vp_(inx)))) ])


            
%             inx = find(F__>0) ; 
            TSigma__(inx,:) = TSigma_2(inx,:) ;
            Sigma_t__(inx,:) = Sigma_t_2(inx,:) ;
            P__(inx,:) = P_2(inx,:) ;
            J2__(inx) = J2_2(inx) ;
            dg__(inx) = dg_(inx) ; 
        end

% vectorized stiffness calculation
%         stiffness_vec3 ; 
        
        Residual(1:2:2*numnode) = Residual(1:2:2*numnode) - Sigma_t__(1:numnode,1)./dXX ; 
        Residual(1:2:2*numnode) = Residual(1:2:2*numnode) + Sigma_t__(numnode+1:2*numnode,1)./dXX ; 
        Residual(1:2:2*numnode) = Residual(1:2:2*numnode) - Sigma_t__(2*numnode+1:3*numnode,3)./dYY ; 
        Residual(1:2:2*numnode) = Residual(1:2:2*numnode) + Sigma_t__(3*numnode+1:4*numnode,3)./dYY ; 

        Residual(2:2:2*numnode) = Residual(2:2:2*numnode) - Sigma_t__(1:numnode,3)./dXX ; 
        Residual(2:2:2*numnode) = Residual(2:2:2*numnode) + Sigma_t__(numnode+1:2*numnode,3)./dXX ; 
        Residual(2:2:2*numnode) = Residual(2:2:2*numnode) - Sigma_t__(2*numnode+1:3*numnode,2)./dYY ; 
        Residual(2:2:2*numnode) = Residual(2:2:2*numnode) + Sigma_t__(3*numnode+1:4*numnode,2)./dYY ; 

%% unnecessary loop over nodes 
% tic
% for iii = 1 : 100
%     iii

Converged_last_step = rho*a0*u0(1:2:end) + rho*a2*vu0(1:2:end) + rho*a4*au0(1:2:end) ;
Residual(1:2:end) = Residual(1:2:end) + Converged_last_step - a0*rho*dU(1:2:end);
Converged_last_step = rho*a0*u0(2:2:end) + rho*a2*vu0(2:2:end) + rho*a4*au0(2:2:end) ;
Residual(2:2:end) = Residual(2:2:end) + Converged_last_step - a0*rho*dU(2:2:end);

app =  1*( -uy + a3*vu0(botNodes*2-1) + a5*au0(botNodes*2-1) ) / a1 ; 
F(botNodes*2-1,1) = F(botNodes*2-1,1) + app ;

app =  1*( uy + a3*vu0(topNodes*2-1) + a5*au0(topNodes*2-1) ) / a1 ;
F(topNodes*2-1,1) = F(topNodes*2-1,1) + app ;
% sparse(F)

% addV ;


            Residual = Residual + K_BC_ * (dU-dU0)-F ;
            
            
            
            