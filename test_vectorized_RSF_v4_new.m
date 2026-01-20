%% Clean version 
% % scp -r  ~//Downloads/Test3_2/*.m goudarzi@eejit.geo.uu.nl:/quanta1/home/goudarzi/Run3/

clear all ; close all ;
Ex = 'Herr' ;
pstep1 = 100 ; % postprocess step 1
pstep2 = 10 ; % postprocess step 2
% mesh resolution 
nnx = 180 ; 
ifprofile = 0 ; % for debugging
rewrite_mesh = 1 ;

Input_new ; % input parameters and grid generation
Create_geom ; 
Initialize_ ; % initailize parameters

fault_width = deltaY/2 ;
% fault_width = 200 ; 
% fault_width = L/1000 ;
%%
% cload all.mat
% pstep1 = 50 ; 
% pstep2 = 10 ; 
while inc < ninc  % time stepping 
    inc = inc + 1 ; 

    Update_dt ;

    eta_vp = 2.5e15;    
    eta_vp2 = eta_vp ; 

%     get current applied displacement
    app =  1*( uy + a3*vu0(topNodes(1)*2-1) + a5*au0(topNodes(1)*2-1) ) / a1 ; 
    AP_ = AP_0 +  app ;
    AP_0 = AP_ ; 
    time = time + dt ; 
    disp(['step is -- > ' num2str(inc)  ' and Time is --> ' num2str(time/60/60/24/365) ' years'])

% Global Newton iterations : 
    for plast_it = 1 : gitmax  % until convergence
        cc_ij = 1 ;        cc_ij2 = 1 ; 
        K = sparse(2*numnode,2*numnode);
        Residual = zeros(2*numnode,1) ;  maxF = 0 ;  
        F = zeros(2*numnode,1) ;

        c_ = 1 ; 
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
for iii_  = 1 : length(inx)
% dg_c(inx(iii_))
    if dg_c(inx(iii_))<0
    iii_/length(inx)
        disp('Bisection correction due to negative root') 
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
            dQdsigma_2(:,1) = TSigma_2(:,1) .* (J2_2).^ (-0.5) / 2 + sin_psi / 3;
            dQdsigma_2(:,2) = TSigma_2(:,2) .* (J2_2).^ (-0.5) / 2 + sin_psi / 3;
            dQdsigma_2(:,3) = TSigma_2(:,3) .* (J2_2).^ (-0.5);
            dQdsigma_2(:,4) = TSigma_2(:,4) .* (J2_2).^ (-0.5) / 2 + sin_psi / 3;

            dFdsigma_(:,1) = TSigma_2(:,1) .* (J2_2).^ (-0.5) / 2 + sin_phi_./ 3;
            dFdsigma_(:,2) = TSigma_2(:,2) .* (J2_2).^ (-0.5) / 2 + sin_phi_./ 3;
            dFdsigma_(:,3) = TSigma_2(:,3) .* (J2_2).^ (-0.5);
            dFdsigma_(:,4) = TSigma_2(:,4) .* (J2_2).^ (-0.5) / 2 + sin_phi_./ 3;

            ddY = fault_width ; 
            lambda_dot = dg_/dt ;  
            sxx = TSigma_2(:,1);
            syy = TSigma_2(:,2);
            sxy = TSigma_2(:,3);
            szz = TSigma_2(:,4);
            get_true_dFdsigma2;
            dFdsigma_(:,1) = f_xx; 
            dFdsigma_(:,2) = f_yy; 
            dFdsigma_(:,3) = f_xy; 
            dFdsigma_(:,4) = f_zz; 

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
            FF = Multiply_ ( M_inv_ , D_) ;
            F2 = [ dQdsigma_2(:,1).*dFdsigma_(:,1) dQdsigma_2(:,1).*dFdsigma_(:,2) dQdsigma_2(:,1).*dFdsigma_(:,3)  dQdsigma_2(:,1).*dFdsigma_(:,4) ...
                dQdsigma_2(:,2).*dFdsigma_(:,1) dQdsigma_2(:,2).*dFdsigma_(:,2) dQdsigma_2(:,2).*dFdsigma_(:,3)  dQdsigma_2(:,2).*dFdsigma_(:,4)...
                dQdsigma_2(:,3).*dFdsigma_(:,1) dQdsigma_2(:,3).*dFdsigma_(:,2) dQdsigma_2(:,3).*dFdsigma_(:,3)  dQdsigma_2(:,3).*dFdsigma_(:,4)...
                dQdsigma_2(:,4).*dFdsigma_(:,1) dQdsigma_2(:,4).*dFdsigma_(:,2) dQdsigma_2(:,4).*dFdsigma_(:,3)  dQdsigma_2(:,4).*dFdsigma_(:,4)];
            F3 = Multiply_ ( FF , F2 ) ;
            F4 = Multiply_ ( F3 , FF ) ;
            
%             dFdsigma' * Mi * Dmat_t * dQdsigma
            G = [dFdsigma_(:,1).*FF(:,1) + dFdsigma_(:,2).*FF(:,5) + dFdsigma_(:,3).*FF(:,9) + dFdsigma_(:,4).*FF(:,13) ... 
                 dFdsigma_(:,1).*FF(:,2) + dFdsigma_(:,2).*FF(:,6) + dFdsigma_(:,3).*FF(:,10) + dFdsigma_(:,4).*FF(:,14) ...
                 dFdsigma_(:,1).*FF(:,3) + dFdsigma_(:,2).*FF(:,7) + dFdsigma_(:,3).*FF(:,11) + dFdsigma_(:,4).*FF(:,15) ...
                 dFdsigma_(:,1).*FF(:,4) + dFdsigma_(:,2).*FF(:,8) + dFdsigma_(:,3).*FF(:,12) + dFdsigma_(:,4).*FF(:,16)] ;
            G2=G(:,1).*dQdsigma_2(:,1)+G(:,2).*dQdsigma_2(:,2)+G(:,3).*dQdsigma_2(:,3)+G(:,4).*dQdsigma_2(:,4);

%             F5 = eta_vp ; 
            dQdsigma_xx = dQdsigma_2(:,1) ; 
            dQdsigma_yy = dQdsigma_2(:,2) ; 
            dQdsigma_xy = dQdsigma_2(:,3) ; 
            dQdsigma_zz = dQdsigma_2(:,4) ; 
            ddY = fault_width ; 
            lambda_dot = dg_/dt ;  
            get_F5_ ; 
            F5 = -1*F5 ; 

            Dvp_ = FF - F4./(F5/dt + 0 + G2);
            
%             inx = find(F__>0) ; 
            TSigma__(inx,:) = TSigma_2(inx,:) ;
            Sigma_t__(inx,:) = Sigma_t_2(inx,:) ;
            P__(inx,:) = P_2(inx,:) ;
            J2__(inx) = J2_2(inx) ;
            dg__(inx) = dg_(inx) ; 
            dQdsigma__(inx,:) = dQdsigma_2(inx,:) ; 
        end

% vectorized stiffness calculation
        stiffness_vec3 ; 
        
        Residual(1:2:2*numnode) = Residual(1:2:2*numnode) - Sigma_t__(1:numnode,1)/deltaX ; 
        Residual(1:2:2*numnode) = Residual(1:2:2*numnode) + Sigma_t__(numnode+1:2*numnode,1)/deltaX ; 
        Residual(1:2:2*numnode) = Residual(1:2:2*numnode) - Sigma_t__(2*numnode+1:3*numnode,3)/deltaY ; 
        Residual(1:2:2*numnode) = Residual(1:2:2*numnode) + Sigma_t__(3*numnode+1:4*numnode,3)/deltaY ; 

        Residual(2:2:2*numnode) = Residual(2:2:2*numnode) - Sigma_t__(1:numnode,3)/deltaX ; 
        Residual(2:2:2*numnode) = Residual(2:2:2*numnode) + Sigma_t__(numnode+1:2*numnode,3)/deltaX ; 
        Residual(2:2:2*numnode) = Residual(2:2:2*numnode) - Sigma_t__(2*numnode+1:3*numnode,2)/deltaY ; 
        Residual(2:2:2*numnode) = Residual(2:2:2*numnode) + Sigma_t__(3*numnode+1:4*numnode,2)/deltaY ; 

%% unnecessary loop over nodes 
% tic
% for iii = 1 : 100
%     iii
if inc == 1 && plast_it == 1 
        for ij = 1 : numnode % loop over nodes 
            if inc==1
%                 ij/numnode
            end
            node_c = node(ij,:) ; 

            if dynamic == 1 
                if inc == 1 && plast_it == 1 
                    MM(2*ij-1,2*ij-1) = MM(2*ij-1,2*ij-1) - 1*rho;
                end
                    Converged_last_step = rho*a0*u0(2*ij-1) + rho*a2*vu0(2*ij-1) + rho*a4*au0(2*ij-1) ;
                    Residual(2*ij-1) = Residual(2*ij-1) + Converged_last_step - a0*rho*dU(2*ij-1);

                if inc == 1 && plast_it == 1 
                    MM(2*ij-0,2*ij-0) = MM(2*ij-0,2*ij-0) - 1*rho;
                end
                    Converged_last_step = rho*a0*u0(2*ij-0) + rho*a2*vu0(2*ij-0) + rho*a4*au0(2*ij-0) ;
                    Residual(2*ij-0) = Residual(2*ij-0) + Converged_last_step - a0*rho*dU(2*ij-0);

            end
            get_BC_v4_vec_FEM () ;
            
        end % on nodes
%         pause
end
% toc
% else
% 
if inc > 1 || plast_it > 1  
Converged_last_step = rho*a0*u0(1:2:end) + rho*a2*vu0(1:2:end) + rho*a4*au0(1:2:end) ;
Residual(1:2:end) = Residual(1:2:end) + Converged_last_step - a0*rho*dU(1:2:end);
Converged_last_step = rho*a0*u0(2:2:end) + rho*a2*vu0(2:2:end) + rho*a4*au0(2:2:end) ;
Residual(2:2:end) = Residual(2:2:end) + Converged_last_step - a0*rho*dU(2:2:end);

app =  1*( -uy + a3*vu0(botNodes*2-1) + a5*au0(botNodes*2-1) ) / a1 ; 
F(botNodes*2-1,1) = F(botNodes*2-1,1) + app ;

app =  1*( uy + a3*vu0(topNodes*2-1) + a5*au0(topNodes*2-1) ) / a1 ;
F(topNodes*2-1,1) = F(topNodes*2-1,1) + app ;
% sparse(F)

end
            Mass = MM*a0 ; 
            K = K + Mass ; 
            K = K + K_BC_ ; 
            Residual = Residual + K_BC_ * (dU-dU0) ;

% Update history 
        HH_(:,1*6-5:1*6-2) = TSigma__(1:numnode,:) ; 
        HH_(:,2*6-5:2*6-2) = TSigma__(numnode+1:2*numnode,:) ; 
        HH_(:,3*6-5:3*6-2) = TSigma__(2*numnode+1:3*numnode,:) ; 
        HH_(:,4*6-5:4*6-2) = TSigma__(3*numnode+1:4*numnode,:) ; 

        HH_(:,1*6-1) = P__(1:numnode,:) ;
        HH_(:,2*6-1) = P__(numnode+1:2*numnode,:) ;
        HH_(:,3*6-1) = P__(2*numnode+1:3*numnode,:) ;
        HH_(:,4*6-1) = P__(3*numnode+1:4*numnode,:) ; 

        HH_2_(:,1*4-3:1*4) = Eps_t_(1:numnode,:) ; 
        HH_2_(:,2*4-3:2*4) = Eps_t_(numnode+1:2*numnode,:) ; 
        HH_2_(:,3*4-3:3*4) = Eps_t_(2*numnode+1:3*numnode,:) ; 
        HH_2_(:,4*4-3:4*4) = Eps_t_(3*numnode+1:4*numnode,:) ; 


%         Residual = K*dU ;
%         norm(K*ddU-F)
        % check convergence 
        ddU    = K\(F-Residual) ;
        if any(isnan(ddU))
            error('ri')
        end
        dU     = dU + ddU;
        
        nR = norm(ddU) ; 
%         nR = norm(F-Residual,1)/length(Residual); % residual norm
        
        disp(['iter: ' num2str(plast_it) '  norm: ' num2str(nR) ' maxF ' num2str(max(F__)) ])

    if nR < tol_plast && plast_it > 1 % converged 
        break;
    end

        au = a0*(dU - u0) - (a2*vu0) - (a4*au0) ;
        vu = a1*(dU - u0) - (a3*vu0) - (a5*au0) ;
    
    end % on global Newton iterations 
    mu_gp = sin_phi_ ; 
    u0 = dU ;
    au0 = au ;
    vu0 = vu ;
    JJ = sqrt(J2__) ; JJ(JJ==0) = [] ;
    LD = [LD ; inc*incr0 mean(JJ) time dt max(Vp)];
    dlmwrite(['out/LD.log' ],[LD(:,1) LD(:,2) LD(:,3) LD(:,4) LD(:,5)],'Delimiter',' ')            
    
% get slip
    Sn = Vp*dt+S0;
    S0 = Sn ; 

% save output
    OutPutData

%% update history parameters:
%     P0 = P; 
%     TSigma0 = TSigma ; 
%     C0        = C;
    HH = HH_;    HH_2 = HH_2_ ;   HH_3 = HH_3_ ;      HH_4 = HH_4_ ;  
    dU0 = dU ;
    theta_ = theta_n ;

    if rem(inc,pstep1)==0
        save all.mat 
    end
    if rem(inc,10)==0
%         save all_tmp.mat 
    end
%     pause

end % on time stepping
if ifprofile == 1
profile viewer
profsave
end
return
%% symbolic calculation of derivation wrt dL
        syms J2 P mu Gve dL eta DT  a_ Vp V0 deltaY mu0 L_ a_ b_ theta_n theta_
        syms dQdsigma_xx dQdsigma_yy dQdsigma_zz dQdsigma_xy
        depsp = dL/DT*[dQdsigma_xx dQdsigma_yy dQdsigma_xy dQdsigma_zz];
        Ekkc_ = depsp(1) + depsp(2) + depsp(4) ;
        
        depsp_d = depsp ;
        depsp_d([1 2 4]) = depsp([1 2 4]) - 1/3*Ekkc_ ; 
        Eii_p  = sqrt(1/2*(depsp_d(1)).^2 + 1/2*(depsp_d(2)).^2 + 1/2*(depsp_d(4)).^2 + depsp_d(3).^2);
        Vp = Eii_p*deltaY ;

        param = Vp./2/V0.*exp((mu0+b_.*log(theta_n*V0/L_))/a_);
        mu = a_*asinh(param) ;
        res2 = - sqrt(J2) + P.*mu + dL.*Gve + dL*eta/DT  %+ min(0,1e14*(dg_2));
        dres2 = diff(res2,dL)

        
        
        
        S1 = 1+theta_/dt;
        S2 = (1/DT + Vp/L_);
        theta_n = S1./S2;        
        param = Vp./2/V0.*exp((mu0+b_.*log(theta_n*V0/L_))/a_);
        mu = a_*asinh(param) ;
        res2 = - sqrt(J2) + P.*mu + dL.*Gve + dL*eta/DT  %+ min(0,1e14*(dg_2));
        dres2 = diff(res2,dL)

        
%% F der with respect to lambda_dot
syms J2 P mu Gve dL eta DT  a_ Vp V0 deltaY mu0 L_ a_ b_ theta_n theta_
syms dQdsigma_xx dQdsigma_yy dQdsigma_zz dQdsigma_xy lambda_dot mu
depsp = lambda_dot*[dQdsigma_xx dQdsigma_yy dQdsigma_xy dQdsigma_zz];
Ekkc_ = depsp(1) + depsp(2) + depsp(4) ;

depsp_d = depsp ;
depsp_d([1 2 4]) = depsp([1 2 4]) - 1/3*Ekkc_ ; 
Eii_p  = sqrt(1/2*(depsp_d(1)).^2 + 1/2*(depsp_d(2)).^2 + 1/2*(depsp_d(4)).^2 + depsp_d(3).^2);
Vp = 2*Eii_p*deltaY ;
param = Vp./2/V0.*exp((mu0+b_.*log(theta_n*V0/L_))/a_);
mu = a_*asinh(param) ;

        
F = sqrt(J2)-P*mu - eta*lambda_dot ; 
diff(F,lambda_dot)

%% F der with respect to lambda
syms J2 P mu Gve dL eta DT  a_ Vp V0 deltaY mu0 L_ a_ b_ theta_n theta_
syms dQdsigma_xx dQdsigma_yy dQdsigma_zz dQdsigma_xy  mu

depsp = dL/DT*[dQdsigma_xx dQdsigma_yy dQdsigma_xy dQdsigma_zz];
Ekkc_ = depsp(1) + depsp(2) + depsp(4) ;
depsp_d = depsp ;
depsp_d([1 2 4]) = depsp([1 2 4]) - 1/3*Ekkc_ ; 
Eii_p  = sqrt(1/2*(depsp_d(1)).^2 + 1/2*(depsp_d(2)).^2 + 1/2*(depsp_d(4)).^2 + depsp_d(3).^2);
Vp = 2*Eii_p*deltaY ;

S1 = 1+theta_/DT;
S2 = (1/DT + Vp/L_);
theta_n = S1./S2;

param = Vp./2/V0.*exp((mu0+b_.*log(theta_n*V0/L_))/a_);
mu = a_*asinh(param) ;

res2 = 0 - sqrt(J2) + P*mu + dL.*Gve + dL*eta/DT ; 
dres2ddL = diff(res2,dL)


%% F der with respect to lambda_dot 2
syms J2 P mu Gve dL eta DT  a_ Vp V0 deltaY mu0 L_ a_ b_ theta_n theta_
syms dQdsigma_xx dQdsigma_yy dQdsigma_zz dQdsigma_xy lambda_dot mu
depsp = lambda_dot*[dQdsigma_xx dQdsigma_yy dQdsigma_xy dQdsigma_zz];
Ekkc_ = depsp(1) + depsp(2) + depsp(4) ;

depsp_d = depsp ;
depsp_d([1 2 4]) = depsp([1 2 4]) - 1/3*Ekkc_ ; 
Eii_p  = sqrt(1/2*(depsp_d(1)).^2 + 1/2*(depsp_d(2)).^2 + 1/2*(depsp_d(4)).^2 + depsp_d(3).^2);
Vp = 2*Eii_p*deltaY ;
S1 = 1+theta_/DT;
S2 = (1/DT + Vp/L_);
theta_n = S1./S2;

param = Vp./2/V0.*exp((mu0+b_.*log(theta_n*V0/L_))/a_);
mu = a_*asinh(param) ;

        
F = sqrt(J2) - P*mu - eta*lambda_dot ; 
diff(F,lambda_dot)


%%
 syms J2 P mu  eta_vp lambda_dot
hhhh = 1e-15;
% J2 = 10 ; P = 50 ; mu = 0.2 ; eta_vp = 1000000  ; lambda_dot = 100; 
             r1 = sqrt(J2)-P.*mu - eta_vp*(lambda_dot+hhhh*i) ; 
%              [mu,~,~,~] = get_mu_2 ( dg_ , dQdsigma_2 , dt , inx , fault_width , b_ , theta_n , 0);
             r2 = sqrt(J2)-P.*mu - eta_vp*lambda_dot ; 
             J_ = imag(r1)/hhhh;
             F5 = J_.\r2 ;
F5

diff(r2,lambda_dot)

%%
clear all
syms P_ dg_2 Gve_  eta_vp dt theta_ L_ deltaY V0 mu0 b_ a_ J2_
syms dQdsigma_xx dQdsigma_yy dQdsigma_zz dQdsigma_xy lambda_dot mu
depsp = dg_2/dt*[dQdsigma_xx dQdsigma_yy dQdsigma_xy dQdsigma_zz];
Ekkc_ = depsp(1) + depsp(2) + depsp(4) ;

depsp_d = depsp ;
depsp_d([1 2 4]) = depsp([1 2 4]) - 1/3*Ekkc_ ; 
Eii_p  = sqrt(1/2*(depsp_d(1)).^2 + 1/2*(depsp_d(2)).^2 + 1/2*(depsp_d(4)).^2 + depsp_d(3).^2);
Vp = 2*Eii_p*deltaY ;
S1 = 1+theta_/dt;
S2 = (1/dt + Vp/L_);
theta_n = S1./S2;

param = Vp./2/V0.*exp((mu0+b_.*log(theta_n*V0/L_))/a_);
mu = a_*asinh(param) ;


F= sqrt(J2_) - P_.*mu - dg_2.*Gve_ - dg_2*eta_vp/dt ==0;
ss=vpasolve(F,dg_2)
%%
close all 
h1 = figure
hold on

h2 = figure
hold on
dd_ = 0  ;
for i = [ 10:10:2100 ]
    dd_ = dd_ + 1; 
ff = dlmread(['out_5/slip_center/jj' num2str(i) '.log']);
[cc,ii] = sort(ff(:,1)); 
ff2 = dlmread(['out_5/slip_rate/jj' num2str(i) '.log']);
maxV = max(ff2(:,2))
if maxV < 1e-8
    cc = 'r-' ;
    cc2 = 'rsq' ;
elseif maxV>1e-8 && maxV<1e-2
    cc = 'b-' ;
    cc2 = 'bsq' ;
elseif maxV>1e-2 
    cc = 'k-';
    cc2 = 'ksq' ;
end

figure(h2)
plot(ff(ii,1),ff(ii,2),cc,'LineWidth',1)
ylim([0 5])

figure(h1)
subplot(1,3,1)
% hold on
plot(ff(ii,1),ff(ii,2),cc,'LineWidth',1)
ylim([0 5])
ylabel('Slip [m]')
xlabel('X [Km]')

ff_ = dlmread(['out_5/LD.log']);
subplot(1,3,2)
hold on
plot(ff_(1:2100,3)/60/60/24/365,ff_(1:2100,2),'k-','LineWidth',1)
plot(ff_(i,3)/60/60/24/365,ff_(i,2),'bsq','LineWidth',0.2)
ylabel('average shear stress [N/m2]')
xlabel('Time [Year]')


figure(h1)
subplot(1,3,3)
% plot(ff_(1:2100,3)/60/60/24/365,ff_(1:2100,2),'k-','LineWidth',1)
semilogy(ff_(i,3)/60/60/24/365,(maxV),cc2,'LineWidth',0.2)
hold on
ylim([1e-16 10])
xlim([0 90])
ylabel('maximum slip rate [m/s]')
xlabel('Time [Year]')

% eval(['print -djpeg out_5/anim/jj' num2str(dd_) , '.jpeg'])


% pause
end


%%
close all
mu0 = 0.2 ; 
V0 = 4e-9 ;  % m/s
a_ = 0.01 ;
L_ = 0.01 ;
b1_ = 0.001 ; 
b2_ = 0.017 ; 

theta_n = L_/V0*exp(40) ; 
Vp = 1e-30 ;
bb_ = [b1_:0.001:b2_];
% bb_ = b2_ ; 
Vp_ = [1e-30 1e-25 1e-20 1e-15 1e-10 1e-5 1e-2 1];

figure
hold on
for ii = 1 : length(Vp_)
Vp = Vp_(ii)
param = Vp./2/V0.*exp((mu0+bb_.*log(theta_n*V0/L_))/a_);
mu_gp = a_*asinh(param) ;

mu_gp2 = mu0 + a_*log(Vp/V0) + bb_*log(theta_n*V0/L)

plot(bb_,mu_gp,'ro-')
hold on
plot(bb_,mu_gp2,'bsq-')
end

% figure
% semilogx(Vp,mu_gp,'rsq-')
% hold on
% semilogx(Vp,mu_gp2,'bsq-')

%%
LD1 = dlmread('~/Desktop/Run_cluster/LD1.log');
LD2 = dlmread('~/Desktop/Run_cluster/LD2.log');
figure
plot(LD1(:,3)/60/60/24/365,log(LD1(:,end)/1),'rsq-','LineWidth',4)
hold on
plot(LD2(:,3)/60/60/24/365,log(LD2(:,end)/1),'bsq','MarkerSize',1)

% LD3 = dlmread('~/Downloads/Test3_2/out/LD.log');
LD3 = dlmread('~/Desktop/Run_cluster/out/LD.log');
plot(LD3(:,3)/60/60/24/365,log(LD3(:,end)/1),'gsq','MarkerSize',1)
LD3 = dlmread('~/Downloads/Test3_2/out/LD.log');
plot(LD3(:,3)/60/60/24/365,log(LD3(:,end)/1),'sq','MarkerSize',1)

%%
LD0 = dlmread('~/Downloads/Test3_3/out_0/LD.log');
LD1 = dlmread('~/Downloads/Test3_3/out_1/LD.log');
LD2 = dlmread('~/Downloads/Test3_3/out_2/LD.log');
figure
plot(LD1(:,3)/60/60/24/365,log(LD1(:,end)/1),'r.-','LineWidth',1)
hold on
plot(LD2(:,3)/60/60/24/365,log(LD2(:,end)/1),'b.','MarkerSize',1)
plot(LD0(:,3)/60/60/24/365,log(LD0(:,end)/1),'k.','MarkerSize',1)
LD3 = dlmread('~/Desktop/Run_cluster/out/LD.log');
plot(LD3(:,3)/60/60/24/365,log(LD3(:,end)/1),'g-','LineWidth',1)
LD4 = dlmread('~/Downloads/Test3_3/out/LD.log');
plot(LD4(:,3)/60/60/24/365,log(LD4(:,end)/1),'c.-','LineWidth',1)
