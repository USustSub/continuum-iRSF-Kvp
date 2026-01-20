
clear all ;     close all ;
load Mesh_data_240.mat 
Ex = 'Duretz' ;
ifprofile = 0 ; 
pstep1 = 1 ;
pstep2 = 1 ; 

inc       = 0; % number of increment
ninc          = 20 ;     % Number of increments
gitmax        = 50;     % Max. number of global iterations
tol_plast     = 1e-10;  % Tolerance of global iterations
nnx = 240 ;
Input_ ; % input parameters and grid generation
LD = [] ; 

% initialize stresses at all nodes 
TSigma0 = 0 ; P0 = 0 ; C0 = coh0;

HH = zeros(size(node,1),(4+1+1)*4); % history parameter for L(1), B(2), R(3), T(4)
HH(:,[6 12 18 24]) = C0;  
HH_ = 0*HH ; 
HH_TT = 0*HH ; 
HH_2 = zeros(size(node,1),4*4); % history parameter
HH_2_ = 0*HH_2 ; 
HH_3 = zeros(size(node,1),4*1); % history parameter
HH_3_ = 0*HH_3+sin_phi ; 
sin_phi
HH_4 = zeros(size(node,1),4*1); % history parameter
HH_4_ = 0*HH_4+theta_ ; 

% get initial elastic solution matrix
stype = 'elastic'; dU_e = zeros(2*numnode,1); 
eta_vp
profile off 
if ifprofile == 1
profile on
end

%%
% Pre_process_ ; 
% % phi_SS = [] ; index_SS = [] ; 
% % for ii__ = 1  :size(SS,1)
% %    xx = SS(ii__ , : ) ; 
% %    [phi,B,en] = get_data ( xx , node , di , form ) ;
% %     phi_SS{ii__} = phi; 
% %     index_SS{ii__} = en ; 
% % end
%%
dU_e = 0*dU_e ;
dU = dU_e;
dU0 = 0*dU_e ; 
plastic = 0 ; % set initial flag
u = sparse(2*numnode,1);
au = u ;vu = u ;u0 = u ;au0 = u ;vu0 = u ;
time = 0 ; 
C_BC = [] ; 
while inc < ninc  % time stepping 
    inc = inc + 1 ; 
%     Update_dt ; 
    disp(['step ---> ' num2str(inc)])
    time = time + dt ; 
% Global Newton iterations : 
    for plast_it = 1 : gitmax  % until convergence
        cc_ij = 1 ;        cc_ij2 = 1 ; 
        Slip = [] ; 
        K = sparse(2*numnode,2*numnode);
        Residual = zeros(2*numnode,1) ;  maxF = 0 ;  
        F = zeros(2*numnode,1) ;

        stepF = 10000 ;
        Cell = cell(stepF ,1);
        Cell2 = cell(stepF ,1);
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

        get_dSigma_vec ;         
%         dSigma_ = [dSigma_1 dSigma_2 dSigma_3 dSigma_4];

        dP_   =-1/3*(dSigma_(:,1)+dSigma_(:,2)+dSigma_(:,4)); 

        P_   = P0_ + dP_; % get pressure
        TSigma0_ = [ HH(:,1*6-5:1*6-2);  HH(:,2*6-5:2*6-2);  HH(:,3*6-5:3*6-2);  HH(:,4*6-5:4*6-2)];
        Sigma_t_ = VE1* TSigma0_ - [P0_ P0_ P0_*0 P0_] + dSigma_; % compute total trial stress (eq. 1)

        TSigma_ = [P_ P_ P_*0 P_] + Sigma_t_; % get deviatoric components
    
        J2_ = 1/2*(TSigma_(:,1).^2+TSigma_(:,2).^2+TSigma_(:,4).^2)+TSigma_(:,3).^2;
        C     = C0; 
        sin_phi_ = [HH_3(:,1);HH_3(:,2);HH_3(:,3);HH_3(:,4)];
        F__     =  sqrt(J2_) - C.*cos_phi - P_.*sin_phi_; % trial yield function
        TSigma__ = TSigma_ ;
        P__ = P_ ; 
        J2__ = J2_ ; 
        Sigma_t__ = Sigma_t_ ;
        inx = find(F__>0) ; 
        
% plasticity occured!
        if any(F__>0)
%             disp('here')
            plastic = 1;
            dQdsigma_(:,1) = TSigma_(:,1) .* (J2_).^ (-0.5) / 2 + sin_psi / 3;
            dQdsigma_(:,2) = TSigma_(:,2) .* (J2_).^ (-0.5) / 2 + sin_psi / 3;
            dQdsigma_(:,3) = TSigma_(:,3) .* (J2_).^ (-0.5);
            dQdsigma_(:,4) = TSigma_(:,4) .* (J2_).^ (-0.5) / 2 + sin_psi / 3;

            h1 = 0;%cos_phi.*h.*sqrt(2/3).*sqrt((dQdsigma(1).^2 + 2*dQdsigma(3).^2 + dQdsigma(2).^2 + dQdsigma(4).^2));
            dg_ = (F__./(Gve_ + K__.*sin_phi_.*sin_psi + eta_vp/dt + h1 ));% increment of plastic multiplier
            
            
            normmm = 100 ; 
            dg_2 = 0*F__ ; hhhh = 0.001 ;
            while(normmm)>1e-12
        %         Res_vp ( dg_ , F__ , Gve_ , K__ , sin_phi_ , sin_psi , eta_vp , dt)
                res1 = 0 - F__ + (dg_2+hhhh).*Gve_ + (dg_2+hhhh).*K__.*sin_phi_.*sin_psi + (dg_2+hhhh)*eta_vp/dt;
                res2 = 0 - F__ + dg_2.*Gve_ + dg_2.*K__.*sin_phi_.*sin_psi + dg_2*eta_vp/dt;
        %         J = (Res(dg_+hhhh)-Res(dg_))/hhhh;
                J = (res1-res2)/hhhh; 
                ddg_ = J.\res2 ;
                dg_2 = dg_2 - ddg_ ;
                normmm = norm(ddg_(inx));
            end

            if norm(dg_2(inx)-dg_(inx))>1e-12
                error('hhh')
            end
%             dg_ = dg_2 ; 
            
            dEps_0_ = [ dEps_0_l ; dEps_0_r; dEps_0_b; dEps_0_t ] ; 
            dEps_t_2 = [dEps_0_(:,1)  dEps_0_(:,2)  dEps_0_(:,3)  dEps_0_(:,4)]  - ...
                        dg_.*[dQdsigma_(:,1) dQdsigma_(:,2) dQdsigma_(:,3)/1 dQdsigma_(:,4)]; % compute corrected strain
                        dEps_t_l2 = dEps_t_2(1:numnode,:);
                        dEps_t_r2 = dEps_t_2(numnode+1:2*numnode,:);
                        dEps_t_b2 = dEps_t_2(2*numnode+1:3*numnode,:);
                        dEps_t_t2 = dEps_t_2(3*numnode+1:4*numnode,:);

            get_dSigma_vec2 ;
        
            dP_2   =-1/3*(dSigma_2(:,1)+dSigma_2(:,2)+dSigma_2(:,4)); 
            P_2   = P0_ + dP_2; % get pressure
            Sigma_t_2 = VE1* TSigma0_ - [P0_ P0_ P0_*0 P0_] + dSigma_2; % compute total trial stress (eq. 1)

            TSigma_2 = [P_2 P_2 P_2*0 P_2] + Sigma_t_2; % get deviatoric components

            J2_2 = 1/2*(TSigma_2(:,1).^2+TSigma_2(:,2).^2+TSigma_2(:,4).^2)+TSigma_2(:,3).^2;
            F__2_     =  sqrt(J2_2) - C.*cos_phi - P_2.*sin_phi_; % trial yield function
            F_vp_   =  sqrt(J2_2) - (C.*cos_phi + P_2.*sin_phi_) - (dg_./dt).*eta_vp;  % Consistency viscoplasticity model
            
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

            TSigma__(inx,:) = TSigma_2(inx,:) ;
            Sigma_t__(inx,:) = Sigma_t_2(inx,:) ;
            P__(inx,:) = P_2(inx,:) ;
            J2__(inx) = J2_2(inx) ; 
        end
        
% vectorized stiffness calculation
        stiffness_vec ;
        Residual(1:2:2*numnode) = Residual(1:2:2*numnode) - Sigma_t__(1:numnode,1)/deltaX ; 
        Residual(1:2:2*numnode) = Residual(1:2:2*numnode) + Sigma_t__(numnode+1:2*numnode,1)/deltaX ; 
        Residual(1:2:2*numnode) = Residual(1:2:2*numnode) - Sigma_t__(2*numnode+1:3*numnode,3)/deltaY ; 
        Residual(1:2:2*numnode) = Residual(1:2:2*numnode) + Sigma_t__(3*numnode+1:4*numnode,3)/deltaY ; 

        Residual(2:2:2*numnode) = Residual(2:2:2*numnode) - Sigma_t__(1:numnode,3)/deltaX ; 
        Residual(2:2:2*numnode) = Residual(2:2:2*numnode) + Sigma_t__(numnode+1:2*numnode,3)/deltaX ; 
        Residual(2:2:2*numnode) = Residual(2:2:2*numnode) - Sigma_t__(2*numnode+1:3*numnode,2)/deltaY ; 
        Residual(2:2:2*numnode) = Residual(2:2:2*numnode) + Sigma_t__(3*numnode+1:4*numnode,2)/deltaY ; 

        
        
% unnecessary loop over nodes 
        for ij = 1 : numnode % loop over nodes 
%             ij/numnode
            node_c = node(ij,:) ; 

%             % only for nodes inside domain
%             if (node_c(1) > 0 && node_c(1) < L && node_c(2) < D && node_c(2) > 0 )
%                 x_l = node_c - [deltaX/2 0 ] ;
%                 x_r = node_c + [deltaX/2 0 ] ;    
% 
%                 x_b = node_c - [0 deltaY/2 ] ;
%                 x_t = node_c + [0 deltaY/2 ] ;    
% %%
% % knowing P0, TSigma0, C0 ---> P, TSigma, C, Dvp,
% % %                 side = 'L'; 
% % %                 get_equations_vp_vec2 ;
% % %                 side = 'R'; 
% % %                 get_equations_vp_vec2 ;
% % %                 side = 'B'; 
% % %                 get_equations_vp_vec2 ;
% % %                 side = 'T'; 
% % %                 get_equations_vp_vec2 ;
%             end
if dynamic == 1 
    K(2*ij-1,2*ij-1) = K(2*ij-1,2*ij-1) - a0*rho;
    Converged_last_step = rho*a0*u0(2*ij-1) + rho*a2*vu0(2*ij-1) + rho*a4*au0(2*ij-1) ;
    Residual(2*ij-1) = Residual(2*ij-1) + Converged_last_step - a0*rho*dU(2*ij-1);

    K(2*ij-0,2*ij-0) = K(2*ij-0,2*ij-0) - a0*rho;
    Converged_last_step = rho*a0*u0(2*ij-0) + rho*a2*vu0(2*ij-0) + rho*a4*au0(2*ij-0) ;
    Residual(2*ij-0) = Residual(2*ij-0) + Converged_last_step - a0*rho*dU(2*ij-0);
end
%             get_BCs_v4 () ;
            get_BCs_v5_2 () ;
        
        end % on nodes 
        if strcmp(Ex,'Duretz')
            K = K + sparse(C_BC(:,1),C_BC(:,2),C_BC(:,3),size(K,1),size(K,1));
        end
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

% add terms to stiffness

        nR = norm(F-Residual,1)/length(Residual); % residual norm
        disp(['iter: ' num2str(plast_it) '  norm: ' num2str(nR) ' maxF ' num2str(max(F__)) ])
% return
%         Residual = K*dU ;
%         norm(K*ddU-F)
        % check convergence 
        ddU    = K\(F-Residual) ;
        if any(isnan(ddU))
            error('ri')
        end
        dU     = dU + ddU;

    if nR < tol_plast % converged 
        break;
    end

        au = a0*(dU - u0) - (a2*vu0) - (a4*au0) ;
        vu = a1*(dU - u0) - (a3*vu0) - (a5*au0) ;

    end % on global Newton iterations 

    u0 = dU ;
    au0 = au ;
    vu0 = vu ;
    JJ = sqrt(J2__) ; JJ(JJ==0) = [] ;
    LD = [LD ; inc*incr0 mean(JJ) time];


    OutPutData
% update history parameters:
%     P0 = P; 
%     TSigma0 = TSigma ; 
%     C0        = C;
    HH = HH_;    HH_2 = HH_2_ ;   HH_3 = HH_3_ ;      HH_4 = HH_4_ ;  
    dU0 = dU ;
% pause
end % on time stepping
if ifprofile == 1
profile viewer
profsave
end

