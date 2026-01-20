
clear all ;     close all ;
% load Mesh_data_440.mat 
Ex = 'Herr' ;
ifprofile = 0 ; 
pstep1 = 1000 ;
pstep2 = 10 ; 

inc       = 0; % number of increment
ninc          = 20000000 ;     % Number of increments
gitmax        = 2500;     % Max. number of global iterations
tol_plast     = 1e-12;  % Tolerance of global iterations
nnx = 80 ;
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
sin_phi_ = [HH_3_(:,1);HH_3_(:,2);HH_3_(:,3);HH_3_(:,4)];

HH_4 = zeros(size(node,1),4*1); % history parameter
HH_4_ = 0*HH_4+theta_ ; 

% get initial elastic solution matrix
stype = 'elastic'; dU_e = zeros(2*numnode,1); 
% 566
profile off 
if ifprofile == 1
profile on
end
    theta_n = theta_ ; 
    Vp = 0 ; 
%%
Pre_process_ ; 
% load RSF/Mesh_data_80.mat 
sin_phi_ = 0*sin_phi_+0.000000;  
mu_gp = sin_phi_ ; 
% eta_vp = 0*2.5e15;       
if strcmp(Ex,'Herr')
    b_ = zeros(length(SS0),1)+b2_;
    b_((SS0(:,1)<LL1 | SS0(:,1)>LL2))=b1_;
end
%%
dU_e = 0*dU_e ;
dU = dU_e;
dU0 = 0*dU_e ; 
plastic = 0 ; % set initial flag
u = sparse(2*numnode,1);
au = u ;vu = u ;u0 = u ;au0 = u ;vu0 = u ;
time = 0 ; 
C_BC = [] ; 
K_BC_ = sparse(2*numnode,2*numnode)  ; 
topNodes = find(node(:,2)==D); 
botNodes = find(node(:,2)==D); 
AP_0 = 0 ; S0 = 0 ; 
uy = 2e-9 ;
%%
 MM = sparse(2*numnode,2*numnode);
 load all.mat
% gitmax        = 125;
% tol_pltast     = 1e-17;
% pstep1 = 1; 
% pstep2 = 1; 
while inc < ninc  % time stepping 
    inc = inc + 1 ; 

    Update_dt ;  

%     if inc > 250

%          eta_vp = 2.5e9;
% %         eta_vp2 = 2.5e16;
        eta_vp2 = eta_vp ; 
%     end
    
%     get current applied displacement
app =  1*( uy + a3*vu0(topNodes(1)*2-1) + a5*au0(topNodes(1)*2-1) ) / a1 ; 
AP_ = AP_0 +  app ;
% vu0(topNodes(1)*2-1)
AP_0 = AP_ ; 
% pause
%     disp(['step ---> ' num2str(inc)])
    time = time + dt ; 
    disp(['step is -- > ' num2str(inc)  ' and Time is --> ' num2str(time/60/60/24/365) ' years'])

% Global Newton iterations : 
    for plast_it = 1 : gitmax  % until convergence
        cc_ij = 1 ;        cc_ij2 = 1 ; 
        Slip = [] ; 
        K = sparse(2*numnode,2*numnode);
        Residual = zeros(2*numnode,1) ;  maxF = 0 ;  
        F = zeros(2*numnode,1) ;
       
                
%         if inc > 1
%             S1 = 1+theta_/dt;
%             S2 = (1/dt + Vp/L_);
%             theta_n = S1./S2;
%         end
        
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
        F__     =  sqrt(J2_) - C.*cos_phi - P_.*sin_phi_; % trial yield function
        TSigma__ = TSigma_ ;
        P__ = P_ ; 
        J2__ = J2_ ; 
        Sigma_t__ = Sigma_t_ ;
        dg__ = 0*F__ ; 
        dQdsigma__ = [0*F__ 0*F__ 0*F__ 0*F__] ; 
%         max(sqrt(J2_(inx)))
%         if plast_it == 1 
%              F_vp_ = F__ ; 
%         end
        inx = find(F__>0) ; 
        inx2 = find(C~=0);
        inx = find(C==0) ; 
%         F__(inx) = 0.001;
%         if plast_it == 1 
%             F__ = 0*F__; 
%             inx = find(F__>0) ; 
%         end
%         Vp = 0*F__ ;
%         if any(F__>0) 
        if plast_it > 1
%         if nnz(J2_) > 0 
% plasticity occured!
            plastic = 1;
            dQdsigma_(:,1) = TSigma_(:,1) .* (J2_).^ (-0.5) / 2 + sin_psi / 3;
            dQdsigma_(:,2) = TSigma_(:,2) .* (J2_).^ (-0.5) / 2 + sin_psi / 3;
            dQdsigma_(:,3) = TSigma_(:,3) .* (J2_).^ (-0.5);
            dQdsigma_(:,4) = TSigma_(:,4) .* (J2_).^ (-0.5) / 2 + sin_psi / 3;

            h1 = 0;%cos_phi.*h.*sqrt(2/3).*sqrt((dQdsigma(1).^2 + 2*dQdsigma(3).^2 + dQdsigma(2).^2 + dQdsigma(4).^2));
            dg_ = (F__./(Gve_ + K__.*sin_phi_.*sin_psi + eta_vp/dt + h1 ));% increment of plastic multiplier
% %{
            normmm = 100 ; 
            if inc == 1 
            dg_2 = 0*dg_ ;
            end
            hhhh = 0.00000000000000000000000001 ;
            ittt = 0 ; 
            while (normmm)>1e-16 || ittt < 1
                ittt = ittt + 1; 
%                 Update_mu  ; 
                [sin_phi_0,Vp,Eii_p,theta_new] = get_mu ( dg_2 , dQdsigma_ , dt , inx , deltaY/2 , b_ , theta_n );
                sin_phi_(inx) = sin_phi_0(inx) ; 
                [sin_phi_hhh,~,~,~] = get_mu ( (dg_2+hhhh) , dQdsigma_ , dt , inx , deltaY/2 , b_ , theta_n );
                [sin_phi_hhh2,~,~,~] = get_mu ( (dg_2+i*hhhh) , dQdsigma_ , dt , inx , deltaY/2 , b_ , theta_n );
                theta_n = theta_new ; 
        %         Res_vp ( dg_ , F__ , Gve_ , K__ , sin_phi_ , sin_psi , eta_vp , dt)
                res1 = 0 - sqrt(J2_) + P_.*sin_phi_hhh + (dg_2+hhhh).*Gve_ + (dg_2+hhhh).*K__.*sin_phi_hhh.*sin_psi + (dg_2+hhhh)*eta_vp/dt ;%+ min(0,1e14*(dg_2+hhhh));
                res2 = 0 - sqrt(J2_) + P_.*sin_phi_ + dg_2.*Gve_ + dg_2.*K__.*sin_phi_.*sin_psi + dg_2*eta_vp/dt ; %+ min(0,1e14*(dg_2));
                res3 = 0 - sqrt(J2_) + P_.*sin_phi_hhh2 + (dg_2+i*hhhh).*Gve_ + (dg_2+i*hhhh).*K__.*sin_phi_hhh.*sin_psi + (dg_2+i*hhhh)*eta_vp/dt ;%+ min(0,1e14*(dg_2+hhhh));
    %         J = (Res(dg_+hhhh)-Res(dg_))/hhhh;
                J = (res1-res2)/hhhh;
                J_ = imag(res3)/hhhh;
%                 [J(J>0) J_(J>0)];
                J = J_ ; 
%                 pause
% dres2 = Gve_ + eta_vp/DT + (P*a_*deltaY*exp((mu0 + b*log((V0*theta_n)/L_))/a_)*((dQdsigma_xx/(3*DT) + dQdsigma_yy/(3*DT) - (2*dQdsigma_zz)/(3*DT))*((dL*dQdsigma_xx)/(3*DT) + (dL*dQdsigma_yy)/(3*DT) - (2*dL*dQdsigma_zz)/(3*DT)) + (dQdsigma_xx/(3*DT) - (2*dQdsigma_yy)/(3*DT) + dQdsigma_zz/(3*DT))*((dL*dQdsigma_xx)/(3*DT) - (2*dL*dQdsigma_yy)/(3*DT) + (dL*dQdsigma_zz)/(3*DT)) + (dQdsigma_yy/(3*DT) - (2*dQdsigma_xx)/(3*DT) + dQdsigma_zz/(3*DT))*((dL*dQdsigma_yy)/(3*DT) - (2*dL*dQdsigma_xx)/(3*DT) + (dL*dQdsigma_zz)/(3*DT)) + (2*dL*dQdsigma_xy^2)/DT^2))/(4*V0*((deltaY^2*exp((2*mu0 + 2*b*log((V0*theta_n)/L_))/a_)*(((dL*dQdsigma_xx)/(3*DT) + (dL*dQdsigma_yy)/(3*DT) - (2*dL*dQdsigma_zz)/(3*DT))^2/2 + ((dL*dQdsigma_xx)/(3*DT) - (2*dL*dQdsigma_yy)/(3*DT) + (dL*dQdsigma_zz)/(3*DT))^2/2 + ((dL*dQdsigma_yy)/(3*DT) - (2*dL*dQdsigma_xx)/(3*DT) + (dL*dQdsigma_zz)/(3*DT))^2/2 + (dL^2*dQdsigma_xy^2)/DT^2))/(4*V0^2) + 1)^(1/2)*(((dL*dQdsigma_xx)/(3*DT) + (dL*dQdsigma_yy)/(3*DT) - (2*dL*dQdsigma_zz)/(3*DT))^2/2 + ((dL*dQdsigma_xx)/(3*DT) - (2*dL*dQdsigma_yy)/(3*DT) + (dL*dQdsigma_zz)/(3*DT))^2/2 + ((dL*dQdsigma_yy)/(3*DT) - (2*dL*dQdsigma_xx)/(3*DT) + (dL*dQdsigma_zz)/(3*DT))^2/2 + (dL^2*dQdsigma_xy^2)/DT^2)^(1/2));
                
                
                get_der_F;
                if inc > 1
                    J = Der ;
                end
%                 pause
                
                
%                 J(inx(1))
%                 pause
                ddg_ = J.\res2 ;
                dg_2 = dg_2 - ddg_ ;
                normmm = norm(ddg_(inx));
%                 ddg_(inx)
% max(ddg_(inx))
%                 min(abs(res2(inx)))
%                 normmm 
%                 max(dg_2(inx))

            if ittt > 10
                warning('ittt>10')
                normmm
%                 break ;
            end
            end
%             ittt 
%             pause
%             sin_phi_(inx) = sin_phi_t(inx) ; 
%             length(inx)
%             max(dg_2(inx))
%             figure
% % plot(SS0(:,1),sin_phi_(:,1),'rsq')
% plot(SS0(sin_phi_>0,1),sin_phi_(sin_phi_>0,1),'rsq')

%}
            %%
%{
%             for ii = 1 %: length(inx)
                dg_correct = 0*dg_ ;
                idx_n = inx ; 
                for dpar = 0 : 1e-18 : 0.01
%                     dpar
                    dg_2(idx_n) = 0*dg_(idx_n) + dpar;
    %                 Update_mu;
                    [sin_phi_,mu_gp] = get_mu ( dg_2 , dQdsigma_ , dt , inx , deltaY , b_ , theta_n);
                    res2 = 0 - sqrt(J2_) + P_.*sin_phi_ + dg_2.*Gve_ + dg_2*eta_vp/dt ; %+ min(0,1e14*(dg_2));
                    res2(idx_n)
                    idx=find(res2(inx)>0);
                    dg_correct(inx(idx)) = dg_2(inx(idx));
                    idx_n = setdiff(inx,inx(idx));
                    if isempty(idx_n)
                        break
                    end
                end
%             end
            dg_ =dg_correct; 
            max(dg_(inx))
            
%%
        F2 =[ ] ;
        for dpar =  7.6881e-14 : 1e-5 : 0.3
%             dpar
            dg_2 = 0*dg_ + dpar;
            [sin_phi_,mu_gp] = get_mu ( dg_2 , dQdsigma_ , dt , inx , deltaY , b_ , theta_n);
            res2 = 0 - sqrt(J2_) + P_.*sin_phi_ + dg_2.*Gve_ + dg_2*eta_vp/dt ; %+ min(0,1e14*(dg_2));
            F2 = [F2 ; dg_2(inx(1)) res2(inx(1)) ];
            res2(inx(1))
        end
        
        figure
        hold on
        plot(F(:,1),F(:,2),'r-')
        plot(F2(:,1),F2(:,2),'r-')
        
        %% 
                hhhh = 0.00000000000000000000001 ;
                dg_2 = 0*dg_2; 
                [sin_phi_,~] = get_mu ( dg_2 , dQdsigma_ , dt , inx , deltaY , b_ , theta_n);
                [sin_phi_hhh,~] = get_mu ( (dg_2+hhhh) , dQdsigma_ , dt , inx , deltaY , b_ , theta_n);
                res1 = 0 - sqrt(J2_) + P_.*sin_phi_hhh + (dg_2+hhhh).*Gve_ + (dg_2+hhhh).*K__.*sin_phi_hhh.*sin_psi + (dg_2+hhhh)*eta_vp/dt ;%+ min(0,1e14*(dg_2+hhhh));
                res2 = 0 - sqrt(J2_) + P_.*sin_phi_ + dg_2.*Gve_ + dg_2.*K__.*sin_phi_.*sin_psi + dg_2*eta_vp/dt ; %+ min(0,1e14*(dg_2));
        %         J = (Res(dg_+hhhh)-Res(dg_))/hhhh;
                J = (res1-res2)/hhhh; 
                J(inx(1))
                res2(inx(1))
%                 ddg_ = J.\res2 ;
%                 dg_2 = dg_2 - ddg_ ;
%                 normmm = norm(ddg_(inx));
Gve = Gve_(inx(1)) ; 
DT = dt ; 
P = P_(inx(1));  
b = b_(inx(1)) ; 
dQdsigma_xx = dQdsigma_(inx(1),1);
dQdsigma_yy = dQdsigma_(inx(1),2);
dQdsigma_xy = dQdsigma_(inx(1),3);
dQdsigma_zz = dQdsigma_(inx(1),4); 
dL = dg_2(inx(1));
eta = eta_vp ;
J2 = J2_(inx(1)); 
res2_ = Gve*dL - sqrt(J2) + P*a_*asinh((deltaY*exp((mu0 + b*log((V0*theta_n)/L_))/a_)*(((dL*dQdsigma_xx)/(3*DT) + (dL*dQdsigma_yy)/(3*DT) - (2*dL*dQdsigma_zz)/(3*DT))^2/2 + ((dL*dQdsigma_xx)/(3*DT) - (2*dL*dQdsigma_yy)/(3*DT) + (dL*dQdsigma_zz)/(3*DT))^2/2 + ((dL*dQdsigma_yy)/(3*DT) - (2*dL*dQdsigma_xx)/(3*DT) + (dL*dQdsigma_zz)/(3*DT))^2/2 + (dL^2*dQdsigma_xy^2)/DT^2)^(1/2))/(2*V0)) + (dL*eta)/DT;
dL = dL+hhhh ; 
res1_= Gve*dL - sqrt(J2) + P*a_*asinh((deltaY*exp((mu0 + b*log((V0*theta_n)/L_))/a_)*(((dL*dQdsigma_xx)/(3*DT) + (dL*dQdsigma_yy)/(3*DT) - (2*dL*dQdsigma_zz)/(3*DT))^2/2 + ((dL*dQdsigma_xx)/(3*DT) - (2*dL*dQdsigma_yy)/(3*DT) + (dL*dQdsigma_zz)/(3*DT))^2/2 + ((dL*dQdsigma_yy)/(3*DT) - (2*dL*dQdsigma_xx)/(3*DT) + (dL*dQdsigma_zz)/(3*DT))^2/2 + (dL^2*dQdsigma_xy^2)/DT^2)^(1/2))/(2*V0)) + (dL*eta)/DT;
JJ = (res1_-res2_)/hhhh


res2_/res2(inx(1))
dres2 = Gve + eta/DT + (P*a_*deltaY*exp((mu0 + b*log((V0*theta_n)/L_))/a_)*((dQdsigma_xx/(3*DT) + dQdsigma_yy/(3*DT) - (2*dQdsigma_zz)/(3*DT))*((dL*dQdsigma_xx)/(3*DT) + (dL*dQdsigma_yy)/(3*DT) - (2*dL*dQdsigma_zz)/(3*DT)) + (dQdsigma_xx/(3*DT) - (2*dQdsigma_yy)/(3*DT) + dQdsigma_zz/(3*DT))*((dL*dQdsigma_xx)/(3*DT) - (2*dL*dQdsigma_yy)/(3*DT) + (dL*dQdsigma_zz)/(3*DT)) + (dQdsigma_yy/(3*DT) - (2*dQdsigma_xx)/(3*DT) + dQdsigma_zz/(3*DT))*((dL*dQdsigma_yy)/(3*DT) - (2*dL*dQdsigma_xx)/(3*DT) + (dL*dQdsigma_zz)/(3*DT)) + (2*dL*dQdsigma_xy^2)/DT^2))/(4*V0*((deltaY^2*exp((2*mu0 + 2*b*log((V0*theta_n)/L_))/a_)*(((dL*dQdsigma_xx)/(3*DT) + (dL*dQdsigma_yy)/(3*DT) - (2*dL*dQdsigma_zz)/(3*DT))^2/2 + ((dL*dQdsigma_xx)/(3*DT) - (2*dL*dQdsigma_yy)/(3*DT) + (dL*dQdsigma_zz)/(3*DT))^2/2 + ((dL*dQdsigma_yy)/(3*DT) - (2*dL*dQdsigma_xx)/(3*DT) + (dL*dQdsigma_zz)/(3*DT))^2/2 + (dL^2*dQdsigma_xy^2)/DT^2))/(4*V0^2) + 1)^(1/2)*(((dL*dQdsigma_xx)/(3*DT) + (dL*dQdsigma_yy)/(3*DT) - (2*dL*dQdsigma_zz)/(3*DT))^2/2 + ((dL*dQdsigma_xx)/(3*DT) - (2*dL*dQdsigma_yy)/(3*DT) + (dL*dQdsigma_zz)/(3*DT))^2/2 + ((dL*dQdsigma_yy)/(3*DT) - (2*dL*dQdsigma_xx)/(3*DT) + (dL*dQdsigma_zz)/(3*DT))^2/2 + (dL^2*dQdsigma_xy^2)/DT^2)^(1/2));


        %% 
%}
%% 
%             dg_ (inx2) = 0 ; 
%             F__(inx2) = -100000000 ; 
            
            dg_ = dg_2 ; 
%             if norm(dg_2(inx)-dg_(inx))>1e-12
%                 error('hhh')
%             end

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
%             disp(['Max F_vp_ --> ' num2str(max(abs(F_vp_(inx)))) ])
            dQdsigma_2(:,1) = TSigma_2(:,1) .* (J2_2).^ (-0.5) / 2 + sin_psi / 3;
            dQdsigma_2(:,2) = TSigma_2(:,2) .* (J2_2).^ (-0.5) / 2 + sin_psi / 3;
            dQdsigma_2(:,3) = TSigma_2(:,3) .* (J2_2).^ (-0.5);
            dQdsigma_2(:,4) = TSigma_2(:,4) .* (J2_2).^ (-0.5) / 2 + sin_psi / 3;
            
            
            
% get dFdsigma
            dFdsigma_(:,1) = TSigma_2(:,1) .* (J2_2).^ (-0.5) / 2 + sin_phi_./ 3;
            dFdsigma_(:,2) = TSigma_2(:,2) .* (J2_2).^ (-0.5) / 2 + sin_phi_./ 3;
            dFdsigma_(:,3) = TSigma_2(:,3) .* (J2_2).^ (-0.5);
            dFdsigma_(:,4) = TSigma_2(:,4) .* (J2_2).^ (-0.5) / 2 + sin_phi_./ 3;            
            ddY = deltaY/2 ; 
            lambda_dot = dg_/dt ;  
            sxx = TSigma_2(:,1);
            syy = TSigma_2(:,2);
            sxy = TSigma_2(:,3);
            szz = TSigma_2(:,4);
            get_true_dFdsigma;
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
%             Dvp = Mi * Dmat_t -  ( Mi * Dmat_t * dQdsigma * dFdsigma' * Mi * Dmat_t ) / (  eta_vp/dt + 0 + dFdsigma' * Mi * Dmat_t * dQdsigma ) ;
            FF = Multiply_ ( M_inv_ , D_) ;
%             F2 = dQdsigma * dFdsigma'
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
            ddY = deltaY/2 ; 
            lambda_dot = dg_/dt ;  
%             F5 =- eta_vp2 -...
%                  (P_2.*a_.*ddY.*exp((mu0 + b_.*log((V0.*theta_n)/L_))/a_).*(2.*lambda_dot.*dQdsigma_xy.^2 + ...
%                  (dQdsigma_xx./3 + dQdsigma_yy./3 - (2.*dQdsigma_zz)/3).*((dQdsigma_xx.*lambda_dot)/3 + ...
%                  (dQdsigma_yy.*lambda_dot)/3 - (2.*dQdsigma_zz.*lambda_dot)/3) + (dQdsigma_xx./3 - (2.*dQdsigma_yy)/3 + dQdsigma_zz./3).*((dQdsigma_xx.*lambda_dot)/3 - (2.*dQdsigma_yy.*lambda_dot)/3 + (dQdsigma_zz.*lambda_dot)/3) + (dQdsigma_yy./3 - (2.*dQdsigma_xx)/3 + dQdsigma_zz./3).*((dQdsigma_yy.*lambda_dot)/3 - (2.*dQdsigma_xx.*lambda_dot)/3 + (dQdsigma_zz.*lambda_dot)/3)))./...
%                  (2.*V0.*((ddY.^2.*exp((2.*mu0 + 2.*b_.*log((V0.*theta_n)./L_))./a_).*(((dQdsigma_xx.*lambda_dot)/3 + (dQdsigma_yy.*lambda_dot)/3 - (2.*dQdsigma_zz.*lambda_dot)/3).^2./2 + ((dQdsigma_xx.*lambda_dot)/3 - (2.*dQdsigma_yy.*lambda_dot)/3 + (dQdsigma_zz.*lambda_dot)/3).^2./2 + ((dQdsigma_yy.*lambda_dot)/3 - (2.*dQdsigma_xx.*lambda_dot)/3 + (dQdsigma_zz.*lambda_dot)/3).^2./2 + dQdsigma_xy.^2.*lambda_dot.^2))./V0.^2 + 1).^(1./2).*(((dQdsigma_xx.*lambda_dot)/3 + (dQdsigma_yy.*lambda_dot)/3 - (2.*dQdsigma_zz.*lambda_dot)/3).^2./2 + ((dQdsigma_xx.*lambda_dot)/3 - (2.*dQdsigma_yy.*lambda_dot)/3 + (dQdsigma_zz.*lambda_dot)/3).^2./2 + ((dQdsigma_yy.*lambda_dot)/3 - (2.*dQdsigma_xx.*lambda_dot)/3 + (dQdsigma_zz.*lambda_dot)/3).^2./2 + dQdsigma_xy.^2.*lambda_dot.^2).^(1./2));
            F5 = - eta_vp2 - (P_2.*a_.*ddY.*exp((mu0 + b_.*log((V0.*theta_n)./L_))./a_).*(2.*lambda_dot.*dQdsigma_xy.^2 + (dQdsigma_xx./3 + dQdsigma_yy./3 - (2.*dQdsigma_zz)./3).*((dQdsigma_xx.*lambda_dot)./3 + (dQdsigma_yy.*lambda_dot)./3 - (2.*dQdsigma_zz.*lambda_dot)./3) + (dQdsigma_xx./3 - (2.*dQdsigma_yy)./3 + dQdsigma_zz./3).*((dQdsigma_xx.*lambda_dot)./3 - (2.*dQdsigma_yy.*lambda_dot)./3 + (dQdsigma_zz.*lambda_dot)./3) + (dQdsigma_yy./3 - (2.*dQdsigma_xx)./3 + dQdsigma_zz./3).*((dQdsigma_yy.*lambda_dot)./3 - (2.*dQdsigma_xx.*lambda_dot)./3 + (dQdsigma_zz.*lambda_dot)./3)))./(2.*V0.*((ddY.^2.*exp((2.*mu0 + 2.*b_.*log((V0.*theta_n)./L_))./a_).*(((dQdsigma_xx.*lambda_dot)./3 + (dQdsigma_yy.*lambda_dot)./3 - (2.*dQdsigma_zz.*lambda_dot)./3).^2./2 + ((dQdsigma_xx.*lambda_dot)./3 - (2.*dQdsigma_yy.*lambda_dot)./3 + (dQdsigma_zz.*lambda_dot)./3).^2./2 + ((dQdsigma_yy.*lambda_dot)./3 - (2.*dQdsigma_xx.*lambda_dot)./3 + (dQdsigma_zz.*lambda_dot)./3).^2./2 + dQdsigma_xy.^2.*lambda_dot.^2))./V0.^2 + 1).^(1./2).*(((dQdsigma_xx.*lambda_dot)./3 + (dQdsigma_yy.*lambda_dot)./3 - (2.*dQdsigma_zz.*lambda_dot)./3).^2./2 + ((dQdsigma_xx.*lambda_dot)./3 - (2.*dQdsigma_yy.*lambda_dot)./3 + (dQdsigma_zz.*lambda_dot)./3).^2./2 + ((dQdsigma_yy.*lambda_dot)./3 - (2.*dQdsigma_xx.*lambda_dot)./3 + (dQdsigma_zz.*lambda_dot)./3).^2./2 + dQdsigma_xy.^2.*lambda_dot.^2).^(1./2));



%             F5 = -eta_vp ... 
%                 - (P_2.*a_.*((ddY.*exp((mu0 + b_.*log((V0.*(theta_./dt + 1))./(L_.*(1./dt + ...
%                 (2.*ddY.*(((dQdsigma_xx.*lambda_dot)./3 + (dQdsigma_yy.*lambda_dot)./3 -...
%                 (2.*dQdsigma_zz.*lambda_dot)./3).^2./2 + ((dQdsigma_xx.*lambda_dot)./3 - ...
%                 (2.*dQdsigma_yy.*lambda_dot)./3 + (dQdsigma_zz.*lambda_dot)./3).^2./2 + ...
%                 ((dQdsigma_yy.*lambda_dot)./3 - (2.*dQdsigma_xx.*lambda_dot)./3 + ...
%                 (dQdsigma_zz.*lambda_dot)./3).^2./2 + dQdsigma_xy.^2.*lambda_dot.^2).^(1./2))./L_))))./a_).*(2.*lambda_dot.*dQdsigma_xy.^2 +...
%                 (dQdsigma_xx./3 + dQdsigma_yy./3 - (2.*dQdsigma_zz)./3).*((dQdsigma_xx.*lambda_dot)./3 + (dQdsigma_yy.*lambda_dot)./3 -...
%                 (2.*dQdsigma_zz.*lambda_dot)./3) + (dQdsigma_xx./3 - (2.*dQdsigma_yy)./3 + dQdsigma_zz./3).*((dQdsigma_xx.*lambda_dot)./3 -...
%                 (2.*dQdsigma_yy.*lambda_dot)./3 + (dQdsigma_zz.*lambda_dot)./3) + (dQdsigma_yy./3 - (2.*dQdsigma_xx)./3 + dQdsigma_zz./3).*((dQdsigma_yy.*lambda_dot)./3 -...
%                 (2.*dQdsigma_xx.*lambda_dot)./3 + (dQdsigma_zz.*lambda_dot)./3)))./(2.*V0.*(((dQdsigma_xx.*lambda_dot)./3 + (dQdsigma_yy.*lambda_dot)./3 - ...
%                 (2.*dQdsigma_zz.*lambda_dot)./3).^2./2 + ((dQdsigma_xx.*lambda_dot)./3 - (2.*dQdsigma_yy.*lambda_dot)./3 + (dQdsigma_zz.*lambda_dot)./3).^2./2 + ...
%                 ((dQdsigma_yy.*lambda_dot)./3 - (2.*dQdsigma_xx.*lambda_dot)./3 +...
%                 (dQdsigma_zz.*lambda_dot)./3).^2./2 + dQdsigma_xy.^2.*lambda_dot.^2).^(1./2)) - ...
%                 (b_.*ddY.^2.*exp((mu0 + b_.*log((V0.*(theta_./dt + 1))./(L_.*(1./dt + (2.*ddY.*(((dQdsigma_xx.*lambda_dot)./3 +...
%                 (dQdsigma_yy.*lambda_dot)./3 - (2.*dQdsigma_zz.*lambda_dot)./3).^2./2 + ...
%                 ((dQdsigma_xx.*lambda_dot)./3 - (2.*dQdsigma_yy.*lambda_dot)./3 + (dQdsigma_zz.*lambda_dot)./3).^2./2 + ...
%                 ((dQdsigma_yy.*lambda_dot)./3 - (2.*dQdsigma_xx.*lambda_dot)./3 + (dQdsigma_zz.*lambda_dot)./3).^2./2 +...
%                 dQdsigma_xy.^2.*lambda_dot.^2).^(1./2))./L_))))./a_).*(2.*lambda_dot.*dQdsigma_xy.^2 + (dQdsigma_xx./3 + dQdsigma_yy./3 - (2.*dQdsigma_zz)./3).*((dQdsigma_xx.*lambda_dot)./3 + (dQdsigma_yy.*lambda_dot)./3 - (2.*dQdsigma_zz.*lambda_dot)./3) + (dQdsigma_xx./3 - (2.*dQdsigma_yy)./3 + dQdsigma_zz./3).*((dQdsigma_xx.*lambda_dot)./3 - (2.*dQdsigma_yy.*lambda_dot)./3 + (dQdsigma_zz.*lambda_dot)./3) + (dQdsigma_yy./3 - (2.*dQdsigma_xx)./3 + dQdsigma_zz./3).*((dQdsigma_yy.*lambda_dot)./3 - (2.*dQdsigma_xx.*lambda_dot)./3 + (dQdsigma_zz.*lambda_dot)./3)))./(L_.*V0.*a_.*(1./dt + (2.*ddY.*(((dQdsigma_xx.*lambda_dot)./3 + (dQdsigma_yy.*lambda_dot)./3 - (2.*dQdsigma_zz.*lambda_dot)./3).^2./2 + ((dQdsigma_xx.*lambda_dot)./3 - (2.*dQdsigma_yy.*lambda_dot)./3 + (dQdsigma_zz.*lambda_dot)./3).^2./2 + ((dQdsigma_yy.*lambda_dot)./3 - (2.*dQdsigma_xx.*lambda_dot)./3 + (dQdsigma_zz.*lambda_dot)./3).^2./2 + dQdsigma_xy.^2.*lambda_dot.^2).^(1./2))./L_))))./((ddY.^2.*exp((2.*mu0 + 2.*b_.*log((V0.*(theta_./dt + 1))./(L_.*(1./dt + (2.*ddY.*(((dQdsigma_xx.*lambda_dot)./3 + (dQdsigma_yy.*lambda_dot)./3 - (2.*dQdsigma_zz.*lambda_dot)./3).^2./2 + ((dQdsigma_xx.*lambda_dot)./3 - (2.*dQdsigma_yy.*lambda_dot)./3 + (dQdsigma_zz.*lambda_dot)./3).^2./2 + ((dQdsigma_yy.*lambda_dot)./3 - (2.*dQdsigma_xx.*lambda_dot)./3 + (dQdsigma_zz.*lambda_dot)./3).^2./2 + dQdsigma_xy.^2.*lambda_dot.^2).^(1./2))./L_))))./a_).*(((dQdsigma_xx.*lambda_dot)./3 + (dQdsigma_yy.*lambda_dot)./3 - (2.*dQdsigma_zz.*lambda_dot)./3).^2./2 + ((dQdsigma_xx.*lambda_dot)./3 - (2.*dQdsigma_yy.*lambda_dot)./3 + (dQdsigma_zz.*lambda_dot)./3).^2./2 + ((dQdsigma_yy.*lambda_dot)./3 - (2.*dQdsigma_xx.*lambda_dot)./3 + (dQdsigma_zz.*lambda_dot)./3).^2./2 + dQdsigma_xy.^2.*lambda_dot.^2))./V0.^2 + 1).^(1./2);
            
            
%             hhhh = 0.00001; 
%              [mu,~,~,~] = get_mu_2 ( dg_ , dQdsigma_2 , dt , inx , deltaY/2 , b_ , theta_n , hhhh);
%              mu = 0.1;  
%              r1 = sqrt(J2_2)-P_2.*mu - eta_vp*(lambda_dot+hhhh*i) ; 
% %              [mu,~,~,~] = get_mu_2 ( dg_ , dQdsigma_2 , dt , inx , deltaY/2 , b_ , theta_n , 0);
%              r2 = sqrt(J2_2)-P_2.*mu - eta_vp*lambda_dot ; 
%              J_ = imag(r1)/hhhh;
%              F5 = J_.\r2 ;

%              max(abs(F5(inx)))
                        
             
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
%         if inc <= 450
%          stiffness_vec2 ; 
%         else
         stiffness_vec3 ; 
%         end
        
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
            Mass = a0*MM  ; 
            K = K + Mass ; 
            K = K + K_BC_ ; 
            Residual = Residual + K_BC_ * (dU-dU0) ;
            rightNodes = find(node(:,1) == L) ;
            leftNodes = find(node(:,1) == 0) ;
            topNodes = find(node(:,2) == D) ;
            botNodes = find(node(:,2) == 0) ;
%             if plast_it == 1 
%                 Residual(2*rightNodes-1) = 0 ;
%                 Residual(2*rightNodes) = 0 ;
%             end
%             mean(dU(2*topNodes-1))
%             mean(dU(2*botNodes-1))
%             mean(dU(2*leftNodes))
%             mean(dU(2*rightNodes))
%             condest(K)
%         if strcmp(Ex,'Duretz')
%             K = K + sparse(C_BC(:,1),C_BC(:,2),C_BC(:,3),size(K,1),size(K,1));
%         end
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

% pause
%         Residual = K*dU ;
%         norm(K*ddU-F)
        % check convergence 
        ddU    = K\(F-Residual) ;
        if any(isnan(ddU))
            error('ri')
        end
        dU     = dU + ddU;
        
        nR = norm(ddU) ; 
        nR2 = norm(F-Residual,1)/length(Residual); % residual norm
        
        disp(['iter: ' num2str(plast_it) '  norm: ' num2str(nR) '  norm2: ' num2str(nR2)  ' maxF ' num2str(max(F__)) ])

  

        au = a0*(dU - u0) - (a2*vu0) - (a4*au0) ;
        vu = a1*(dU - u0) - (a3*vu0) - (a5*au0) ;

    if nR < tol_plast % converged 
        break;
    end
%         pause
        
    end % on global Newton iterations 
    mu_gp = sin_phi_ ; 
    u0 = dU ;
    au0 = au ;
    vu0 = vu ;
    JJ = sqrt(J2__) ; JJ(JJ==0) = [] ;
    LD = [LD ; inc*incr0 mean(JJ) time dt];
    dlmwrite(['out/LD.log' ],[LD(:,1) LD(:,2) LD(:,3) ],'Delimiter',' ')            
    
    Sn = Vp*dt+S0;
    S0 = Sn ; 
%     figure
% % plot(SS0(:,1),sin_phi_(:,1),'rsq')
% plot(SS0(sin_phi_>0,1),sin_phi_(sin_phi_>0,1),'rsq')
% pause
    OutPutData

%% Update friction coefficient
if plastic==1% && rem(inc,10)==0
%         dt__ = 10*dt ;
        S1 = 1+theta_/dt;
        S2 = (1/dt + Vp/L_);
        theta_n = S1./S2;
        theta_n(theta_n<1000) = 1000; 
        theta_ = theta_n ; 
end
max(theta_)
min(theta_)

%         depsp = [dg__.*dQdsigma__(:,1) dg__.*dQdsigma__(:,2) dg__.*dQdsigma__(:,3) dg__.*dQdsigma__(:,4)]/dt;
%         
%         Ekkc_ = depsp(:,1) + depsp(:,2) + depsp(:,4) ;
%         depsp_d = depsp ;
%         depsp_d(:,[1 2 4]) = depsp(:,[1 2 4]) - 1/3*Ekkc_ ; 
%         Eii_p  = sqrt(1/2*(depsp_d(:,1)).^2 + 1/2*(depsp_d(:,2)).^2 + 1/2*(depsp_d(:,4)).^2 + depsp_d(:,3).^2);
%         Vp = Eii_p*deltaY ;
                
%         sin_phi_ = get_mu ( dg__ , dQdsigma__ , dt , inx , deltaY , b_ , theta_n );

%{
if strcmp(Ex,'Herr')
        depsp = [dg__.*dQdsigma__(:,1) dg__.*dQdsigma__(:,2) dg__.*dQdsigma__(:,3) dg__.*dQdsigma__(:,4)]/dt;
        
        Ekkc_ = depsp(:,1) + depsp(:,2) + depsp(:,4) ;
        depsp_d = depsp ;
        depsp_d(:,[1 2 4]) = depsp(:,[1 2 4]) - 1/3*Ekkc_ ; 
        Eii_p  = sqrt(1/2*(depsp_d(:,1)).^2 + 1/2*(depsp_d(:,2)).^2 + 1/2*(depsp_d(:,4)).^2 + depsp_d(:,3).^2);
        Vp = 2*Eii_p*deltaY ;
        
%         Vp(Vp== 0) = 1e-14;
        
% get slip
        u_t = [] ; u_b = [] ; 
        for ij = 1 : size(phi_DD_t,1)
            u_t(ij,:) = [phi_DD_t(ij,:)*dU(en_DD_t(ij,1:2:size(en_DD_t,2))) phi_DD_t(ij,:)*dU(en_DD_t(ij,1:2:size(en_DD_t,2)))];
            u_b(ij,:) = [phi_DD_b(ij,:)*dU(en_DD_b(ij,2:2:size(en_DD_b,2))) phi_DD_b(ij,:)*dU(en_DD_b(ij,2:2:size(en_DD_b,2)))];
        end
        v_t = [] ; v_b = [] ; 
        for ij = 1 : size(phi_DD_t,1)
            v_t(ij,:) = [phi_DD_t(ij,:)*vu(en_DD_t(ij,1:2:size(en_DD_t,2))) phi_DD_t(ij,:)*vu(en_DD_t(ij,1:2:size(en_DD_t,2)))];
            v_b(ij,:) = [phi_DD_b(ij,:)*vu(en_DD_b(ij,2:2:size(en_DD_b,2))) phi_DD_b(ij,:)*vu(en_DD_b(ij,2:2:size(en_DD_b,2)))];
        end
        slip = u_t(:,1) - u_b(:,1) ; 
        slip_rate = (v_t(:,1)) - (v_b(:,1)) ; 
%         Vp = 0*Vp ; 
%         Vp(DD(:,3)) = abs(slip_rate) ; 
        
% update theta 
        S1 = 1+theta_/dt;
        S2 = (1/dt + Vp/L_);
        theta_n = S1./S2;
% steady state 
%         theta_n = L_./Vp ; 
% update friction coefficient
        param = Vp./2/V0.*exp((mu0+b_.*log(theta_n*V0/L_))/a_);
        mu_gp = a_*asinh(param) ;
%         if all(F__<0)
%             mu_gp = 0*mu_gp + 0.08;
%         end
        mu_gp(mu_gp==0) = sin_phi ;
        sin_phi_ = mu_gp ; 
        theta_ = theta_n ; 
end 

%}
% update history parameters:
%     P0 = P; 
%     TSigma0 = TSigma ; 
%     C0        = C;
    HH = HH_;    HH_2 = HH_2_ ;   HH_3 = HH_3_ ;      HH_4 = HH_4_ ;  
    dU0 = dU ;
%     theta_ = theta_n ;

    if rem(inc,50)==0
        save all.mat 
    end
    
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
%              [mu,~,~,~] = get_mu_2 ( dg_ , dQdsigma_2 , dt , inx , deltaY/2 , b_ , theta_n , 0);
             r2 = sqrt(J2)-P.*mu - eta_vp*lambda_dot ; 
             J_ = imag(r1)/hhhh;
             F5 = J_.\r2 ;
F5

diff(r2,lambda_dot)
