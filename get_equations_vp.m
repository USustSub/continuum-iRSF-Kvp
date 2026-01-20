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

    TSigma0 = HH(ij,idx*6-5:idx*6-2)';
    P0 = HH(ij,idx*6-1);
    C0 = HH(ij,idx*6-0); 
    Eps = HH_2(ij,idx*4-3:idx*4);
    sin_phi = HH_3(ij,idx);
    theta_= HH_4(ij,idx);

    % get strain
    if plast_it == 1 && inc == 1 
        [~,B,en] = get_data ( xx , node , di , form ) ;
        Shape{cc_ij}{1} = B ;  
        Shape{cc_ij}{2} = en ;  
    else
        B = Shape{cc_ij}{1};
        en = Shape{cc_ij}{2};
    end
    
    dEps_t = B*(dU(en)-dU0(en)); dEps_t = [dEps_t;0]; % add z strain
    dEps_0 = dEps_t; % save total strain  
        
    Eps_t = Eps' + [dEps_t(1) dEps_t(2) 0.5*dEps_t(3) dEps_t(4)]';
    Ekkc = Eps_t(1) + Eps_t(2) + Eps_t(4) ;
    Epsdc = Eps_t ; 
    Epsdc([1 2 4]) = Eps_t([1 2 4]) - 1/3*Ekkc;
    Eii  = sqrt(1/2*(Epsdc(1)).^2 + 1/2*(Epsdc(2)).^2 + 1/2*(Epsdc(4)).^2 + Epsdc(3).^2);
    
%     if plast_it > 1
%     Eii = Epsdc(1);
%     end
%         % Updat total strains (TO DO)
%         Exxc_t = Exxc + dExxc_t;
%         Eyyc_t = Eyyc + dEyyc_t;
%         Exyv_t = Exyv + dExyv_t;
%         Ezzc_t = Ezzc + 0;

    % get (visco)elastic tangent 
    [Dmat_t,Dmat,Gve,K_] = identify_tangent_v2 ( xx , mat , dynamic , node) ;
    dSigma = Dmat_t * dEps_t ;% total stress increment 
    dP    =-1/3*(sum(dSigma([1 2 4]))); % pressure increment

    Sigma_t = VE1* TSigma0 - [P0;P0;0;P0] + dSigma; % compute total trial stress (eq. 1)

    P   = P0 + dP; % get pressure
    TSigma = [P;P;0;P] + Sigma_t; % get deviatoric components
    J2 = 1/2*(TSigma(1)^2+TSigma(2)^2+TSigma(4)^2)+TSigma(3)^2;
    C     = C0; 
    F_     =  sqrt(J2) - C.*cos_phi - P.*sin_phi; % trial yield function
    % Plastic corrections
%     %{
    if F_>0 % if trial stresses exceed yield 
%         disp('plasticity')
        plastic = 1;
        dQdsigma = get_dQdSigma ( TSigma , sin_psi , J2 );
        h1 = cos_phi.*h.*sqrt(2/3).*sqrt((dQdsigma(1).^2 + 2*dQdsigma(3).^2 + dQdsigma(2).^2 + dQdsigma(4).^2));
        dg = (F_./(Gve + K_.*sin_phi*sin_psi + eta_vp/dt + h1 ));% increment of plastic multiplier
        dEps_t = [dEps_0(1); dEps_0(2) ; dEps_0(3); dEps_0(4)]  - ...
            dg*[dQdsigma(1);dQdsigma(2);dQdsigma(3)/1;dQdsigma(4)]; % compute corrected strain
        dSigma = Dmat_t * dEps_t ;% total stress increment 
        dP    =-1/3*(sum(dSigma([1 2 4]))); % pressure increment
        Sigma_t = VE1 .* TSigma0 - [P0;P0;0;P0] + dSigma; % compute total stresses
        P   = P0 + dP; % get pressure
        TSigma = [P;P;0;P] + Sigma_t; % get deviatoric components

        % Check yield function (F_vp=0)
        J2 = 1/2*(TSigma(1).^2+TSigma(2).^2+TSigma(4).^2)+TSigma(3).^2;
        dep=sqrt(2/3).*sqrt((dg.*dQdsigma(1)).^2+(dg.*dQdsigma(2)).^2+(dg.*dQdsigma(4)).^2+2*(dg.*dQdsigma(3)).^2); % eq. 

        C     = C0 + h.*dep; 
        F_     =  sqrt(J2) - (C.*cos_phi + P.*sin_phi) ;                       % Backbone plastic yield function
        F_vp   =  sqrt(J2) - (C.*cos_phi + P.*sin_phi) - (dg./dt).*eta_vp;  % Consistency viscoplasticity model
% %         fprintf('    Backbone EP max. Fc = %2.2e \n', max(F_(:)))
% %         fprintf('    Viscoplast. max. Fc = %2.2e \n', max(F_vp(:)))


        % now get the tangent operator:
        dQdsigma = get_dQdSigma ( TSigma , sin_psi , J2 )';
        dFdsigma = get_dQdSigma ( TSigma , sin_phi , J2 )';
        d2Qdsigma2 = get_dQ2dSigma2 ( TSigma  , J2 );
% %         [Dmat_t,~,Gve,K_] = identify_tangent_v2 ( xx , mat , dynamic , node) ;

        I = [ 1 0 0 0 ; 0 1 0 0 ; 0 0 1 0 ; 0 0 0 1] ;    
        M = I + dg * Dmat_t * d2Qdsigma2 ;
        Mi = inv(M) ; 
        Dvp = Mi * Dmat_t -  ( Mi * Dmat_t * dQdsigma * dFdsigma' * Mi * Dmat_t ) / (  eta_vp/dt + h1 + dFdsigma' * Mi * Dmat_t * dQdsigma ) ; 
        
        Dmat = Dvp(1:3,1:3); % use viscoplastic tangent
%         maxFvp = max(abs(F_vp),maxFvp);
    end % F_ > 0
%}
        maxF = max(abs(F_),maxF);
    
        kx = 0 ; ky = 0 ;
if idx == 1
    % x momentum equation - left
        Rx = Dmat(1,:)*B/deltaX  ;
        kx = kx - Rx ;
        Residual(2*ij-1) = Residual(2*ij-1) - Sigma_t(1)/deltaX;

    % y momentum equation - left
        Ry = Dmat(3,:)*B/deltaX  ;
        ky = ky - Ry ;
        Residual(2*ij  ) = Residual(2*ij  ) - Sigma_t(3)/deltaX;
end


if idx == 2
    % x momentum equation - right
        Rx = Dmat(1,:)*B/deltaX  ;
        kx = kx + Rx ;
        Residual(2*ij-1) = Residual(2*ij-1) + Sigma_t(1)/deltaX;

    % y momentum equation - right
        Ry = Dmat(3,:)*B/deltaX  ;
        ky = ky + Ry ;
        Residual(2*ij  ) = Residual(2*ij  ) + Sigma_t(3)/deltaX;
end


if idx == 3
    % x momentum equation - bot
        Rx = Dmat(3,:)*B/deltaY  ;
        kx = kx - Rx ;
        Residual(2*ij-1) = Residual(2*ij-1) - Sigma_t(3)/deltaY;

    % y momentum equation - bot
        Ry = Dmat(2,:)*B/deltaY  ;
        ky = ky - Ry ;
        Residual(2*ij  ) = Residual(2*ij  ) - Sigma_t(2)/deltaY;
end

if idx == 4
    % x momentum equation - top
        Rx = Dmat(3,:)*B/deltaY  ;
        kx = kx + Rx ;
        Residual(2*ij-1) = Residual(2*ij-1) + Sigma_t(3)/deltaY;

    % y momentum equation - top
        Ry = Dmat(2,:)*B/deltaY  ;
        ky = ky + Ry ; 
        Residual(2*ij  ) = Residual(2*ij  ) + Sigma_t(2)/deltaY;
end

% %     K(2*ij-1,en) = K(2*ij-1,en) + kx ; 
% %     K(2*ij  ,en) = K(2*ij  ,en) + ky ;
    
    [i j s] = find(kx); 
    Cell{c_} = [en(i)'*0+2*ij-1 ,en(j)', s'];
    [i j s] = find(ky); 
    Cell2{c_} = [en(i)'*0+2*ij ,en(j)', s'];
    c_ = c_ + 1 ; 

% add inertia terms 
if dynamic == 1 
    K(2*ij-1,2*ij-1) = K(2*ij-1,2*ij-1) - a0*rho;
    Converged_last_step = rho*a0*u0(2*ij-1) + rho*a2*vu0(2*ij-1) + rho*a4*au0(2*ij-1) ;
    Residual(2*ij-1) = Residual(2*ij-1) + Converged_last_step - a0*rho*dU(2*ij-1);

    K(2*ij-0,2*ij-0) = K(2*ij-0,2*ij-0) - a0*rho;
    Converged_last_step = rho*a0*u0(2*ij-0) + rho*a2*vu0(2*ij-0) + rho*a4*au0(2*ij-0) ;
    Residual(2*ij-0) = Residual(2*ij-0) + Converged_last_step - a0*rho*dU(2*ij-0);
end
    
% Residual(2*ij-1)
    HH_(ij,idx*6-5:idx*6-2) = TSigma' ;
    HH_(ij,idx*6-1) = P ;
    HH_(ij,idx*6-0) = C ; 
    HH_2_(ij,idx*4-3:idx*4) = Eps_t' ; 
%     SS = [SS ; xx [Eps_t(1:3)' log10(abs(Eii))] sqrt(J2)];
    SS(cc_ij,:) = [ xx [Eps_t(1:3)' log10(abs(Eii))] sqrt(J2) ] ;
    cc_ij = cc_ij + 1 ; 
    
    
if strcmp(Ex,'Herr')
% get Vp
%     sc = find(abs(node(:,2)-max(node(:,2))/2)<dY);
%     sc1 = sc(node(sc,2)>max(node(:,2))/2);
%     sc2 = sc(node(sc,2)<max(node(:,2))/2);
% %     plot(node(sc1,1),node(sc1,2),'gsq')
% %     plot(node(sc2,1),node(sc2,2),'ksq')
% %     Y_t = max(node(sc1,2));
% %     Y_b = max(node(sc2,2));
%     Y_t = D/2+25000 ; 
%     Y_b = D/2-25000 ; 
%     
%     xx_t = [xx(1) Y_t] ;
%     [index] = define_support(node,xx_t,di);
%     [phi,~,~] = MLS_ShapeFunction(xx_t,index,node,di,form);
%     u_t = phi*dU(index*2-1); 
%     v_t = phi*vu(index*2-1); 
%     
%     xx_b = [xx(1) Y_b] ;
%     [index] = define_support(node,xx_b,di);
%     [phi,~,~] = MLS_ShapeFunction(xx_b,index,node,di,form);
%     u_b = phi*dU(index*2-1); 
%     v_b = phi*vu(index*2-1); 

%     slip = u_t - u_b ; 
%     slip_rate_ = v_t - v_b ; 
    slip_rate_ = 2e-9*2 ; 
    slip = 0 ; 
    u_t = 0 ; u_b = 0 ; 
    Vp = slip_rate_ ;
% get b parameter
    if xx(1) < LL1 || xx(1) > LL2
        b_ = b1_; 
    else
        b_ = b2_ ; 
    end
        
% update theta 
    S1 = 1+theta_/dt;
    S2 = (1/dt + Vp/L_);
    theta_n = S1/S2;
% update friction coefficient
    param = Vp/2/V0*exp((mu0+b_*log(theta_n*V0/L_))/a_);
    mu_gp = a_*asinh(param) ;
    
    HH_3_(ij,idx) = mu_gp ; 
    HH_4_(ij,idx) = theta_n ; 
    Slip = [Slip ; xx u_t u_b slip mu_gp theta_n slip_rate_] ; 
end