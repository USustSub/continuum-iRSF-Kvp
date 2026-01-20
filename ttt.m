%%
clc
clear all; 
load test.mat
               
Exx = 0 ; 
Eyy = 0 ; 
Ezz = 0 ; 
Exy = 0 ;
Txx0 =0  ;
Tyy0 = 0 ; 
Tzz0 = 0; 
Txy0 = 0 ; 
%
% get strain increments using dUxc_t (TO DO)
        [~,B,en] = get_data ( x_l , node , di , form ) ;
        [Dmat_t,Dmat,Gve,K_] = identify_tangent_v2 ( x_l , mat) ;
        dExx_t  = B(1,:)*dU(en);
        dEyy_t  = B(2,:)*dU(en);
        dEzz_t  = 0; 
        dExy_t  = 1*( B(3,:)*dU(en) );
%         [ dExx_t dEyy_t dExy_t dEzz_t]'
        % Updat total strains (TO DO)
        Exx_t = Exx + dExx_t;
        Eyy_t = Eyy + dEyy_t;
        Exy_t = Exy + dExy_t;
        Ezz_t = Ezz + 0;
        % Save strains 
        dExx_0 = dExx_t;    
        dEyy_0 = dEyy_t;  
        dExy_0 = dExy_t;
        dEzz_0 = dEzz_t; 
        % Total stress increments from total strain increment:
        dSxx = (K_ + 4/3*Gve).*dExx_t + (K_ - 2/3*Gve).*dEyy_t + (K_ - 2/3*Gve).*dEzz_t;
        dSyy = (K_ + 4/3*Gve).*dEyy_t + (K_ - 2/3*Gve).*dExx_t + (K_ - 2/3*Gve).*dEzz_t;
        dSzz = (K_ + 4/3*Gve).*dEzz_t + (K_ - 2/3*Gve).*dExx_t + (K_ - 2/3*Gve).*dEyy_t;
        dTxy = 1*Gve.*dExy_t ;
%         [dSxx dSyy dTxy dSzz]'
        % Total stresses
        dP    =-1/3*(dSxx + dSyy + dSzz);
        
        Sxx_t = VE1 .* Txx0 - P0 + dSxx;
        Syy_t = VE1 .* Tyy0 - P0 + dSyy;
        Txy_t = VE1 .* Txy0       + dTxy;
        Szz_t = VE1 .* Tzz0 - P0 + dSzz;
%         [Sxx_t Syy_t Txy_t Szz_t]'
        % Compute pressure and deviatoric components
        P   = P0 + dP;
        Txx = P + Sxx_t; 
        Tyy = P + Syy_t; 
        Tzz = P + Szz_t; 
        Txy = Txy_t;
        
        % Check yield function
        J2    = 1/2*(Txx.^2 + Tyy.^2 + Tzz.^2) + Txy.^2;
        Cc     = C0; 
        syc    = Cc.*cos_phi + (P).*sin_phi;
        F     =  sqrt(J2) - syc;
%         F
        fprintf('max. trial     F = %2.2e\n', max(F(:)))
        % Plastic corrections
%         if max(F(:))>0 
            % Plastic flags
            pl = F>0;
            % dQds centers
            txx = P+Sxx_t; 
            tyy = P+Syy_t;
            txy = Txy_t; 
            tzz = P+Szz_t; 
            dQdsxx = 1 * txx .* (1 * J2) .^ (-0.1e1 / 0.2e1) / 0.2e1 + sin_psi / 0.3e1;
            dQdsyy = 1 * tyy .* (1 * J2) .^ (-0.1e1 / 0.2e1) / 0.2e1 + sin_psi / 0.3e1;
            dQdsxy = 1 * txy .* (1 * J2) .^ (-0.1e1 / 0.2e1);
            dQdszz = 1 * tzz .* (1 * J2) .^ (-0.1e1 / 0.2e1) / 0.2e1 + sin_psi / 0.3e1;
%             [dQdsxx dQdsyy dQdsxy dQdszz ]'
            % dgamma - increment of plastic multiplier
            h1 = cos_phi.*h.*sqrt(2/3).*sqrt((dQdsxx.^2 + 2*dQdsxy.^2 + dQdsyy.^2 + dQdszz.^2));
            dg = (F./(1.*Gve + K_.*sin_phi*sin_psi + eta_vp/dt + h1 ));
%             dg
            % Corrected strain on centers
            dExx_t = dExx_0 - 1.*dg.*dQdsxx;
            dEyy_t = dEyy_0 - 1.*dg.*dQdsyy;
            dExy_t = dExy_0 - 1.*dg.*dQdsxy/1;
            dEzz_t = dEzz_0 - 1.*dg.*dQdszz;
%             [dExy_0]'-dg*dQdsxy/2
%             [dQdsxx dQdsyy dQdsxy dQdszz ]'
            
            % Total stress increments
            dSxx = (K_ + 4/3*Gve).*dExx_t + (K_ - 2/3*Gve).*dEyy_t + (K_ - 2/3*Gve).*dEzz_t;
            dSyy = (K_ + 4/3*Gve).*dEyy_t + (K_ - 2/3*Gve).*dExx_t + (K_ - 2/3*Gve).*dEzz_t;
            dSzz = (K_ + 4/3*Gve).*dEzz_t + (K_ - 2/3*Gve).*dExx_t + (K_ - 2/3*Gve).*dEyy_t;
            dTxy = 1*Gve.*dExy_t ;
            % Total stresses
            dP    =-1/3*(dSxx + dSyy + dSzz);
            Sxx_t = VE1 .* Txx0 - P0 + dSxx;
            Syy_t = VE1 .* Tyy0 - P0 + dSyy;
            Txy_t = VE1 .* Txy0       + dTxy;
            Szz_t = VE1 .* Tzz0 - P0 + dSzz;
            % Compute pressure and deviatoric components
            P     = P0 + dP;
            Txx   = P + Sxx_t; Tyy = P + Syy_t; Tzz = P + Szz_t; Txy = Txy_t;
            % Check yield function
            J2    = 1/2*(Txx.^2 + Tyy.^2 + Tzz.^2) + Txy.^2;
            dep   = sqrt(2/3).*sqrt((dg.*dQdsxx).^2+(dg.*dQdsyy).^2+(dg.*dQdszz).^2+2*(dg.*dQdsxy).^2);
            C     = C0 + 1*h.*dep;
            syc    = C.*cos_phi + (P).*sin_phi;
            F     =  sqrt(1*J2) - syc;                       % Backbone plastic yield function
%             gamdotc = plc.*(sqrt(1*J2c) - syc) ./ eta_vp; % plastic strain rate
            %             Fc_vp   =  sqrt(1*J2c) - syc - gamdotc.*eta_vp; % equivalent to line below
            %             Fv_vp   =  sqrt(1*J2v) - syv - gamdotv.*eta_vp;
            F_vp   =  sqrt(1*J2) - syc - (dg./dt).*eta_vp;  % Consistency viscoplasticity model
%             sqrt(1*J2)
%             syc
%             (dg./dt).*eta_vp
            fprintf('    Backbone EP max. Fc = %2.2e \n', max(F(:)))
            fprintf('    Viscoplast. max. Fc = %2.2e \n', max(F_vp(:)))
%         end
        

            txx = P+Sxx_t; 
            tyy = P+Syy_t;
            txy = Txy_t; 
            tzz = P+Szz_t; 

dQdsxx = 1 * txx .* (1 * J2) .^ (-0.1e1 / 0.2e1) / 0.2e1 + sin_psi / 0.3e1;
dQdsyy = 1 * tyy .* (1 * J2) .^ (-0.1e1 / 0.2e1) / 0.2e1 + sin_psi / 0.3e1;
dQdsxy = 1 * 1 * txy .* (1 * J2) .^ (-0.1e1 / 0.2e1);
dQdszz = 1 * tzz .* (1 * J2) .^ (-0.1e1 / 0.2e1) / 0.2e1 + sin_psi / 0.3e1;
% Reduced elastic stiffness matrix
    d2Qdsxxdsxx = 1 .* (1 .* J2) .^ (-0.1e1 ./ 0.2e1) ./ 0.3e1 - 1 .^ 2 .* txx .^ 2 .* (1 .* J2) .^ (-0.3e1 ./ 0.2e1) ./ 0.4e1;
    d2Qdsxxdsyy = -1 .* (1 .* J2) .^ (-0.1e1 ./ 0.2e1) ./ 0.6e1 - 1 .^ 2 .* txx .* tyy .* (1 .* J2) .^ (-0.3e1 ./ 0.2e1) ./ 0.4e1;
    d2Qdsxxdsxy = -1 .^ 2 .* txx .* txy .* (1 .* J2) .^ (-0.3e1 ./ 0.2e1) ./ 0.2e1;
    d2Qdsxxdszz = -1 .* (1 .* J2) .^ (-0.1e1 ./ 0.2e1) ./ 0.6e1 - 1 .^ 2 .* txx .* tzz .* (1 .* J2) .^ (-0.3e1 ./ 0.2e1) ./ 0.4e1;
    d2Qdsyydsxx = -1 .* (1 .* J2) .^ (-0.1e1 ./ 0.2e1) ./ 0.6e1 - 1 .^ 2 .* txx .* tyy .* (1 .* J2) .^ (-0.3e1 ./ 0.2e1) ./ 0.4e1;
    d2Qdsyydsyy = 1 .* (1 .* J2) .^ (-0.1e1 ./ 0.2e1) ./ 0.3e1 - 1 .^ 2 .* tyy .^ 2 .* (1 .* J2) .^ (-0.3e1 ./ 0.2e1) ./ 0.4e1;
    d2Qdsyydsxy = -1 .^ 2 .* tyy .* txy .* (1 .* J2) .^ (-0.3e1 ./ 0.2e1) ./ 0.2e1;
    d2Qdsyydszz = -1 .* (1 .* J2) .^ (-0.1e1 ./ 0.2e1) ./ 0.6e1 - 1 .^ 2 .* tyy .* tzz .* (1 .* J2) .^ (-0.3e1 ./ 0.2e1) ./ 0.4e1;
    d2Qdsxydsxx = -1 .^ 2 .* txx .* txy .* (1 .* J2) .^ (-0.3e1 ./ 0.2e1) ./ 0.2e1;
    d2Qdsxydsyy = -1 .^ 2 .* tyy .* txy .* (1 .* J2) .^ (-0.3e1 ./ 0.2e1) ./ 0.2e1;
    d2Qdsxydsxy = 1 .* (1 .* J2) .^ (-0.1e1 ./ 0.2e1) - 1 .^ 2 .* txy .^ 2 .* (1 .* J2) .^ (-0.3e1 ./ 0.2e1);
    d2Qdsxydszz = -1 .^ 2 .* txy .* tzz .* (1 .* J2) .^ (-0.3e1 ./ 0.2e1) ./ 0.2e1;
    d2Qdszzdsxx = -1 .* (1 .* J2) .^ (-0.1e1 ./ 0.2e1) ./ 0.6e1 - 1 .^ 2 .* txx .* tzz .* (1 .* J2) .^ (-0.3e1 ./ 0.2e1) ./ 0.4e1;
    d2Qdszzdsyy = -1 .* (1 .* J2) .^ (-0.1e1 ./ 0.2e1) ./ 0.6e1 - 1 .^ 2 .* tyy .* tzz .* (1 .* J2) .^ (-0.3e1 ./ 0.2e1) ./ 0.4e1;
    d2Qdszzdsxy = -1 .^ 2 .* txy .* tzz .* (1 .* J2) .^ (-0.3e1 ./ 0.2e1) ./ 0.2e1;
    d2Qdszzdszz = 1 .* (1 .* J2) .^ (-0.1e1 ./ 0.2e1) ./ 0.3e1 - 1 .^ 2 .* tzz .^ 2 .* (1 .* J2) .^ (-0.3e1 ./ 0.2e1) ./ 0.4e1;
    
    d2QdSigma2 =  [
                    d2Qdsxxdsxx d2Qdsxxdsyy d2Qdsxxdsxy d2Qdsxxdszz  
                    d2Qdsyydsxx d2Qdsyydsyy d2Qdsyydsxy d2Qdsyydszz  
                    d2Qdsxydsxx d2Qdsxydsyy d2Qdsxydsxy d2Qdsxydszz  
                    d2Qdszzdsxx d2Qdszzdsyy d2Qdszzdsxy d2Qdszzdszz  
                    ];
% dQds centers
dQdsxx = 1 * txx .* (1 * J2) .^ (-0.1e1 / 0.2e1) / 0.2e1 + sin_psi / 0.3e1;
dQdsyy = 1 * tyy .* (1 * J2) .^ (-0.1e1 / 0.2e1) / 0.2e1 + sin_psi / 0.3e1;
dQdsxy = 1 * 1 * txy .* (1 * J2) .^ (-0.1e1 / 0.2e1);
dQdszz = 1 * tzz .* (1 * J2) .^ (-0.1e1 / 0.2e1) / 0.2e1 + sin_psi / 0.3e1;
% dFds centers
dFdsxx = 1 * txx .* (1 * J2) .^ (-0.1e1 / 0.2e1) / 0.2e1 + sin_phi / 0.3e1;
dFdsyy = 1 * tyy .* (1 * J2) .^ (-0.1e1 / 0.2e1) / 0.2e1 + sin_phi / 0.3e1;
dFdsxy = 1 * 1 * txy .* (1 * J2) .^ (-0.1e1 / 0.2e1);
dFdszz = 1 * tzz .* (1 * J2) .^ (-0.1e1 / 0.2e1) / 0.2e1 + sin_phi / 0.3e1;
% [dFdsxx dFdsyy dFdsxy dFdszz]'

De11c    = K_ + 4/3.*Gve;  De12c    = K_ - 2/3.*Gve;  De13c    = 0*K_;  De14c    = K_ - 2/3.*Gve;
De21c    = K_ - 2/3.*Gve;  De22c    = K_ + 4/3.*Gve;  De23c    = 0*K_;  De24c    = K_ - 2/3.*Gve;
De31c    = 0.*Gve;         De32c    = 0.*Gve;         De33c    = 1.*Gve; De34c    = 0.*Gve;
De41c    = K_ - 2/3.*Gve;  De42c    = K_ - 2/3.*Gve;  De43c    = 0*K_;  De44c    = K_ + 4/3.*Gve;
dlam = dg ;
    D11 = De11c; D12 = De12c; D13 = De13c; D14 = De14c;
    D21 = De21c; D22 = De22c; D23 = De23c; D24 = De24c;
    D31 = De31c; D32 = De32c; D33 = De33c; D34 = De34c;
    D41 = De41c; D42 = De42c; D43 = De43c; D44 = De44c;
    M11 = 1 + dlam .* (D11 .* d2Qdsxxdsxx + D12 .* d2Qdsyydsxx + D14 .* d2Qdszzdsxx);
    M12 = dlam .* (D11 .* d2Qdsxxdsyy + D12 .* d2Qdsyydsyy + D14 .* d2Qdszzdsyy);
    M13 = dlam .* (D11 .* d2Qdsxxdsxy + D12 .* d2Qdsyydsxy + D14 .* d2Qdszzdsxy);
    M14 = dlam .* (D11 .* d2Qdsxxdszz + D12 .* d2Qdsyydszz + D14 .* d2Qdszzdszz);
    M21 = dlam .* (D21 .* d2Qdsxxdsxx + D22 .* d2Qdsyydsxx + D24 .* d2Qdszzdsxx);
    M22 = 1 + dlam .* (D21 .* d2Qdsxxdsyy + D22 .* d2Qdsyydsyy + D24 .* d2Qdszzdsyy);
    M23 = dlam .* (D21 .* d2Qdsxxdsxy + D22 .* d2Qdsyydsxy + D24 .* d2Qdszzdsxy);
    M24 = dlam .* (D21 .* d2Qdsxxdszz + D22 .* d2Qdsyydszz + D24 .* d2Qdszzdszz);
    M31 = dlam .* D33 .* d2Qdsxydsxx;
    M32 = dlam .* D33 .* d2Qdsxydsyy;
    M33 = 1 + dlam .* D33 .* d2Qdsxydsxy;
    M34 = dlam .* D33 .* d2Qdsxydszz;
    M41 = dlam .* (D41 .* d2Qdsxxdsxx + D42 .* d2Qdsyydsxx + D44 .* d2Qdszzdsxx);
    M42 = dlam .* (D41 .* d2Qdsxxdsyy + D42 .* d2Qdsyydsyy + D44 .* d2Qdszzdsyy);
    M43 = dlam .* (D41 .* d2Qdsxxdsxy + D42 .* d2Qdsyydsxy + D44 .* d2Qdszzdsxy);
    M44 = 1 + dlam .* (D41 .* d2Qdsxxdszz + D42 .* d2Qdsyydszz + D44 .* d2Qdszzdszz);
    [M11 M12 M13 M14  ;
     M21 M22 M23 M24  ;
     M31 M32 M33 M34  ;
     M41 M42 M43 M44  ]
 
    [ Mi11,Mi12,Mi13,Mi14, Mi21,Mi22,Mi23,Mi24, Mi31,Mi32,Mi33,Mi34, Mi41,Mi42,Mi43,Mi44  ] = M2Di_EP5_inverse_4x4( M11,M12,M13,M14, M21,M22,M23,M24, M31,M32,M33,M34, M41,M42,M43,M44  );
    den = (dFdsxx .* Mi11 + dFdsyy .* Mi21 + dFdsxy .* Mi31 + dFdszz .* Mi41) .* (D11 .* dQdsxx + D12 .* dQdsyy + D14 .* dQdszz) + (dFdsxx .* Mi12 + dFdsyy .* Mi22 + dFdsxy .* Mi32 + dFdszz .* Mi42) .* (D21 .* dQdsxx + D22 .* dQdsyy + D24 .* dQdszz) + (dFdsxx .* Mi13 + dFdsyy .* Mi23 + dFdsxy .* Mi33 + dFdszz .* Mi43) .* D33 .* dQdsxy + (dFdsxx .* Mi14 + dFdsyy .* Mi24 + dFdsxy .* Mi34 + dFdszz .* Mi44) .* (D41 .* dQdsxx + D42 .* dQdsyy + D44 .* dQdszz);;
    den = den + 1.*eta_vp/dt + h1;
    Dep11c = Mi11 .* (D11 .* (1 - 1 ./ den .* (dQdsxx .* dFdsxx .* (Mi11 .* D11 + Mi12 .* D21 + Mi14 .* D41) + dQdsxx .* dFdsyy .* (Mi21 .* D11 + Mi22 .* D21 + Mi24 .* D41) + dQdsxx .* dFdsxy .* (Mi31 .* D11 + Mi32 .* D21 + Mi34 .* D41) + dQdsxx .* dFdszz .* (Mi41 .* D11 + Mi42 .* D21 + Mi44 .* D41))) - D12 ./ den .* (dQdsyy .* dFdsxx .* (Mi11 .* D11 + Mi12 .* D21 + Mi14 .* D41) + dQdsyy .* dFdsyy .* (Mi21 .* D11 + Mi22 .* D21 + Mi24 .* D41) + dQdsyy .* dFdsxy .* (Mi31 .* D11 + Mi32 .* D21 + Mi34 .* D41) + dQdsyy .* dFdszz .* (Mi41 .* D11 + Mi42 .* D21 + Mi44 .* D41)) - D14 ./ den .* (dQdszz .* dFdsxx .* (Mi11 .* D11 + Mi12 .* D21 + Mi14 .* D41) + dQdszz .* dFdsyy .* (Mi21 .* D11 + Mi22 .* D21 + Mi24 .* D41) + dQdszz .* dFdsxy .* (Mi31 .* D11 + Mi32 .* D21 + Mi34 .* D41) + dQdszz .* dFdszz .* (Mi41 .* D11 + Mi42 .* D21 + Mi44 .* D41))) + Mi12 .* (D21 .* (1 - 1 ./ den .* (dQdsxx .* dFdsxx .* (Mi11 .* D11 + Mi12 .* D21 + Mi14 .* D41) + dQdsxx .* dFdsyy .* (Mi21 .* D11 + Mi22 .* D21 + Mi24 .* D41) + dQdsxx .* dFdsxy .* (Mi31 .* D11 + Mi32 .* D21 + Mi34 .* D41) + dQdsxx .* dFdszz .* (Mi41 .* D11 + Mi42 .* D21 + Mi44 .* D41))) - D22 ./ den .* (dQdsyy .* dFdsxx .* (Mi11 .* D11 + Mi12 .* D21 + Mi14 .* D41) + dQdsyy .* dFdsyy .* (Mi21 .* D11 + Mi22 .* D21 + Mi24 .* D41) + dQdsyy .* dFdsxy .* (Mi31 .* D11 + Mi32 .* D21 + Mi34 .* D41) + dQdsyy .* dFdszz .* (Mi41 .* D11 + Mi42 .* D21 + Mi44 .* D41)) - D24 ./ den .* (dQdszz .* dFdsxx .* (Mi11 .* D11 + Mi12 .* D21 + Mi14 .* D41) + dQdszz .* dFdsyy .* (Mi21 .* D11 + Mi22 .* D21 + Mi24 .* D41) + dQdszz .* dFdsxy .* (Mi31 .* D11 + Mi32 .* D21 + Mi34 .* D41) + dQdszz .* dFdszz .* (Mi41 .* D11 + Mi42 .* D21 + Mi44 .* D41))) - Mi13 .* D33 ./ den .* (dQdsxy .* dFdsxx .* (Mi11 .* D11 + Mi12 .* D21 + Mi14 .* D41) + dQdsxy .* dFdsyy .* (Mi21 .* D11 + Mi22 .* D21 + Mi24 .* D41) + dQdsxy .* dFdsxy .* (Mi31 .* D11 + Mi32 .* D21 + Mi34 .* D41) + dQdsxy .* dFdszz .* (Mi41 .* D11 + Mi42 .* D21 + Mi44 .* D41)) + Mi14 .* (D41 .* (1 - 1 ./ den .* (dQdsxx .* dFdsxx .* (Mi11 .* D11 + Mi12 .* D21 + Mi14 .* D41) + dQdsxx .* dFdsyy .* (Mi21 .* D11 + Mi22 .* D21 + Mi24 .* D41) + dQdsxx .* dFdsxy .* (Mi31 .* D11 + Mi32 .* D21 + Mi34 .* D41) + dQdsxx .* dFdszz .* (Mi41 .* D11 + Mi42 .* D21 + Mi44 .* D41))) - D42 ./ den .* (dQdsyy .* dFdsxx .* (Mi11 .* D11 + Mi12 .* D21 + Mi14 .* D41) + dQdsyy .* dFdsyy .* (Mi21 .* D11 + Mi22 .* D21 + Mi24 .* D41) + dQdsyy .* dFdsxy .* (Mi31 .* D11 + Mi32 .* D21 + Mi34 .* D41) + dQdsyy .* dFdszz .* (Mi41 .* D11 + Mi42 .* D21 + Mi44 .* D41)) - D44 ./ den .* (dQdszz .* dFdsxx .* (Mi11 .* D11 + Mi12 .* D21 + Mi14 .* D41) + dQdszz .* dFdsyy .* (Mi21 .* D11 + Mi22 .* D21 + Mi24 .* D41) + dQdszz .* dFdsxy .* (Mi31 .* D11 + Mi32 .* D21 + Mi34 .* D41) + dQdszz .* dFdszz .* (Mi41 .* D11 + Mi42 .* D21 + Mi44 .* D41)));
    Dep12c = Mi11 .* (-D11 ./ den .* (dQdsxx .* dFdsxx .* (Mi11 .* D12 + Mi12 .* D22 + Mi14 .* D42) + dQdsxx .* dFdsyy .* (Mi21 .* D12 + Mi22 .* D22 + Mi24 .* D42) + dQdsxx .* dFdsxy .* (Mi31 .* D12 + Mi32 .* D22 + Mi34 .* D42) + dQdsxx .* dFdszz .* (Mi41 .* D12 + Mi42 .* D22 + Mi44 .* D42)) + D12 .* (1 - 1 ./ den .* (dQdsyy .* dFdsxx .* (Mi11 .* D12 + Mi12 .* D22 + Mi14 .* D42) + dQdsyy .* dFdsyy .* (Mi21 .* D12 + Mi22 .* D22 + Mi24 .* D42) + dQdsyy .* dFdsxy .* (Mi31 .* D12 + Mi32 .* D22 + Mi34 .* D42) + dQdsyy .* dFdszz .* (Mi41 .* D12 + Mi42 .* D22 + Mi44 .* D42))) - D14 ./ den .* (dQdszz .* dFdsxx .* (Mi11 .* D12 + Mi12 .* D22 + Mi14 .* D42) + dQdszz .* dFdsyy .* (Mi21 .* D12 + Mi22 .* D22 + Mi24 .* D42) + dQdszz .* dFdsxy .* (Mi31 .* D12 + Mi32 .* D22 + Mi34 .* D42) + dQdszz .* dFdszz .* (Mi41 .* D12 + Mi42 .* D22 + Mi44 .* D42))) + Mi12 .* (-D21 ./ den .* (dQdsxx .* dFdsxx .* (Mi11 .* D12 + Mi12 .* D22 + Mi14 .* D42) + dQdsxx .* dFdsyy .* (Mi21 .* D12 + Mi22 .* D22 + Mi24 .* D42) + dQdsxx .* dFdsxy .* (Mi31 .* D12 + Mi32 .* D22 + Mi34 .* D42) + dQdsxx .* dFdszz .* (Mi41 .* D12 + Mi42 .* D22 + Mi44 .* D42)) + D22 .* (1 - 1 ./ den .* (dQdsyy .* dFdsxx .* (Mi11 .* D12 + Mi12 .* D22 + Mi14 .* D42) + dQdsyy .* dFdsyy .* (Mi21 .* D12 + Mi22 .* D22 + Mi24 .* D42) + dQdsyy .* dFdsxy .* (Mi31 .* D12 + Mi32 .* D22 + Mi34 .* D42) + dQdsyy .* dFdszz .* (Mi41 .* D12 + Mi42 .* D22 + Mi44 .* D42))) - D24 ./ den .* (dQdszz .* dFdsxx .* (Mi11 .* D12 + Mi12 .* D22 + Mi14 .* D42) + dQdszz .* dFdsyy .* (Mi21 .* D12 + Mi22 .* D22 + Mi24 .* D42) + dQdszz .* dFdsxy .* (Mi31 .* D12 + Mi32 .* D22 + Mi34 .* D42) + dQdszz .* dFdszz .* (Mi41 .* D12 + Mi42 .* D22 + Mi44 .* D42))) - Mi13 .* D33 ./ den .* (dQdsxy .* dFdsxx .* (Mi11 .* D12 + Mi12 .* D22 + Mi14 .* D42) + dQdsxy .* dFdsyy .* (Mi21 .* D12 + Mi22 .* D22 + Mi24 .* D42) + dQdsxy .* dFdsxy .* (Mi31 .* D12 + Mi32 .* D22 + Mi34 .* D42) + dQdsxy .* dFdszz .* (Mi41 .* D12 + Mi42 .* D22 + Mi44 .* D42)) + Mi14 .* (-D41 ./ den .* (dQdsxx .* dFdsxx .* (Mi11 .* D12 + Mi12 .* D22 + Mi14 .* D42) + dQdsxx .* dFdsyy .* (Mi21 .* D12 + Mi22 .* D22 + Mi24 .* D42) + dQdsxx .* dFdsxy .* (Mi31 .* D12 + Mi32 .* D22 + Mi34 .* D42) + dQdsxx .* dFdszz .* (Mi41 .* D12 + Mi42 .* D22 + Mi44 .* D42)) + D42 .* (1 - 1 ./ den .* (dQdsyy .* dFdsxx .* (Mi11 .* D12 + Mi12 .* D22 + Mi14 .* D42) + dQdsyy .* dFdsyy .* (Mi21 .* D12 + Mi22 .* D22 + Mi24 .* D42) + dQdsyy .* dFdsxy .* (Mi31 .* D12 + Mi32 .* D22 + Mi34 .* D42) + dQdsyy .* dFdszz .* (Mi41 .* D12 + Mi42 .* D22 + Mi44 .* D42))) - D44 ./ den .* (dQdszz .* dFdsxx .* (Mi11 .* D12 + Mi12 .* D22 + Mi14 .* D42) + dQdszz .* dFdsyy .* (Mi21 .* D12 + Mi22 .* D22 + Mi24 .* D42) + dQdszz .* dFdsxy .* (Mi31 .* D12 + Mi32 .* D22 + Mi34 .* D42) + dQdszz .* dFdszz .* (Mi41 .* D12 + Mi42 .* D22 + Mi44 .* D42)));
    Dep13c = Mi11 .* (-D11 ./ den .* (dQdsxx .* dFdsxx .* Mi13 .* D33 + dQdsxx .* dFdsyy .* Mi23 .* D33 + dQdsxx .* dFdsxy .* Mi33 .* D33 + dQdsxx .* dFdszz .* Mi43 .* D33) - D12 ./ den .* (dQdsyy .* dFdsxx .* Mi13 .* D33 + dQdsyy .* dFdsyy .* Mi23 .* D33 + dQdsyy .* dFdsxy .* Mi33 .* D33 + dQdsyy .* dFdszz .* Mi43 .* D33) - D14 ./ den .* (dQdszz .* dFdsxx .* Mi13 .* D33 + dQdszz .* dFdsyy .* Mi23 .* D33 + dQdszz .* dFdsxy .* Mi33 .* D33 + dQdszz .* dFdszz .* Mi43 .* D33)) + Mi12 .* (-D21 ./ den .* (dQdsxx .* dFdsxx .* Mi13 .* D33 + dQdsxx .* dFdsyy .* Mi23 .* D33 + dQdsxx .* dFdsxy .* Mi33 .* D33 + dQdsxx .* dFdszz .* Mi43 .* D33) - D22 ./ den .* (dQdsyy .* dFdsxx .* Mi13 .* D33 + dQdsyy .* dFdsyy .* Mi23 .* D33 + dQdsyy .* dFdsxy .* Mi33 .* D33 + dQdsyy .* dFdszz .* Mi43 .* D33) - D24 ./ den .* (dQdszz .* dFdsxx .* Mi13 .* D33 + dQdszz .* dFdsyy .* Mi23 .* D33 + dQdszz .* dFdsxy .* Mi33 .* D33 + dQdszz .* dFdszz .* Mi43 .* D33)) + Mi13 .* D33 .* (1 - 1 ./ den .* (dQdsxy .* dFdsxx .* Mi13 .* D33 + dQdsxy .* dFdsyy .* Mi23 .* D33 + dQdsxy .* dFdsxy .* Mi33 .* D33 + dQdsxy .* dFdszz .* Mi43 .* D33)) + Mi14 .* (-D41 ./ den .* (dQdsxx .* dFdsxx .* Mi13 .* D33 + dQdsxx .* dFdsyy .* Mi23 .* D33 + dQdsxx .* dFdsxy .* Mi33 .* D33 + dQdsxx .* dFdszz .* Mi43 .* D33) - D42 ./ den .* (dQdsyy .* dFdsxx .* Mi13 .* D33 + dQdsyy .* dFdsyy .* Mi23 .* D33 + dQdsyy .* dFdsxy .* Mi33 .* D33 + dQdsyy .* dFdszz .* Mi43 .* D33) - D44 ./ den .* (dQdszz .* dFdsxx .* Mi13 .* D33 + dQdszz .* dFdsyy .* Mi23 .* D33 + dQdszz .* dFdsxy .* Mi33 .* D33 + dQdszz .* dFdszz .* Mi43 .* D33));
    Dep21c = Mi21 .* (D11 .* (1 - 1 ./ den .* (dQdsxx .* dFdsxx .* (Mi11 .* D11 + Mi12 .* D21 + Mi14 .* D41) + dQdsxx .* dFdsyy .* (Mi21 .* D11 + Mi22 .* D21 + Mi24 .* D41) + dQdsxx .* dFdsxy .* (Mi31 .* D11 + Mi32 .* D21 + Mi34 .* D41) + dQdsxx .* dFdszz .* (Mi41 .* D11 + Mi42 .* D21 + Mi44 .* D41))) - D12 ./ den .* (dQdsyy .* dFdsxx .* (Mi11 .* D11 + Mi12 .* D21 + Mi14 .* D41) + dQdsyy .* dFdsyy .* (Mi21 .* D11 + Mi22 .* D21 + Mi24 .* D41) + dQdsyy .* dFdsxy .* (Mi31 .* D11 + Mi32 .* D21 + Mi34 .* D41) + dQdsyy .* dFdszz .* (Mi41 .* D11 + Mi42 .* D21 + Mi44 .* D41)) - D14 ./ den .* (dQdszz .* dFdsxx .* (Mi11 .* D11 + Mi12 .* D21 + Mi14 .* D41) + dQdszz .* dFdsyy .* (Mi21 .* D11 + Mi22 .* D21 + Mi24 .* D41) + dQdszz .* dFdsxy .* (Mi31 .* D11 + Mi32 .* D21 + Mi34 .* D41) + dQdszz .* dFdszz .* (Mi41 .* D11 + Mi42 .* D21 + Mi44 .* D41))) + Mi22 .* (D21 .* (1 - 1 ./ den .* (dQdsxx .* dFdsxx .* (Mi11 .* D11 + Mi12 .* D21 + Mi14 .* D41) + dQdsxx .* dFdsyy .* (Mi21 .* D11 + Mi22 .* D21 + Mi24 .* D41) + dQdsxx .* dFdsxy .* (Mi31 .* D11 + Mi32 .* D21 + Mi34 .* D41) + dQdsxx .* dFdszz .* (Mi41 .* D11 + Mi42 .* D21 + Mi44 .* D41))) - D22 ./ den .* (dQdsyy .* dFdsxx .* (Mi11 .* D11 + Mi12 .* D21 + Mi14 .* D41) + dQdsyy .* dFdsyy .* (Mi21 .* D11 + Mi22 .* D21 + Mi24 .* D41) + dQdsyy .* dFdsxy .* (Mi31 .* D11 + Mi32 .* D21 + Mi34 .* D41) + dQdsyy .* dFdszz .* (Mi41 .* D11 + Mi42 .* D21 + Mi44 .* D41)) - D24 ./ den .* (dQdszz .* dFdsxx .* (Mi11 .* D11 + Mi12 .* D21 + Mi14 .* D41) + dQdszz .* dFdsyy .* (Mi21 .* D11 + Mi22 .* D21 + Mi24 .* D41) + dQdszz .* dFdsxy .* (Mi31 .* D11 + Mi32 .* D21 + Mi34 .* D41) + dQdszz .* dFdszz .* (Mi41 .* D11 + Mi42 .* D21 + Mi44 .* D41))) - Mi23 .* D33 ./ den .* (dQdsxy .* dFdsxx .* (Mi11 .* D11 + Mi12 .* D21 + Mi14 .* D41) + dQdsxy .* dFdsyy .* (Mi21 .* D11 + Mi22 .* D21 + Mi24 .* D41) + dQdsxy .* dFdsxy .* (Mi31 .* D11 + Mi32 .* D21 + Mi34 .* D41) + dQdsxy .* dFdszz .* (Mi41 .* D11 + Mi42 .* D21 + Mi44 .* D41)) + Mi24 .* (D41 .* (1 - 1 ./ den .* (dQdsxx .* dFdsxx .* (Mi11 .* D11 + Mi12 .* D21 + Mi14 .* D41) + dQdsxx .* dFdsyy .* (Mi21 .* D11 + Mi22 .* D21 + Mi24 .* D41) + dQdsxx .* dFdsxy .* (Mi31 .* D11 + Mi32 .* D21 + Mi34 .* D41) + dQdsxx .* dFdszz .* (Mi41 .* D11 + Mi42 .* D21 + Mi44 .* D41))) - D42 ./ den .* (dQdsyy .* dFdsxx .* (Mi11 .* D11 + Mi12 .* D21 + Mi14 .* D41) + dQdsyy .* dFdsyy .* (Mi21 .* D11 + Mi22 .* D21 + Mi24 .* D41) + dQdsyy .* dFdsxy .* (Mi31 .* D11 + Mi32 .* D21 + Mi34 .* D41) + dQdsyy .* dFdszz .* (Mi41 .* D11 + Mi42 .* D21 + Mi44 .* D41)) - D44 ./ den .* (dQdszz .* dFdsxx .* (Mi11 .* D11 + Mi12 .* D21 + Mi14 .* D41) + dQdszz .* dFdsyy .* (Mi21 .* D11 + Mi22 .* D21 + Mi24 .* D41) + dQdszz .* dFdsxy .* (Mi31 .* D11 + Mi32 .* D21 + Mi34 .* D41) + dQdszz .* dFdszz .* (Mi41 .* D11 + Mi42 .* D21 + Mi44 .* D41)));
    Dep22c = Mi21 .* (-D11 ./ den .* (dQdsxx .* dFdsxx .* (Mi11 .* D12 + Mi12 .* D22 + Mi14 .* D42) + dQdsxx .* dFdsyy .* (Mi21 .* D12 + Mi22 .* D22 + Mi24 .* D42) + dQdsxx .* dFdsxy .* (Mi31 .* D12 + Mi32 .* D22 + Mi34 .* D42) + dQdsxx .* dFdszz .* (Mi41 .* D12 + Mi42 .* D22 + Mi44 .* D42)) + D12 .* (1 - 1 ./ den .* (dQdsyy .* dFdsxx .* (Mi11 .* D12 + Mi12 .* D22 + Mi14 .* D42) + dQdsyy .* dFdsyy .* (Mi21 .* D12 + Mi22 .* D22 + Mi24 .* D42) + dQdsyy .* dFdsxy .* (Mi31 .* D12 + Mi32 .* D22 + Mi34 .* D42) + dQdsyy .* dFdszz .* (Mi41 .* D12 + Mi42 .* D22 + Mi44 .* D42))) - D14 ./ den .* (dQdszz .* dFdsxx .* (Mi11 .* D12 + Mi12 .* D22 + Mi14 .* D42) + dQdszz .* dFdsyy .* (Mi21 .* D12 + Mi22 .* D22 + Mi24 .* D42) + dQdszz .* dFdsxy .* (Mi31 .* D12 + Mi32 .* D22 + Mi34 .* D42) + dQdszz .* dFdszz .* (Mi41 .* D12 + Mi42 .* D22 + Mi44 .* D42))) + Mi22 .* (-D21 ./ den .* (dQdsxx .* dFdsxx .* (Mi11 .* D12 + Mi12 .* D22 + Mi14 .* D42) + dQdsxx .* dFdsyy .* (Mi21 .* D12 + Mi22 .* D22 + Mi24 .* D42) + dQdsxx .* dFdsxy .* (Mi31 .* D12 + Mi32 .* D22 + Mi34 .* D42) + dQdsxx .* dFdszz .* (Mi41 .* D12 + Mi42 .* D22 + Mi44 .* D42)) + D22 .* (1 - 1 ./ den .* (dQdsyy .* dFdsxx .* (Mi11 .* D12 + Mi12 .* D22 + Mi14 .* D42) + dQdsyy .* dFdsyy .* (Mi21 .* D12 + Mi22 .* D22 + Mi24 .* D42) + dQdsyy .* dFdsxy .* (Mi31 .* D12 + Mi32 .* D22 + Mi34 .* D42) + dQdsyy .* dFdszz .* (Mi41 .* D12 + Mi42 .* D22 + Mi44 .* D42))) - D24 ./ den .* (dQdszz .* dFdsxx .* (Mi11 .* D12 + Mi12 .* D22 + Mi14 .* D42) + dQdszz .* dFdsyy .* (Mi21 .* D12 + Mi22 .* D22 + Mi24 .* D42) + dQdszz .* dFdsxy .* (Mi31 .* D12 + Mi32 .* D22 + Mi34 .* D42) + dQdszz .* dFdszz .* (Mi41 .* D12 + Mi42 .* D22 + Mi44 .* D42))) - Mi23 .* D33 ./ den .* (dQdsxy .* dFdsxx .* (Mi11 .* D12 + Mi12 .* D22 + Mi14 .* D42) + dQdsxy .* dFdsyy .* (Mi21 .* D12 + Mi22 .* D22 + Mi24 .* D42) + dQdsxy .* dFdsxy .* (Mi31 .* D12 + Mi32 .* D22 + Mi34 .* D42) + dQdsxy .* dFdszz .* (Mi41 .* D12 + Mi42 .* D22 + Mi44 .* D42)) + Mi24 .* (-D41 ./ den .* (dQdsxx .* dFdsxx .* (Mi11 .* D12 + Mi12 .* D22 + Mi14 .* D42) + dQdsxx .* dFdsyy .* (Mi21 .* D12 + Mi22 .* D22 + Mi24 .* D42) + dQdsxx .* dFdsxy .* (Mi31 .* D12 + Mi32 .* D22 + Mi34 .* D42) + dQdsxx .* dFdszz .* (Mi41 .* D12 + Mi42 .* D22 + Mi44 .* D42)) + D42 .* (1 - 1 ./ den .* (dQdsyy .* dFdsxx .* (Mi11 .* D12 + Mi12 .* D22 + Mi14 .* D42) + dQdsyy .* dFdsyy .* (Mi21 .* D12 + Mi22 .* D22 + Mi24 .* D42) + dQdsyy .* dFdsxy .* (Mi31 .* D12 + Mi32 .* D22 + Mi34 .* D42) + dQdsyy .* dFdszz .* (Mi41 .* D12 + Mi42 .* D22 + Mi44 .* D42))) - D44 ./ den .* (dQdszz .* dFdsxx .* (Mi11 .* D12 + Mi12 .* D22 + Mi14 .* D42) + dQdszz .* dFdsyy .* (Mi21 .* D12 + Mi22 .* D22 + Mi24 .* D42) + dQdszz .* dFdsxy .* (Mi31 .* D12 + Mi32 .* D22 + Mi34 .* D42) + dQdszz .* dFdszz .* (Mi41 .* D12 + Mi42 .* D22 + Mi44 .* D42)));
    Dep23c = Mi21 .* (-D11 ./ den .* (dQdsxx .* dFdsxx .* Mi13 .* D33 + dQdsxx .* dFdsyy .* Mi23 .* D33 + dQdsxx .* dFdsxy .* Mi33 .* D33 + dQdsxx .* dFdszz .* Mi43 .* D33) - D12 ./ den .* (dQdsyy .* dFdsxx .* Mi13 .* D33 + dQdsyy .* dFdsyy .* Mi23 .* D33 + dQdsyy .* dFdsxy .* Mi33 .* D33 + dQdsyy .* dFdszz .* Mi43 .* D33) - D14 ./ den .* (dQdszz .* dFdsxx .* Mi13 .* D33 + dQdszz .* dFdsyy .* Mi23 .* D33 + dQdszz .* dFdsxy .* Mi33 .* D33 + dQdszz .* dFdszz .* Mi43 .* D33)) + Mi22 .* (-D21 ./ den .* (dQdsxx .* dFdsxx .* Mi13 .* D33 + dQdsxx .* dFdsyy .* Mi23 .* D33 + dQdsxx .* dFdsxy .* Mi33 .* D33 + dQdsxx .* dFdszz .* Mi43 .* D33) - D22 ./ den .* (dQdsyy .* dFdsxx .* Mi13 .* D33 + dQdsyy .* dFdsyy .* Mi23 .* D33 + dQdsyy .* dFdsxy .* Mi33 .* D33 + dQdsyy .* dFdszz .* Mi43 .* D33) - D24 ./ den .* (dQdszz .* dFdsxx .* Mi13 .* D33 + dQdszz .* dFdsyy .* Mi23 .* D33 + dQdszz .* dFdsxy .* Mi33 .* D33 + dQdszz .* dFdszz .* Mi43 .* D33)) + Mi23 .* D33 .* (1 - 1 ./ den .* (dQdsxy .* dFdsxx .* Mi13 .* D33 + dQdsxy .* dFdsyy .* Mi23 .* D33 + dQdsxy .* dFdsxy .* Mi33 .* D33 + dQdsxy .* dFdszz .* Mi43 .* D33)) + Mi24 .* (-D41 ./ den .* (dQdsxx .* dFdsxx .* Mi13 .* D33 + dQdsxx .* dFdsyy .* Mi23 .* D33 + dQdsxx .* dFdsxy .* Mi33 .* D33 + dQdsxx .* dFdszz .* Mi43 .* D33) - D42 ./ den .* (dQdsyy .* dFdsxx .* Mi13 .* D33 + dQdsyy .* dFdsyy .* Mi23 .* D33 + dQdsyy .* dFdsxy .* Mi33 .* D33 + dQdsyy .* dFdszz .* Mi43 .* D33) - D44 ./ den .* (dQdszz .* dFdsxx .* Mi13 .* D33 + dQdszz .* dFdsyy .* Mi23 .* D33 + dQdszz .* dFdsxy .* Mi33 .* D33 + dQdszz .* dFdszz .* Mi43 .* D33));

    Dep31c = Mi31 .* (D11 .* (1 - 1 ./ den .* (dQdsxx .* dFdsxx .* (Mi11 .* D11 + Mi12 .* D21 + Mi14 .* D41) + dQdsxx .* dFdsyy .* (Mi21 .* D11 + Mi22 .* D21 + Mi24 .* D41) + dQdsxx .* dFdsxy .* (Mi31 .* D11 + Mi32 .* D21 + Mi34 .* D41) + dQdsxx .* dFdszz .* (Mi41 .* D11 + Mi42 .* D21 + Mi44 .* D41))) - D12 ./ den .* (dQdsyy .* dFdsxx .* (Mi11 .* D11 + Mi12 .* D21 + Mi14 .* D41) + dQdsyy .* dFdsyy .* (Mi21 .* D11 + Mi22 .* D21 + Mi24 .* D41) + dQdsyy .* dFdsxy .* (Mi31 .* D11 + Mi32 .* D21 + Mi34 .* D41) + dQdsyy .* dFdszz .* (Mi41 .* D11 + Mi42 .* D21 + Mi44 .* D41)) - D14 ./ den .* (dQdszz .* dFdsxx .* (Mi11 .* D11 + Mi12 .* D21 + Mi14 .* D41) + dQdszz .* dFdsyy .* (Mi21 .* D11 + Mi22 .* D21 + Mi24 .* D41) + dQdszz .* dFdsxy .* (Mi31 .* D11 + Mi32 .* D21 + Mi34 .* D41) + dQdszz .* dFdszz .* (Mi41 .* D11 + Mi42 .* D21 + Mi44 .* D41))) + Mi32 .* (D21 .* (1 - 1 ./ den .* (dQdsxx .* dFdsxx .* (Mi11 .* D11 + Mi12 .* D21 + Mi14 .* D41) + dQdsxx .* dFdsyy .* (Mi21 .* D11 + Mi22 .* D21 + Mi24 .* D41) + dQdsxx .* dFdsxy .* (Mi31 .* D11 + Mi32 .* D21 + Mi34 .* D41) + dQdsxx .* dFdszz .* (Mi41 .* D11 + Mi42 .* D21 + Mi44 .* D41))) - D22 ./ den .* (dQdsyy .* dFdsxx .* (Mi11 .* D11 + Mi12 .* D21 + Mi14 .* D41) + dQdsyy .* dFdsyy .* (Mi21 .* D11 + Mi22 .* D21 + Mi24 .* D41) + dQdsyy .* dFdsxy .* (Mi31 .* D11 + Mi32 .* D21 + Mi34 .* D41) + dQdsyy .* dFdszz .* (Mi41 .* D11 + Mi42 .* D21 + Mi44 .* D41)) - D24 ./ den .* (dQdszz .* dFdsxx .* (Mi11 .* D11 + Mi12 .* D21 + Mi14 .* D41) + dQdszz .* dFdsyy .* (Mi21 .* D11 + Mi22 .* D21 + Mi24 .* D41) + dQdszz .* dFdsxy .* (Mi31 .* D11 + Mi32 .* D21 + Mi34 .* D41) + dQdszz .* dFdszz .* (Mi41 .* D11 + Mi42 .* D21 + Mi44 .* D41))) - Mi33 .* D33 ./ den .* (dQdsxy .* dFdsxx .* (Mi11 .* D11 + Mi12 .* D21 + Mi14 .* D41) + dQdsxy .* dFdsyy .* (Mi21 .* D11 + Mi22 .* D21 + Mi24 .* D41) + dQdsxy .* dFdsxy .* (Mi31 .* D11 + Mi32 .* D21 + Mi34 .* D41) + dQdsxy .* dFdszz .* (Mi41 .* D11 + Mi42 .* D21 + Mi44 .* D41)) + Mi34 .* (D41 .* (1 - 1 ./ den .* (dQdsxx .* dFdsxx .* (Mi11 .* D11 + Mi12 .* D21 + Mi14 .* D41) + dQdsxx .* dFdsyy .* (Mi21 .* D11 + Mi22 .* D21 + Mi24 .* D41) + dQdsxx .* dFdsxy .* (Mi31 .* D11 + Mi32 .* D21 + Mi34 .* D41) + dQdsxx .* dFdszz .* (Mi41 .* D11 + Mi42 .* D21 + Mi44 .* D41))) - D42 ./ den .* (dQdsyy .* dFdsxx .* (Mi11 .* D11 + Mi12 .* D21 + Mi14 .* D41) + dQdsyy .* dFdsyy .* (Mi21 .* D11 + Mi22 .* D21 + Mi24 .* D41) + dQdsyy .* dFdsxy .* (Mi31 .* D11 + Mi32 .* D21 + Mi34 .* D41) + dQdsyy .* dFdszz .* (Mi41 .* D11 + Mi42 .* D21 + Mi44 .* D41)) - D44 ./ den .* (dQdszz .* dFdsxx .* (Mi11 .* D11 + Mi12 .* D21 + Mi14 .* D41) + dQdszz .* dFdsyy .* (Mi21 .* D11 + Mi22 .* D21 + Mi24 .* D41) + dQdszz .* dFdsxy .* (Mi31 .* D11 + Mi32 .* D21 + Mi34 .* D41) + dQdszz .* dFdszz .* (Mi41 .* D11 + Mi42 .* D21 + Mi44 .* D41)));
    Dep32c = Mi31 .* (-D11 ./ den .* (dQdsxx .* dFdsxx .* (Mi11 .* D12 + Mi12 .* D22 + Mi14 .* D42) + dQdsxx .* dFdsyy .* (Mi21 .* D12 + Mi22 .* D22 + Mi24 .* D42) + dQdsxx .* dFdsxy .* (Mi31 .* D12 + Mi32 .* D22 + Mi34 .* D42) + dQdsxx .* dFdszz .* (Mi41 .* D12 + Mi42 .* D22 + Mi44 .* D42)) + D12 .* (1 - 1 ./ den .* (dQdsyy .* dFdsxx .* (Mi11 .* D12 + Mi12 .* D22 + Mi14 .* D42) + dQdsyy .* dFdsyy .* (Mi21 .* D12 + Mi22 .* D22 + Mi24 .* D42) + dQdsyy .* dFdsxy .* (Mi31 .* D12 + Mi32 .* D22 + Mi34 .* D42) + dQdsyy .* dFdszz .* (Mi41 .* D12 + Mi42 .* D22 + Mi44 .* D42))) - D14 ./ den .* (dQdszz .* dFdsxx .* (Mi11 .* D12 + Mi12 .* D22 + Mi14 .* D42) + dQdszz .* dFdsyy .* (Mi21 .* D12 + Mi22 .* D22 + Mi24 .* D42) + dQdszz .* dFdsxy .* (Mi31 .* D12 + Mi32 .* D22 + Mi34 .* D42) + dQdszz .* dFdszz .* (Mi41 .* D12 + Mi42 .* D22 + Mi44 .* D42))) + Mi32 .* (-D21 ./ den .* (dQdsxx .* dFdsxx .* (Mi11 .* D12 + Mi12 .* D22 + Mi14 .* D42) + dQdsxx .* dFdsyy .* (Mi21 .* D12 + Mi22 .* D22 + Mi24 .* D42) + dQdsxx .* dFdsxy .* (Mi31 .* D12 + Mi32 .* D22 + Mi34 .* D42) + dQdsxx .* dFdszz .* (Mi41 .* D12 + Mi42 .* D22 + Mi44 .* D42)) + D22 .* (1 - 1 ./ den .* (dQdsyy .* dFdsxx .* (Mi11 .* D12 + Mi12 .* D22 + Mi14 .* D42) + dQdsyy .* dFdsyy .* (Mi21 .* D12 + Mi22 .* D22 + Mi24 .* D42) + dQdsyy .* dFdsxy .* (Mi31 .* D12 + Mi32 .* D22 + Mi34 .* D42) + dQdsyy .* dFdszz .* (Mi41 .* D12 + Mi42 .* D22 + Mi44 .* D42))) - D24 ./ den .* (dQdszz .* dFdsxx .* (Mi11 .* D12 + Mi12 .* D22 + Mi14 .* D42) + dQdszz .* dFdsyy .* (Mi21 .* D12 + Mi22 .* D22 + Mi24 .* D42) + dQdszz .* dFdsxy .* (Mi31 .* D12 + Mi32 .* D22 + Mi34 .* D42) + dQdszz .* dFdszz .* (Mi41 .* D12 + Mi42 .* D22 + Mi44 .* D42))) - Mi33 .* D33 ./ den .* (dQdsxy .* dFdsxx .* (Mi11 .* D12 + Mi12 .* D22 + Mi14 .* D42) + dQdsxy .* dFdsyy .* (Mi21 .* D12 + Mi22 .* D22 + Mi24 .* D42) + dQdsxy .* dFdsxy .* (Mi31 .* D12 + Mi32 .* D22 + Mi34 .* D42) + dQdsxy .* dFdszz .* (Mi41 .* D12 + Mi42 .* D22 + Mi44 .* D42)) + Mi34 .* (-D41 ./ den .* (dQdsxx .* dFdsxx .* (Mi11 .* D12 + Mi12 .* D22 + Mi14 .* D42) + dQdsxx .* dFdsyy .* (Mi21 .* D12 + Mi22 .* D22 + Mi24 .* D42) + dQdsxx .* dFdsxy .* (Mi31 .* D12 + Mi32 .* D22 + Mi34 .* D42) + dQdsxx .* dFdszz .* (Mi41 .* D12 + Mi42 .* D22 + Mi44 .* D42)) + D42 .* (1 - 1 ./ den .* (dQdsyy .* dFdsxx .* (Mi11 .* D12 + Mi12 .* D22 + Mi14 .* D42) + dQdsyy .* dFdsyy .* (Mi21 .* D12 + Mi22 .* D22 + Mi24 .* D42) + dQdsyy .* dFdsxy .* (Mi31 .* D12 + Mi32 .* D22 + Mi34 .* D42) + dQdsyy .* dFdszz .* (Mi41 .* D12 + Mi42 .* D22 + Mi44 .* D42))) - D44 ./ den .* (dQdszz .* dFdsxx .* (Mi11 .* D12 + Mi12 .* D22 + Mi14 .* D42) + dQdszz .* dFdsyy .* (Mi21 .* D12 + Mi22 .* D22 + Mi24 .* D42) + dQdszz .* dFdsxy .* (Mi31 .* D12 + Mi32 .* D22 + Mi34 .* D42) + dQdszz .* dFdszz .* (Mi41 .* D12 + Mi42 .* D22 + Mi44 .* D42)));
    Dep33c = Mi31 .* (-D11 ./ den .* (dQdsxx .* dFdsxx .* Mi13 .* D33 + dQdsxx .* dFdsyy .* Mi23 .* D33 + dQdsxx .* dFdsxy .* Mi33 .* D33 + dQdsxx .* dFdszz .* Mi43 .* D33) - D12 ./ den .* (dQdsyy .* dFdsxx .* Mi13 .* D33 + dQdsyy .* dFdsyy .* Mi23 .* D33 + dQdsyy .* dFdsxy .* Mi33 .* D33 + dQdsyy .* dFdszz .* Mi43 .* D33) - D14 ./ den .* (dQdszz .* dFdsxx .* Mi13 .* D33 + dQdszz .* dFdsyy .* Mi23 .* D33 + dQdszz .* dFdsxy .* Mi33 .* D33 + dQdszz .* dFdszz .* Mi43 .* D33)) + Mi32 .* (-D21 ./ den .* (dQdsxx .* dFdsxx .* Mi13 .* D33 + dQdsxx .* dFdsyy .* Mi23 .* D33 + dQdsxx .* dFdsxy .* Mi33 .* D33 + dQdsxx .* dFdszz .* Mi43 .* D33) - D22 ./ den .* (dQdsyy .* dFdsxx .* Mi13 .* D33 + dQdsyy .* dFdsyy .* Mi23 .* D33 + dQdsyy .* dFdsxy .* Mi33 .* D33 + dQdsyy .* dFdszz .* Mi43 .* D33) - D24 ./ den .* (dQdszz .* dFdsxx .* Mi13 .* D33 + dQdszz .* dFdsyy .* Mi23 .* D33 + dQdszz .* dFdsxy .* Mi33 .* D33 + dQdszz .* dFdszz .* Mi43 .* D33)) + Mi33 .* D33 .* (1 - 1 ./ den .* (dQdsxy .* dFdsxx .* Mi13 .* D33 + dQdsxy .* dFdsyy .* Mi23 .* D33 + dQdsxy .* dFdsxy .* Mi33 .* D33 + dQdsxy .* dFdszz .* Mi43 .* D33)) + Mi34 .* (-D41 ./ den .* (dQdsxx .* dFdsxx .* Mi13 .* D33 + dQdsxx .* dFdsyy .* Mi23 .* D33 + dQdsxx .* dFdsxy .* Mi33 .* D33 + dQdsxx .* dFdszz .* Mi43 .* D33) - D42 ./ den .* (dQdsyy .* dFdsxx .* Mi13 .* D33 + dQdsyy .* dFdsyy .* Mi23 .* D33 + dQdsyy .* dFdsxy .* Mi33 .* D33 + dQdsyy .* dFdszz .* Mi43 .* D33) - D44 ./ den .* (dQdszz .* dFdsxx .* Mi13 .* D33 + dQdszz .* dFdsyy .* Mi23 .* D33 + dQdszz .* dFdsxy .* Mi33 .* D33 + dQdszz .* dFdszz .* Mi43 .* D33));

    [Dep11c Dep12c Dep13c ;
     Dep21c Dep22c Dep23c ;
     Dep31c Dep32c Dep33c ;   
    ]

%%
clear all; 
load test.mat
                [~,B,en] = get_data ( x_l , node , di , form ) ;
                dEps_t = B*dU(en); dEps_t = [dEps_t;0]; % add z strain
%                 dEps_t(3) = dEps_t(3)/2 ;
                dEps_0 = dEps_t; % save total strain  
                
                [Dmat_t,Dmat,Gve,K_] = identify_tangent_v2 ( x_l , mat) ;
%                 Dmat_t(3,3) = Dmat_t(3,3)*2;
                dSigma = Dmat_t * [dEps_t];% total stress increment 
                dP    =-1/3*(sum(dSigma([1 2 4]))); % pressure increment
                
                Sigma_t = VE1* TSigma0 - P0 + dSigma; % compute total trial stress (eq. 1)
                
                P   = P0 + dP; % get pressure
                TSigma = [P;P;0;P] + Sigma_t; % get deviatoric components
                J2 = 1/2*(TSigma(1)^2+TSigma(2)^2+TSigma(4)^2)+TSigma(3)^2;
                C     = C0; 
                F     =  sqrt(J2) - C.*cos_phi - P.*sin_phi; % trial yield function
        fprintf('max. trial     F = %2.2e\n', max(F(:)))
                % Plastic corrections
%                 if F>0 % if trial stresses exceed yield 
                    dQdsigma = get_dQdSigma ( TSigma , sin_psi , J2 );
                    h1 = cos_phi.*h.*sqrt(2/3).*sqrt((dQdsigma(1).^2 + 2*dQdsigma(3).^2 + dQdsigma(2).^2 + dQdsigma(4).^2));
                    dg = (F./(Gve + K_.*sin_phi*sin_psi + eta_vp/dt + h1 ));% increment of plastic multiplier
                    dEps_t = [dEps_0(1); dEps_0(2) ; dEps_0(3); dEps_0(4)]  - ...
                        dg*[dQdsigma(1);dQdsigma(2);dQdsigma(3)/1;dQdsigma(4)]; % compute corrected strain
                    dSigma = Dmat_t * dEps_t ;% total stress increment 
                    dP    =-1/3*(sum(dSigma([1 2 4]))); % pressure increment
                    Sigma_t = VE1 .* TSigma0 - P0 + dSigma; % compute total stresses
                    P   = P0 + dP; % get pressure
                    TSigma = [P;P;0;P] + Sigma_t; % get deviatoric components
                
                    % Check yield function (F_vp=0)
                    J2 = 1/2*(TSigma(1).^2+TSigma(2).^2+TSigma(4).^2)+TSigma(3).^2;
                    dep=sqrt(2/3).*sqrt((dg.*dQdsigma(1)).^2+(dg.*dQdsigma(2)).^2+(dg.*dQdsigma(4)).^2+2*(dg.*dQdsigma(3)).^2); % eq. 

                    C     = C0 + h.*dep; 
                    F     =  sqrt(J2) - (C.*cos_phi + P.*sin_phi) ;                       % Backbone plastic yield function
                    F_vp   =  sqrt(J2) - (C.*cos_phi + P.*sin_phi) - (dg./dt).*eta_vp;  % Consistency viscoplasticity model
%                     sqrt(J2)
%                     (C.*cos_phi + P.*sin_phi)
%                     (dg./dt).*eta_vp
            fprintf('    Backbone EP max. Fc = %2.2e \n', max(F(:)))
            fprintf('    Viscoplast. max. Fc = %2.2e \n', max(F_vp(:)))
%                 end % F > 0
                    

                    dQdSigma = get_dQdSigma ( TSigma , sin_psi , J2 )';
                    dFdSigma = get_dQdSigma ( TSigma , sin_phi , J2 )';
                    d2QdSigma2 = get_dQ2dSigma2 ( TSigma  , J2 );
                    [Dmat_t,Dmat,Gve,K_] = identify_tangent_v2 ( x_l , mat) ;

                    I = [ 1 0 0 0 ; 0 1 0 0 ; 0 0 1 0 ; 0 0 0 1] ;    
                    M = I + dg * Dmat_t * d2QdSigma2 ;
                    Mi = inv(M) ; 
                    Dvp = Mi * Dmat_t -  ( Mi * Dmat_t * dQdSigma * dFdSigma' * Mi * Dmat_t ) / (  eta_vp/dt + h1 + dFdSigma' * Mi * Dmat_t * dQdSigma ) ; 
Dvp(1:3,1:3)
% M