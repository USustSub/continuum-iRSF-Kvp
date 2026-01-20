
K = sparse(2*numnode,2*numnode);
F = zeros(2*numnode,1) ; 
Residual = zeros(2*numnode,1) ; 

if strcmp(stype,'elastic')
%%
% loop over nodes and assemble stiffness matrix
for ij = 1 : size(node,1)
    if rem(ij,10) == 0 
        ij/size(node,1)
    end
    node_c = node(ij,:) ; 
    
    
    [Dmat_t,Dmat,G,K_] = identify_tangent_v2 ( node_c , mat) ;

% only for nodes inside domain
if (node_c(1) > 0 && node_c(1) < L && node_c(2) < D && node_c(2) > 0 )
    x_l = node_c - [deltaX/2 0 ] ;
    x_r = node_c + [deltaX/2 0 ] ;    
    
    x_b = node_c - [0 deltaY/2 ] ;
    x_t = node_c + [0 deltaY/2 ] ;    

    get_equations ();
end

    get_BCs () ;

end

elseif strcmp(stype,'viscoplastic')
%%

Dep11v    = Kv + 4/3.*Gv;  Dep12v    = Kv - 2/3.*Gv;  Dep13v    = 0*Kv;  Dep14v    = Kv - 2/3.*Gv;
Dep21v    = Kv - 2/3.*Gv;  Dep22v    = Kv + 4/3.*Gv;  Dep23v    = 0*Kv;  Dep24v    = Kv - 2/3.*Gv;
Dep31v    = 0.*Gv;         Dep32v    = 0.*Gv;         Dep33v    = 1.*Gv; Dep34v    = 0.*Gv;
Dep41v    = Kv - 2/3.*Gv;  Dep42v    = Kv - 2/3.*Gv;  Dep43v    = 0*Kv;  Dep44v    = Kv + 4/3.*Gv;


dlam = dlamc;
Pc   = -1/3*(Sxxc_t+Syyc_t+Szzc_t); 
Txxc = Pc + Sxxc_t; 
Tyyc = Pc + Syyc_t; 
Tzzc = Pc + Szzc_t; 
Txyc = Txyc_t;
J2c  = 1/2*(Txxc.^2 + Tyyc.^2 + Tzzc.^2) + Txyc.^2;
txx = Pc+Sxxc_t;
tyy = Pc+Syyc_t;
txy = Txyc_t;
tzz = Pc+Szzc_t; J2 = J2c;
% dQds centers
dQdsxx = vm * txx .* (vm * J2) .^ (-0.1e1 / 0.2e1) / 0.2e1 + sin_psi / 0.3e1;
dQdsyy = vm * tyy .* (vm * J2) .^ (-0.1e1 / 0.2e1) / 0.2e1 + sin_psi / 0.3e1;
dQdsxy = 1 * vm * txy .* (vm * J2) .^ (-0.1e1 / 0.2e1);
dQdszz = vm * tzz .* (vm * J2) .^ (-0.1e1 / 0.2e1) / 0.2e1 + sin_psi / 0.3e1;
% dFds centers
dFdsxx = vm * txx .* (vm * J2) .^ (-0.1e1 / 0.2e1) / 0.2e1 + sin_phi / 0.3e1;
dFdsyy = vm * tyy .* (vm * J2) .^ (-0.1e1 / 0.2e1) / 0.2e1 + sin_phi / 0.3e1;
dFdsxy = 0.1000000000e1 * vm * txy .* (vm * J2) .^ (-0.1e1 / 0.2e1);
dFdszz = vm * tzz .* (vm * J2) .^ (-0.1e1 / 0.2e1) / 0.2e1 + sin_phi / 0.3e1;
% Reduced elastic stiffness matrix
if consistent==1
    d2Qdsxxdsxx = vm .* (vm .* J2) .^ (-0.1e1 ./ 0.2e1) ./ 0.3e1 - vm .^ 2 .* txx .^ 2 .* (vm .* J2) .^ (-0.3e1 ./ 0.2e1) ./ 0.4e1;
    d2Qdsxxdsyy = -vm .* (vm .* J2) .^ (-0.1e1 ./ 0.2e1) ./ 0.6e1 - vm .^ 2 .* txx .* tyy .* (vm .* J2) .^ (-0.3e1 ./ 0.2e1) ./ 0.4e1;
    d2Qdsxxdsxy = -vm .^ 2 .* txx .* txy .* (vm .* J2) .^ (-0.3e1 ./ 0.2e1) ./ 0.2e1;
    d2Qdsxxdszz = -vm .* (vm .* J2) .^ (-0.1e1 ./ 0.2e1) ./ 0.6e1 - vm .^ 2 .* txx .* tzz .* (vm .* J2) .^ (-0.3e1 ./ 0.2e1) ./ 0.4e1;
    d2Qdsyydsxx = -vm .* (vm .* J2) .^ (-0.1e1 ./ 0.2e1) ./ 0.6e1 - vm .^ 2 .* txx .* tyy .* (vm .* J2) .^ (-0.3e1 ./ 0.2e1) ./ 0.4e1;
    d2Qdsyydsyy = vm .* (vm .* J2) .^ (-0.1e1 ./ 0.2e1) ./ 0.3e1 - vm .^ 2 .* tyy .^ 2 .* (vm .* J2) .^ (-0.3e1 ./ 0.2e1) ./ 0.4e1;
    d2Qdsyydsxy = -vm .^ 2 .* tyy .* txy .* (vm .* J2) .^ (-0.3e1 ./ 0.2e1) ./ 0.2e1;
    d2Qdsyydszz = -vm .* (vm .* J2) .^ (-0.1e1 ./ 0.2e1) ./ 0.6e1 - vm .^ 2 .* tyy .* tzz .* (vm .* J2) .^ (-0.3e1 ./ 0.2e1) ./ 0.4e1;
    d2Qdsxydsxx = -vm .^ 2 .* txx .* txy .* (vm .* J2) .^ (-0.3e1 ./ 0.2e1) ./ 0.2e1;
    d2Qdsxydsyy = -vm .^ 2 .* tyy .* txy .* (vm .* J2) .^ (-0.3e1 ./ 0.2e1) ./ 0.2e1;
    d2Qdsxydsxy = vm .* (vm .* J2) .^ (-0.1e1 ./ 0.2e1) - vm .^ 2 .* txy .^ 2 .* (vm .* J2) .^ (-0.3e1 ./ 0.2e1);
    d2Qdsxydszz = -vm .^ 2 .* txy .* tzz .* (vm .* J2) .^ (-0.3e1 ./ 0.2e1) ./ 0.2e1;
    d2Qdszzdsxx = -vm .* (vm .* J2) .^ (-0.1e1 ./ 0.2e1) ./ 0.6e1 - vm .^ 2 .* txx .* tzz .* (vm .* J2) .^ (-0.3e1 ./ 0.2e1) ./ 0.4e1;
    d2Qdszzdsyy = -vm .* (vm .* J2) .^ (-0.1e1 ./ 0.2e1) ./ 0.6e1 - vm .^ 2 .* tyy .* tzz .* (vm .* J2) .^ (-0.3e1 ./ 0.2e1) ./ 0.4e1;
    d2Qdszzdsxy = -vm .^ 2 .* txy .* tzz .* (vm .* J2) .^ (-0.3e1 ./ 0.2e1) ./ 0.2e1;
    d2Qdszzdszz = vm .* (vm .* J2) .^ (-0.1e1 ./ 0.2e1) ./ 0.3e1 - vm .^ 2 .* tzz .^ 2 .* (vm .* J2) .^ (-0.3e1 ./ 0.2e1) ./ 0.4e1;
    
    syms d2Qdsxxdsxx d2Qdsxxdsyy d2Qdsxxdsxy d2Qdsxxdszz
    syms d2Qdsyydsxx d2Qdsyydsyy d2Qdsyydsxy d2Qdsyydszz
    syms d2Qdsxydsxx d2Qdsxydsyy d2Qdsxydsxy d2Qdsxydszz
    syms d2Qdszzdsxx d2Qdszzdsyy d2Qdszzdsxy d2Qdszzdszz
    d2QdSigma2 =  [
                    d2Qdsxxdsxx d2Qdsxxdsyy d2Qdsxxdsxy d2Qdsxxdszz  
                    d2Qdsyydsxx d2Qdsyydsyy d2Qdsyydsxy d2Qdsyydszz  
                    d2Qdsxydsxx d2Qdsxydsyy d2Qdsxydsxy d2Qdsxydszz  
                    d2Qdszzdsxx d2Qdszzdsyy d2Qdszzdsxy d2Qdszzdszz  
                    ]; 
	syms dQdsxx dQdsyy dQdsxy dQdszz 
    assume(dQdsxx>0)     
    assume(dQdsyy>0)
    assume(dQdsxy>0)     
    assume(dQdszz>0)
    dQdSigma = [dQdsxx dQdsyy dQdsxy dQdszz]' ; 
	syms dFdsxx dFdsyy dFdsxy dFdszz 
    assume(dFdsxx>0)     
    assume(dFdsyy>0)
    assume(dFdsxy>0)     
    assume(dFdszz>0)
    dFdSigma = [dFdsxx dFdsyy dFdsxy dFdszz]' ; 

    syms D11 D12 D13 D14 D21 D22 D23 D24 D31 D32 D33 D34 D41 D42 D43 D44
    syms dg eta_vp dt
    Dve = [D11 D12 0 D14 ; 
           D21 D22 0 D24 ;
           0 0 D33 0 ;
           D41 D42 0 D44 ; 
           ]; 
   
    I = [ 1 0 0 0 ; 0 1 0 0 ; 0 0 1 0 ; 0 0 0 1] ;    
    M = I + dg * Dve * d2QdSigma2 ;
%     Mi = inv(M) ; 
    syms Mi11 Mi12 Mi13 Mi14 Mi21 Mi22 Mi23 Mi24 Mi31 Mi32 Mi33 Mi34 Mi41 Mi42 Mi43 Mi44
    Mi = [ Mi11 Mi12 Mi13 Mi14 ; 
           Mi21 Mi22 Mi23 Mi24 ;
           Mi31 Mi32 Mi33 Mi34 ;
           Mi41 Mi42 Mi43 Mi44	]; 

    Dvp = Mi * Dve -  ( Mi * Dve * dQdSigma * dFdSigma' * Mi * Dve ) / (  eta_vp/dt + h1c + dFdSigma' * Mi * Dve * dQdSigma ) ; 
%     
%     D11 = De11c; D12 = De12c; D13 = De13c; D14 = De14c;
%     D21 = De21c; D22 = De22c; D23 = De23c; D24 = De24c;
%     D31 = De31c; D32 = De32c; D33 = De33c; D34 = De34c;
%     D41 = De41c; D42 = De42c; D43 = De43c; D44 = De44c;
%     M11 = 1 + dlam .* (D11 .* d2Qdsxxdsxx + D12 .* d2Qdsyydsxx + D14 .* d2Qdszzdsxx);
%     M12 = dlam .* (D11 .* d2Qdsxxdsyy + D12 .* d2Qdsyydsyy + D14 .* d2Qdszzdsyy);
%     M13 = dlam .* (D11 .* d2Qdsxxdsxy + D12 .* d2Qdsyydsxy + D14 .* d2Qdszzdsxy);
%     M14 = dlam .* (D11 .* d2Qdsxxdszz + D12 .* d2Qdsyydszz + D14 .* d2Qdszzdszz);
%     M21 = dlam .* (D21 .* d2Qdsxxdsxx + D22 .* d2Qdsyydsxx + D24 .* d2Qdszzdsxx);
%     M22 = 1 + dlam .* (D21 .* d2Qdsxxdsyy + D22 .* d2Qdsyydsyy + D24 .* d2Qdszzdsyy);
%     M23 = dlam .* (D21 .* d2Qdsxxdsxy + D22 .* d2Qdsyydsxy + D24 .* d2Qdszzdsxy);
%     M24 = dlam .* (D21 .* d2Qdsxxdszz + D22 .* d2Qdsyydszz + D24 .* d2Qdszzdszz);
%     M31 = dlam .* D33 .* d2Qdsxydsxx;
%     M32 = dlam .* D33 .* d2Qdsxydsyy;
%     M33 = 1 + dlam .* D33 .* d2Qdsxydsxy;
%     M34 = dlam .* D33 .* d2Qdsxydszz;
%     M41 = dlam .* (D41 .* d2Qdsxxdsxx + D42 .* d2Qdsyydsxx + D44 .* d2Qdszzdsxx);
%     M42 = dlam .* (D41 .* d2Qdsxxdsyy + D42 .* d2Qdsyydsyy + D44 .* d2Qdszzdsyy);
%     M43 = dlam .* (D41 .* d2Qdsxxdsxy + D42 .* d2Qdsyydsxy + D44 .* d2Qdszzdsxy);
%     M44 = 1 + dlam .* (D41 .* d2Qdsxxdszz + D42 .* d2Qdsyydszz + D44 .* d2Qdszzdszz);
%     [ Mi11,Mi12,Mi13,Mi14, Mi21,Mi22,Mi23,Mi24, Mi31,Mi32,Mi33,Mi34, Mi41,Mi42,Mi43,Mi44  ] = M2Di_EP5_inverse_4x4( M11,M12,M13,M14, M21,M22,M23,M24, M31,M32,M33,M34, M41,M42,M43,M44  );
%     den = (dFdsxx .* Mi11 + dFdsyy .* Mi21 + dFdsxy .* Mi31 + dFdszz .* Mi41) .* (D11 .* dQdsxx + D12 .* dQdsyy + D14 .* dQdszz) + (dFdsxx .* Mi12 + dFdsyy .* Mi22 + dFdsxy .* Mi32 + dFdszz .* Mi42) .* (D21 .* dQdsxx + D22 .* dQdsyy + D24 .* dQdszz) + (dFdsxx .* Mi13 + dFdsyy .* Mi23 + dFdsxy .* Mi33 + dFdszz .* Mi43) .* D33 .* dQdsxy + (dFdsxx .* Mi14 + dFdsyy .* Mi24 + dFdsxy .* Mi34 + dFdszz .* Mi44) .* (D41 .* dQdsxx + D42 .* dQdsyy + D44 .* dQdszz);;
% 
%     den = den + vp.*eta_vp/dt + h1c;
%     Dep11c = Mi11 .* (D11 .* (1 - 1 ./ den .* (dQdsxx .* dFdsxx .* (Mi11 .* D11 + Mi12 .* D21 + Mi14 .* D41) + dQdsxx .* dFdsyy .* (Mi21 .* D11 + Mi22 .* D21 + Mi24 .* D41) + dQdsxx .* dFdsxy .* (Mi31 .* D11 + Mi32 .* D21 + Mi34 .* D41) + dQdsxx .* dFdszz .* (Mi41 .* D11 + Mi42 .* D21 + Mi44 .* D41))) - D12 ./ den .* (dQdsyy .* dFdsxx .* (Mi11 .* D11 + Mi12 .* D21 + Mi14 .* D41) + dQdsyy .* dFdsyy .* (Mi21 .* D11 + Mi22 .* D21 + Mi24 .* D41) + dQdsyy .* dFdsxy .* (Mi31 .* D11 + Mi32 .* D21 + Mi34 .* D41) + dQdsyy .* dFdszz .* (Mi41 .* D11 + Mi42 .* D21 + Mi44 .* D41)) - D14 ./ den .* (dQdszz .* dFdsxx .* (Mi11 .* D11 + Mi12 .* D21 + Mi14 .* D41) + dQdszz .* dFdsyy .* (Mi21 .* D11 + Mi22 .* D21 + Mi24 .* D41) + dQdszz .* dFdsxy .* (Mi31 .* D11 + Mi32 .* D21 + Mi34 .* D41) + dQdszz .* dFdszz .* (Mi41 .* D11 + Mi42 .* D21 + Mi44 .* D41))) + Mi12 .* (D21 .* (1 - 1 ./ den .* (dQdsxx .* dFdsxx .* (Mi11 .* D11 + Mi12 .* D21 + Mi14 .* D41) + dQdsxx .* dFdsyy .* (Mi21 .* D11 + Mi22 .* D21 + Mi24 .* D41) + dQdsxx .* dFdsxy .* (Mi31 .* D11 + Mi32 .* D21 + Mi34 .* D41) + dQdsxx .* dFdszz .* (Mi41 .* D11 + Mi42 .* D21 + Mi44 .* D41))) - D22 ./ den .* (dQdsyy .* dFdsxx .* (Mi11 .* D11 + Mi12 .* D21 + Mi14 .* D41) + dQdsyy .* dFdsyy .* (Mi21 .* D11 + Mi22 .* D21 + Mi24 .* D41) + dQdsyy .* dFdsxy .* (Mi31 .* D11 + Mi32 .* D21 + Mi34 .* D41) + dQdsyy .* dFdszz .* (Mi41 .* D11 + Mi42 .* D21 + Mi44 .* D41)) - D24 ./ den .* (dQdszz .* dFdsxx .* (Mi11 .* D11 + Mi12 .* D21 + Mi14 .* D41) + dQdszz .* dFdsyy .* (Mi21 .* D11 + Mi22 .* D21 + Mi24 .* D41) + dQdszz .* dFdsxy .* (Mi31 .* D11 + Mi32 .* D21 + Mi34 .* D41) + dQdszz .* dFdszz .* (Mi41 .* D11 + Mi42 .* D21 + Mi44 .* D41))) - Mi13 .* D33 ./ den .* (dQdsxy .* dFdsxx .* (Mi11 .* D11 + Mi12 .* D21 + Mi14 .* D41) + dQdsxy .* dFdsyy .* (Mi21 .* D11 + Mi22 .* D21 + Mi24 .* D41) + dQdsxy .* dFdsxy .* (Mi31 .* D11 + Mi32 .* D21 + Mi34 .* D41) + dQdsxy .* dFdszz .* (Mi41 .* D11 + Mi42 .* D21 + Mi44 .* D41)) + Mi14 .* (D41 .* (1 - 1 ./ den .* (dQdsxx .* dFdsxx .* (Mi11 .* D11 + Mi12 .* D21 + Mi14 .* D41) + dQdsxx .* dFdsyy .* (Mi21 .* D11 + Mi22 .* D21 + Mi24 .* D41) + dQdsxx .* dFdsxy .* (Mi31 .* D11 + Mi32 .* D21 + Mi34 .* D41) + dQdsxx .* dFdszz .* (Mi41 .* D11 + Mi42 .* D21 + Mi44 .* D41))) - D42 ./ den .* (dQdsyy .* dFdsxx .* (Mi11 .* D11 + Mi12 .* D21 + Mi14 .* D41) + dQdsyy .* dFdsyy .* (Mi21 .* D11 + Mi22 .* D21 + Mi24 .* D41) + dQdsyy .* dFdsxy .* (Mi31 .* D11 + Mi32 .* D21 + Mi34 .* D41) + dQdsyy .* dFdszz .* (Mi41 .* D11 + Mi42 .* D21 + Mi44 .* D41)) - D44 ./ den .* (dQdszz .* dFdsxx .* (Mi11 .* D11 + Mi12 .* D21 + Mi14 .* D41) + dQdszz .* dFdsyy .* (Mi21 .* D11 + Mi22 .* D21 + Mi24 .* D41) + dQdszz .* dFdsxy .* (Mi31 .* D11 + Mi32 .* D21 + Mi34 .* D41) + dQdszz .* dFdszz .* (Mi41 .* D11 + Mi42 .* D21 + Mi44 .* D41)));
%     Dep12c = Mi11 .* (-D11 ./ den .* (dQdsxx .* dFdsxx .* (Mi11 .* D12 + Mi12 .* D22 + Mi14 .* D42) + dQdsxx .* dFdsyy .* (Mi21 .* D12 + Mi22 .* D22 + Mi24 .* D42) + dQdsxx .* dFdsxy .* (Mi31 .* D12 + Mi32 .* D22 + Mi34 .* D42) + dQdsxx .* dFdszz .* (Mi41 .* D12 + Mi42 .* D22 + Mi44 .* D42)) + D12 .* (1 - 1 ./ den .* (dQdsyy .* dFdsxx .* (Mi11 .* D12 + Mi12 .* D22 + Mi14 .* D42) + dQdsyy .* dFdsyy .* (Mi21 .* D12 + Mi22 .* D22 + Mi24 .* D42) + dQdsyy .* dFdsxy .* (Mi31 .* D12 + Mi32 .* D22 + Mi34 .* D42) + dQdsyy .* dFdszz .* (Mi41 .* D12 + Mi42 .* D22 + Mi44 .* D42))) - D14 ./ den .* (dQdszz .* dFdsxx .* (Mi11 .* D12 + Mi12 .* D22 + Mi14 .* D42) + dQdszz .* dFdsyy .* (Mi21 .* D12 + Mi22 .* D22 + Mi24 .* D42) + dQdszz .* dFdsxy .* (Mi31 .* D12 + Mi32 .* D22 + Mi34 .* D42) + dQdszz .* dFdszz .* (Mi41 .* D12 + Mi42 .* D22 + Mi44 .* D42))) + Mi12 .* (-D21 ./ den .* (dQdsxx .* dFdsxx .* (Mi11 .* D12 + Mi12 .* D22 + Mi14 .* D42) + dQdsxx .* dFdsyy .* (Mi21 .* D12 + Mi22 .* D22 + Mi24 .* D42) + dQdsxx .* dFdsxy .* (Mi31 .* D12 + Mi32 .* D22 + Mi34 .* D42) + dQdsxx .* dFdszz .* (Mi41 .* D12 + Mi42 .* D22 + Mi44 .* D42)) + D22 .* (1 - 1 ./ den .* (dQdsyy .* dFdsxx .* (Mi11 .* D12 + Mi12 .* D22 + Mi14 .* D42) + dQdsyy .* dFdsyy .* (Mi21 .* D12 + Mi22 .* D22 + Mi24 .* D42) + dQdsyy .* dFdsxy .* (Mi31 .* D12 + Mi32 .* D22 + Mi34 .* D42) + dQdsyy .* dFdszz .* (Mi41 .* D12 + Mi42 .* D22 + Mi44 .* D42))) - D24 ./ den .* (dQdszz .* dFdsxx .* (Mi11 .* D12 + Mi12 .* D22 + Mi14 .* D42) + dQdszz .* dFdsyy .* (Mi21 .* D12 + Mi22 .* D22 + Mi24 .* D42) + dQdszz .* dFdsxy .* (Mi31 .* D12 + Mi32 .* D22 + Mi34 .* D42) + dQdszz .* dFdszz .* (Mi41 .* D12 + Mi42 .* D22 + Mi44 .* D42))) - Mi13 .* D33 ./ den .* (dQdsxy .* dFdsxx .* (Mi11 .* D12 + Mi12 .* D22 + Mi14 .* D42) + dQdsxy .* dFdsyy .* (Mi21 .* D12 + Mi22 .* D22 + Mi24 .* D42) + dQdsxy .* dFdsxy .* (Mi31 .* D12 + Mi32 .* D22 + Mi34 .* D42) + dQdsxy .* dFdszz .* (Mi41 .* D12 + Mi42 .* D22 + Mi44 .* D42)) + Mi14 .* (-D41 ./ den .* (dQdsxx .* dFdsxx .* (Mi11 .* D12 + Mi12 .* D22 + Mi14 .* D42) + dQdsxx .* dFdsyy .* (Mi21 .* D12 + Mi22 .* D22 + Mi24 .* D42) + dQdsxx .* dFdsxy .* (Mi31 .* D12 + Mi32 .* D22 + Mi34 .* D42) + dQdsxx .* dFdszz .* (Mi41 .* D12 + Mi42 .* D22 + Mi44 .* D42)) + D42 .* (1 - 1 ./ den .* (dQdsyy .* dFdsxx .* (Mi11 .* D12 + Mi12 .* D22 + Mi14 .* D42) + dQdsyy .* dFdsyy .* (Mi21 .* D12 + Mi22 .* D22 + Mi24 .* D42) + dQdsyy .* dFdsxy .* (Mi31 .* D12 + Mi32 .* D22 + Mi34 .* D42) + dQdsyy .* dFdszz .* (Mi41 .* D12 + Mi42 .* D22 + Mi44 .* D42))) - D44 ./ den .* (dQdszz .* dFdsxx .* (Mi11 .* D12 + Mi12 .* D22 + Mi14 .* D42) + dQdszz .* dFdsyy .* (Mi21 .* D12 + Mi22 .* D22 + Mi24 .* D42) + dQdszz .* dFdsxy .* (Mi31 .* D12 + Mi32 .* D22 + Mi34 .* D42) + dQdszz .* dFdszz .* (Mi41 .* D12 + Mi42 .* D22 + Mi44 .* D42)));
%     Dep13c = Mi11 .* (-D11 ./ den .* (dQdsxx .* dFdsxx .* Mi13 .* D33 + dQdsxx .* dFdsyy .* Mi23 .* D33 + dQdsxx .* dFdsxy .* Mi33 .* D33 + dQdsxx .* dFdszz .* Mi43 .* D33) - D12 ./ den .* (dQdsyy .* dFdsxx .* Mi13 .* D33 + dQdsyy .* dFdsyy .* Mi23 .* D33 + dQdsyy .* dFdsxy .* Mi33 .* D33 + dQdsyy .* dFdszz .* Mi43 .* D33) - D14 ./ den .* (dQdszz .* dFdsxx .* Mi13 .* D33 + dQdszz .* dFdsyy .* Mi23 .* D33 + dQdszz .* dFdsxy .* Mi33 .* D33 + dQdszz .* dFdszz .* Mi43 .* D33)) + Mi12 .* (-D21 ./ den .* (dQdsxx .* dFdsxx .* Mi13 .* D33 + dQdsxx .* dFdsyy .* Mi23 .* D33 + dQdsxx .* dFdsxy .* Mi33 .* D33 + dQdsxx .* dFdszz .* Mi43 .* D33) - D22 ./ den .* (dQdsyy .* dFdsxx .* Mi13 .* D33 + dQdsyy .* dFdsyy .* Mi23 .* D33 + dQdsyy .* dFdsxy .* Mi33 .* D33 + dQdsyy .* dFdszz .* Mi43 .* D33) - D24 ./ den .* (dQdszz .* dFdsxx .* Mi13 .* D33 + dQdszz .* dFdsyy .* Mi23 .* D33 + dQdszz .* dFdsxy .* Mi33 .* D33 + dQdszz .* dFdszz .* Mi43 .* D33)) + Mi13 .* D33 .* (1 - 1 ./ den .* (dQdsxy .* dFdsxx .* Mi13 .* D33 + dQdsxy .* dFdsyy .* Mi23 .* D33 + dQdsxy .* dFdsxy .* Mi33 .* D33 + dQdsxy .* dFdszz .* Mi43 .* D33)) + Mi14 .* (-D41 ./ den .* (dQdsxx .* dFdsxx .* Mi13 .* D33 + dQdsxx .* dFdsyy .* Mi23 .* D33 + dQdsxx .* dFdsxy .* Mi33 .* D33 + dQdsxx .* dFdszz .* Mi43 .* D33) - D42 ./ den .* (dQdsyy .* dFdsxx .* Mi13 .* D33 + dQdsyy .* dFdsyy .* Mi23 .* D33 + dQdsyy .* dFdsxy .* Mi33 .* D33 + dQdsyy .* dFdszz .* Mi43 .* D33) - D44 ./ den .* (dQdszz .* dFdsxx .* Mi13 .* D33 + dQdszz .* dFdsyy .* Mi23 .* D33 + dQdszz .* dFdsxy .* Mi33 .* D33 + dQdszz .* dFdszz .* Mi43 .* D33));
%     Dep21c = Mi21 .* (D11 .* (1 - 1 ./ den .* (dQdsxx .* dFdsxx .* (Mi11 .* D11 + Mi12 .* D21 + Mi14 .* D41) + dQdsxx .* dFdsyy .* (Mi21 .* D11 + Mi22 .* D21 + Mi24 .* D41) + dQdsxx .* dFdsxy .* (Mi31 .* D11 + Mi32 .* D21 + Mi34 .* D41) + dQdsxx .* dFdszz .* (Mi41 .* D11 + Mi42 .* D21 + Mi44 .* D41))) - D12 ./ den .* (dQdsyy .* dFdsxx .* (Mi11 .* D11 + Mi12 .* D21 + Mi14 .* D41) + dQdsyy .* dFdsyy .* (Mi21 .* D11 + Mi22 .* D21 + Mi24 .* D41) + dQdsyy .* dFdsxy .* (Mi31 .* D11 + Mi32 .* D21 + Mi34 .* D41) + dQdsyy .* dFdszz .* (Mi41 .* D11 + Mi42 .* D21 + Mi44 .* D41)) - D14 ./ den .* (dQdszz .* dFdsxx .* (Mi11 .* D11 + Mi12 .* D21 + Mi14 .* D41) + dQdszz .* dFdsyy .* (Mi21 .* D11 + Mi22 .* D21 + Mi24 .* D41) + dQdszz .* dFdsxy .* (Mi31 .* D11 + Mi32 .* D21 + Mi34 .* D41) + dQdszz .* dFdszz .* (Mi41 .* D11 + Mi42 .* D21 + Mi44 .* D41))) + Mi22 .* (D21 .* (1 - 1 ./ den .* (dQdsxx .* dFdsxx .* (Mi11 .* D11 + Mi12 .* D21 + Mi14 .* D41) + dQdsxx .* dFdsyy .* (Mi21 .* D11 + Mi22 .* D21 + Mi24 .* D41) + dQdsxx .* dFdsxy .* (Mi31 .* D11 + Mi32 .* D21 + Mi34 .* D41) + dQdsxx .* dFdszz .* (Mi41 .* D11 + Mi42 .* D21 + Mi44 .* D41))) - D22 ./ den .* (dQdsyy .* dFdsxx .* (Mi11 .* D11 + Mi12 .* D21 + Mi14 .* D41) + dQdsyy .* dFdsyy .* (Mi21 .* D11 + Mi22 .* D21 + Mi24 .* D41) + dQdsyy .* dFdsxy .* (Mi31 .* D11 + Mi32 .* D21 + Mi34 .* D41) + dQdsyy .* dFdszz .* (Mi41 .* D11 + Mi42 .* D21 + Mi44 .* D41)) - D24 ./ den .* (dQdszz .* dFdsxx .* (Mi11 .* D11 + Mi12 .* D21 + Mi14 .* D41) + dQdszz .* dFdsyy .* (Mi21 .* D11 + Mi22 .* D21 + Mi24 .* D41) + dQdszz .* dFdsxy .* (Mi31 .* D11 + Mi32 .* D21 + Mi34 .* D41) + dQdszz .* dFdszz .* (Mi41 .* D11 + Mi42 .* D21 + Mi44 .* D41))) - Mi23 .* D33 ./ den .* (dQdsxy .* dFdsxx .* (Mi11 .* D11 + Mi12 .* D21 + Mi14 .* D41) + dQdsxy .* dFdsyy .* (Mi21 .* D11 + Mi22 .* D21 + Mi24 .* D41) + dQdsxy .* dFdsxy .* (Mi31 .* D11 + Mi32 .* D21 + Mi34 .* D41) + dQdsxy .* dFdszz .* (Mi41 .* D11 + Mi42 .* D21 + Mi44 .* D41)) + Mi24 .* (D41 .* (1 - 1 ./ den .* (dQdsxx .* dFdsxx .* (Mi11 .* D11 + Mi12 .* D21 + Mi14 .* D41) + dQdsxx .* dFdsyy .* (Mi21 .* D11 + Mi22 .* D21 + Mi24 .* D41) + dQdsxx .* dFdsxy .* (Mi31 .* D11 + Mi32 .* D21 + Mi34 .* D41) + dQdsxx .* dFdszz .* (Mi41 .* D11 + Mi42 .* D21 + Mi44 .* D41))) - D42 ./ den .* (dQdsyy .* dFdsxx .* (Mi11 .* D11 + Mi12 .* D21 + Mi14 .* D41) + dQdsyy .* dFdsyy .* (Mi21 .* D11 + Mi22 .* D21 + Mi24 .* D41) + dQdsyy .* dFdsxy .* (Mi31 .* D11 + Mi32 .* D21 + Mi34 .* D41) + dQdsyy .* dFdszz .* (Mi41 .* D11 + Mi42 .* D21 + Mi44 .* D41)) - D44 ./ den .* (dQdszz .* dFdsxx .* (Mi11 .* D11 + Mi12 .* D21 + Mi14 .* D41) + dQdszz .* dFdsyy .* (Mi21 .* D11 + Mi22 .* D21 + Mi24 .* D41) + dQdszz .* dFdsxy .* (Mi31 .* D11 + Mi32 .* D21 + Mi34 .* D41) + dQdszz .* dFdszz .* (Mi41 .* D11 + Mi42 .* D21 + Mi44 .* D41)));
%     Dep22c = Mi21 .* (-D11 ./ den .* (dQdsxx .* dFdsxx .* (Mi11 .* D12 + Mi12 .* D22 + Mi14 .* D42) + dQdsxx .* dFdsyy .* (Mi21 .* D12 + Mi22 .* D22 + Mi24 .* D42) + dQdsxx .* dFdsxy .* (Mi31 .* D12 + Mi32 .* D22 + Mi34 .* D42) + dQdsxx .* dFdszz .* (Mi41 .* D12 + Mi42 .* D22 + Mi44 .* D42)) + D12 .* (1 - 1 ./ den .* (dQdsyy .* dFdsxx .* (Mi11 .* D12 + Mi12 .* D22 + Mi14 .* D42) + dQdsyy .* dFdsyy .* (Mi21 .* D12 + Mi22 .* D22 + Mi24 .* D42) + dQdsyy .* dFdsxy .* (Mi31 .* D12 + Mi32 .* D22 + Mi34 .* D42) + dQdsyy .* dFdszz .* (Mi41 .* D12 + Mi42 .* D22 + Mi44 .* D42))) - D14 ./ den .* (dQdszz .* dFdsxx .* (Mi11 .* D12 + Mi12 .* D22 + Mi14 .* D42) + dQdszz .* dFdsyy .* (Mi21 .* D12 + Mi22 .* D22 + Mi24 .* D42) + dQdszz .* dFdsxy .* (Mi31 .* D12 + Mi32 .* D22 + Mi34 .* D42) + dQdszz .* dFdszz .* (Mi41 .* D12 + Mi42 .* D22 + Mi44 .* D42))) + Mi22 .* (-D21 ./ den .* (dQdsxx .* dFdsxx .* (Mi11 .* D12 + Mi12 .* D22 + Mi14 .* D42) + dQdsxx .* dFdsyy .* (Mi21 .* D12 + Mi22 .* D22 + Mi24 .* D42) + dQdsxx .* dFdsxy .* (Mi31 .* D12 + Mi32 .* D22 + Mi34 .* D42) + dQdsxx .* dFdszz .* (Mi41 .* D12 + Mi42 .* D22 + Mi44 .* D42)) + D22 .* (1 - 1 ./ den .* (dQdsyy .* dFdsxx .* (Mi11 .* D12 + Mi12 .* D22 + Mi14 .* D42) + dQdsyy .* dFdsyy .* (Mi21 .* D12 + Mi22 .* D22 + Mi24 .* D42) + dQdsyy .* dFdsxy .* (Mi31 .* D12 + Mi32 .* D22 + Mi34 .* D42) + dQdsyy .* dFdszz .* (Mi41 .* D12 + Mi42 .* D22 + Mi44 .* D42))) - D24 ./ den .* (dQdszz .* dFdsxx .* (Mi11 .* D12 + Mi12 .* D22 + Mi14 .* D42) + dQdszz .* dFdsyy .* (Mi21 .* D12 + Mi22 .* D22 + Mi24 .* D42) + dQdszz .* dFdsxy .* (Mi31 .* D12 + Mi32 .* D22 + Mi34 .* D42) + dQdszz .* dFdszz .* (Mi41 .* D12 + Mi42 .* D22 + Mi44 .* D42))) - Mi23 .* D33 ./ den .* (dQdsxy .* dFdsxx .* (Mi11 .* D12 + Mi12 .* D22 + Mi14 .* D42) + dQdsxy .* dFdsyy .* (Mi21 .* D12 + Mi22 .* D22 + Mi24 .* D42) + dQdsxy .* dFdsxy .* (Mi31 .* D12 + Mi32 .* D22 + Mi34 .* D42) + dQdsxy .* dFdszz .* (Mi41 .* D12 + Mi42 .* D22 + Mi44 .* D42)) + Mi24 .* (-D41 ./ den .* (dQdsxx .* dFdsxx .* (Mi11 .* D12 + Mi12 .* D22 + Mi14 .* D42) + dQdsxx .* dFdsyy .* (Mi21 .* D12 + Mi22 .* D22 + Mi24 .* D42) + dQdsxx .* dFdsxy .* (Mi31 .* D12 + Mi32 .* D22 + Mi34 .* D42) + dQdsxx .* dFdszz .* (Mi41 .* D12 + Mi42 .* D22 + Mi44 .* D42)) + D42 .* (1 - 1 ./ den .* (dQdsyy .* dFdsxx .* (Mi11 .* D12 + Mi12 .* D22 + Mi14 .* D42) + dQdsyy .* dFdsyy .* (Mi21 .* D12 + Mi22 .* D22 + Mi24 .* D42) + dQdsyy .* dFdsxy .* (Mi31 .* D12 + Mi32 .* D22 + Mi34 .* D42) + dQdsyy .* dFdszz .* (Mi41 .* D12 + Mi42 .* D22 + Mi44 .* D42))) - D44 ./ den .* (dQdszz .* dFdsxx .* (Mi11 .* D12 + Mi12 .* D22 + Mi14 .* D42) + dQdszz .* dFdsyy .* (Mi21 .* D12 + Mi22 .* D22 + Mi24 .* D42) + dQdszz .* dFdsxy .* (Mi31 .* D12 + Mi32 .* D22 + Mi34 .* D42) + dQdszz .* dFdszz .* (Mi41 .* D12 + Mi42 .* D22 + Mi44 .* D42)));
%     Dep23c = Mi21 .* (-D11 ./ den .* (dQdsxx .* dFdsxx .* Mi13 .* D33 + dQdsxx .* dFdsyy .* Mi23 .* D33 + dQdsxx .* dFdsxy .* Mi33 .* D33 + dQdsxx .* dFdszz .* Mi43 .* D33) - D12 ./ den .* (dQdsyy .* dFdsxx .* Mi13 .* D33 + dQdsyy .* dFdsyy .* Mi23 .* D33 + dQdsyy .* dFdsxy .* Mi33 .* D33 + dQdsyy .* dFdszz .* Mi43 .* D33) - D14 ./ den .* (dQdszz .* dFdsxx .* Mi13 .* D33 + dQdszz .* dFdsyy .* Mi23 .* D33 + dQdszz .* dFdsxy .* Mi33 .* D33 + dQdszz .* dFdszz .* Mi43 .* D33)) + Mi22 .* (-D21 ./ den .* (dQdsxx .* dFdsxx .* Mi13 .* D33 + dQdsxx .* dFdsyy .* Mi23 .* D33 + dQdsxx .* dFdsxy .* Mi33 .* D33 + dQdsxx .* dFdszz .* Mi43 .* D33) - D22 ./ den .* (dQdsyy .* dFdsxx .* Mi13 .* D33 + dQdsyy .* dFdsyy .* Mi23 .* D33 + dQdsyy .* dFdsxy .* Mi33 .* D33 + dQdsyy .* dFdszz .* Mi43 .* D33) - D24 ./ den .* (dQdszz .* dFdsxx .* Mi13 .* D33 + dQdszz .* dFdsyy .* Mi23 .* D33 + dQdszz .* dFdsxy .* Mi33 .* D33 + dQdszz .* dFdszz .* Mi43 .* D33)) + Mi23 .* D33 .* (1 - 1 ./ den .* (dQdsxy .* dFdsxx .* Mi13 .* D33 + dQdsxy .* dFdsyy .* Mi23 .* D33 + dQdsxy .* dFdsxy .* Mi33 .* D33 + dQdsxy .* dFdszz .* Mi43 .* D33)) + Mi24 .* (-D41 ./ den .* (dQdsxx .* dFdsxx .* Mi13 .* D33 + dQdsxx .* dFdsyy .* Mi23 .* D33 + dQdsxx .* dFdsxy .* Mi33 .* D33 + dQdsxx .* dFdszz .* Mi43 .* D33) - D42 ./ den .* (dQdsyy .* dFdsxx .* Mi13 .* D33 + dQdsyy .* dFdsyy .* Mi23 .* D33 + dQdsyy .* dFdsxy .* Mi33 .* D33 + dQdsyy .* dFdszz .* Mi43 .* D33) - D44 ./ den .* (dQdszz .* dFdsxx .* Mi13 .* D33 + dQdszz .* dFdsyy .* Mi23 .* D33 + dQdszz .* dFdsxy .* Mi33 .* D33 + dQdszz .* dFdszz .* Mi43 .* D33));
end


%% Elastic or elasto-plastic operators
D11c = plc.*Dep11c + (1-plc).* De11c;
D12c = plc.*Dep12c + (1-plc).* De12c;
D13c = plc.*Dep13c + (1-plc).* De13c;
D21c = plc.*Dep21c + (1-plc).* De21c;
D22c = plc.*Dep22c + (1-plc).* De22c;
D23c = plc.*Dep23c + (1-plc).* De23c;
D31v = plv.*Dep31v + (1-plv).* De31v;
D32v = plv.*Dep32v + (1-plv).* De32v;
D33v = plv.*Dep33v + (1-plv).* De33v;
D11c = (1-plc2).*D11c + 0*plc2;
D12c = (1-plc2).*D12c;
D13c = (1-plc2).*D13c;
D21c = (1-plc2).*D21c + 0*plc2;
D22c = (1-plc2).*D22c;
D23c = (1-plc2).*D23c;
D31v = (1-plv2).*D31v;
D32v = (1-plv2).*D32v;
D33v = (1-plv2).*D33v;


end




