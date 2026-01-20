% function  [alpha,LSsuccess] = M2Di_EP13_LineSearchEP(dU,ddU,NumUx,NumUyG,Gvec,VEc,Kc,Cc,phic,psic,Gvev,VEv,Kv,Cv,phiv,psiv,dx,dy,ncx,ncy,Pc0,Pv0,Txxc0,Tyyc0,Txyc0,Tzzc0,Txxv0,Tyyv0,Txyv0,Tzzv0, vm, BC, sin_phi, cos_phi, sin_psi, safe, res0, rho0, gy, pconf, eta_vp, dt, alpha_min, Ep_acc_c0, Ep_acc_v0, hc, hv, Cc0, Cv0, imp_hard  )

% Line search
dUo      = dU;
alpha    = 1;
alphamax = alpha;

alpha_min     = 0.25;   % Minimum correction step
alphamin = alpha_min;

nls      = 11;

dalpha   = (alphamax-alphamin)/(nls-1);
alphav   = alphamin:dalpha:alphamax;
resi     = zeros(nls,1);
LIsuccessc = 1;
LIsuccessv = 1;

for ils = 1:nls
    
    % Apply Newton step
    dU       = dUo + alphav(ils)*ddU;
    
    
% get residual vector
Residual = 0*Residual ; 
    get_res ;
    
    
    nR = norm(Residual,1)/length(Residual);
    resi(ils) = nR;
%     if ils==1,  
%         fprintf('LS It. %03d - alpha = %2.2f --- ||F|| = %2.8e\n', ils, alphav(ils), nR); 
%     end
end
%     figure
%     plot(1:length(resi),resi,'rsq-')
%     pause

    
[val,ind] = min(resi);
alpha     = alphav(ind);
% fprintf('Selected alpha = %1.2f yields||F|| = %2.8e\n',alpha, resi(ind));
LSsuccess = 1;
