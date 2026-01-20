%% Update friction coefficient
function [sin_phi_,Vp,Eii_p,theta_n] = get_mu_new ( dg_2 , dQdsigma_ , dt , inx , deltaY , b_ ,  theta_ , Mat_)
mu0 = Mat_(1) ; 
V0 = Mat_(2) ;  % m/s
a_ = Mat_(3) ;
L_ = Mat_(4) ;

        sin_phi_ = 0 *dg_2+ NaN; 
        
        depsp = [dg_2.*dQdsigma_(:,1) dg_2.*dQdsigma_(:,2) dg_2.*dQdsigma_(:,3) dg_2.*dQdsigma_(:,4)]/dt;
        
        Ekkc_ = depsp(:,1) + depsp(:,2) + depsp(:,4) ;
        depsp_d = depsp ;
        depsp_d(:,[1 2 4]) = depsp(:,[1 2 4]) - 1/3*Ekkc_ ; 
        Eii_p  = sqrt(1/2*(depsp_d(:,1)).^2 + 1/2*(depsp_d(:,2)).^2 + 1/2*(depsp_d(:,4)).^2 + depsp_d(:,3).^2);
%         Eii_p  = sqrt(depsp_d(:,1).^2 + depsp_d(:,3).^2 ) ;
        
        Vp = 2*Eii_p*deltaY ;
%          max(Eii_p(inx))
%         max(Vp(inx))
%      max(mu_gp(inx))   
%         sparse(depsp)
%         Vp(Vp== 0) = 1e-20;
        
% update theta 
        S1 = 1+theta_/dt;
        S2 = (1/dt + Vp/L_);
        theta_n = S1./S2;
%         theta_n = theta_ ; 
% steady state 
%         theta_n = L_./Vp ; 
% update friction coefficient
        param = Vp./2/V0.*exp((mu0+b_.*log(theta_n*V0/L_))/a_);
        mu_gp = a_*asinh(param) ;
%         max(Eii_p(inx))
%         if all(F__<0)
%             mu_gp = 0*mu_gp + 0.08;
%         end
%         mu_gp(mu_gp==0) = sin_phi ;
        sin_phi_(inx) = mu_gp(inx) ; 
%         sin_phi_(sin_phi_==0) = 0.0001;
%         theta_ = theta_n ; 
end

