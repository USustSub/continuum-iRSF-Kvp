[sin_phi_0_p_n,Vp_p_n,Eii_p_p_n,theta_new_p_n] = get_mu_new2 ( dg_p_n , dQdsigma_(inx,:) , dt , inx , deltaY/2 , b_( inx) , theta_( inx) );
sin_phi_p_n = sin_phi_0_p_n ;
Fp_n = 0 - sqrt(J2_(inx)) + P_(inx).*sin_phi_p_n + dg_p_n.*Gve_(inx)  + dg_p_n*eta_vp/dt ; 
