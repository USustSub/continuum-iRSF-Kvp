[sin_phi_0_a_n,Vp_a_n,Eii_p_a_n,theta_new_a_n] = get_mu_new2 ( dg_a_n , dQdsigma_(inx,:) , dt , inx , deltaY/2 , b_( inx) , theta_( inx) );
sin_phi_a_n = sin_phi_0_a_n ;
Fa_n = 0 - sqrt(J2_(inx)) + P_(inx).*sin_phi_a_n + dg_a_n.*Gve_(inx)  + dg_a_n*eta_vp/dt ; 
[sin_phi_0_b_n,Vp_b,Eii_p_n_b,theta_new_b_n] = get_mu_new2 ( dg_b_n , dQdsigma_(inx,:) , dt , inx , deltaY/2 , b_( inx) , theta_( inx) );
sin_phi_b_n = sin_phi_0_b_n ;
Fb_n = 0 - sqrt(J2_(inx)) + P_(inx).*sin_phi_b_n + dg_b_n.*Gve_(inx)  + dg_b_n*eta_vp/dt ; 
