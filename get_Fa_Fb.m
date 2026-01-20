[sin_phi_0_a,Vp_a,Eii_p_a,theta_new_a] = get_mu_new2 ( dg_a , dQdsigma_(inx(iii_),:) , dt , inx(iii_) , deltaY/2 , b_( inx(iii_)) , theta_( inx(iii_)) );
sin_phi_a = sin_phi_0_a ;
Fa = 0 - sqrt(J2_(inx(iii_))) + P_(inx(iii_)).*sin_phi_a + dg_a.*Gve_(inx(iii_))  + dg_a*eta_vp/dt ; 
[sin_phi_0_b,Vp_b,Eii_p_b,theta_new_b] = get_mu_new2 ( dg_b , dQdsigma_(inx(iii_),:) , dt , inx(iii_) , deltaY/2 , b_( inx(iii_)) , theta_( inx(iii_)) );
sin_phi_b = sin_phi_0_b ;
Fb = 0 - sqrt(J2_(inx(iii_))) + P_(inx(iii_)).*sin_phi_b + dg_b.*Gve_(inx(iii_))  + dg_b*eta_vp/dt ; 
