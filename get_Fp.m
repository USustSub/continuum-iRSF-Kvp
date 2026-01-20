[sin_phi_0_p,Vp_p,Eii_p_p,theta_new_p] = get_mu_new2 ( dg_p , dQdsigma_(inx(iii_),:) , dt , inx(iii_) , deltaY/2 , b_( inx(iii_)) , theta_( inx(iii_)) );
sin_phi_p = sin_phi_0_p ;
Fp = 0 - sqrt(J2_(inx(iii_))) + P_(inx(iii_)).*sin_phi_p + dg_p.*Gve_(inx(iii_))  + dg_p*eta_vp/dt ; 
