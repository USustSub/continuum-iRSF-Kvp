
dQdsigma_xx = dQdsigma_(:,1) ; 
dQdsigma_yy = dQdsigma_(:,2) ; 
dQdsigma_xy = dQdsigma_(:,3) ; 
dQdsigma_zz = dQdsigma_(:,4) ; 
ddY = deltaY/2 ; 


Der = Gve_ + eta_vp./dt + (P_.*a_.*ddY.*exp((mu0 + b_.*log((V0.*theta_n)./L_))./a_).*((dQdsigma_xx./(3.*dt) + dQdsigma_yy./(3.*dt) - (2.*dQdsigma_zz)./(3.*dt)).*((dg_2.*dQdsigma_xx)./(3.*dt) + (dg_2.*dQdsigma_yy)./(3.*dt) - (2.*dg_2.*dQdsigma_zz)./(3.*dt)) + (dQdsigma_xx./(3.*dt) - (2.*dQdsigma_yy)./(3.*dt) + dQdsigma_zz./(3.*dt)).*((dg_2.*dQdsigma_xx)./(3.*dt) - (2.*dg_2.*dQdsigma_yy)./(3.*dt) + (dg_2.*dQdsigma_zz)./(3.*dt)) + (dQdsigma_yy./(3.*dt) - (2.*dQdsigma_xx)./(3.*dt) + dQdsigma_zz./(3.*dt)).*((dg_2.*dQdsigma_yy)./(3.*dt) - (2.*dg_2.*dQdsigma_xx)./(3.*dt) + (dg_2.*dQdsigma_zz)./(3.*dt)) + (2.*dg_2.*dQdsigma_xy.^2)./dt.^2))./(2.*V0.*((ddY.^2.*exp((2.*mu0 + 2.*b_.*log((V0.*theta_n)./L_))./a_).*(((dg_2.*dQdsigma_xx)./(3.*dt) + (dg_2.*dQdsigma_yy)./(3.*dt) - (2.*dg_2.*dQdsigma_zz)./(3.*dt)).^2./2 + ((dg_2.*dQdsigma_xx)./(3.*dt) - (2.*dg_2.*dQdsigma_yy)./(3.*dt) + (dg_2.*dQdsigma_zz)./(3.*dt)).^2./2 + ((dg_2.*dQdsigma_yy)./(3.*dt) - (2.*dg_2.*dQdsigma_xx)./(3.*dt) + (dg_2.*dQdsigma_zz)./(3.*dt)).^2./2 + (dg_2.^2.*dQdsigma_xy.^2)./dt.^2))./V0.^2 + 1).^(1./2).*(((dg_2.*dQdsigma_xx)./(3.*dt) + (dg_2.*dQdsigma_yy)./(3.*dt) - (2.*dg_2.*dQdsigma_zz)./(3.*dt)).^2./2 + ((dg_2.*dQdsigma_xx)./(3.*dt) - (2.*dg_2.*dQdsigma_yy)./(3.*dt) + (dg_2.*dQdsigma_zz)./(3.*dt)).^2./2 + ((dg_2.*dQdsigma_yy)./(3.*dt) - (2.*dg_2.*dQdsigma_xx)./(3.*dt) + (dg_2.*dQdsigma_zz)./(3.*dt)).^2./2 + (dg_2.^2.*dQdsigma_xy.^2)./dt.^2).^(1./2));
