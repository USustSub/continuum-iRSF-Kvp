function dQdsigma = get_dQdSigma ( TSigma , sin_psi , J2 ) 
    dQdsigma(1) = TSigma(1) .* (J2).^ (-0.5) / 2 + sin_psi / 3;
    dQdsigma(2) = TSigma(2) .* (J2).^ (-0.5) / 2 + sin_psi / 3;
    dQdsigma(3) = TSigma(3) .* (J2).^ (-0.5);
    dQdsigma(4) = TSigma(4) .* (J2).^ (-0.5) / 2 + sin_psi / 3;

