%{
ff = dEps_t_l2' ; ff = ff(:); f1_l= Dmat_t1_l2*ff;

ff = dEps_t_r2' ; ff = ff(:); f1_r= Dmat_t1_r2*ff;

ff = dEps_t_b2' ; ff = ff(:); f1_b= Dmat_t1_b2*ff;

ff = dEps_t_t2' ; ff = ff(:); f1_t= Dmat_t1_t2*ff;

dSigma_1_ = [f1_l;f1_r;f1_b;f1_t];

ff = dEps_t_l2' ; ff = ff(:); f1_l= Dmat_t2_l2*ff;

ff = dEps_t_r2' ; ff = ff(:); f1_r= Dmat_t2_r2*ff;

ff = dEps_t_b2' ; ff = ff(:); f1_b= Dmat_t2_b2*ff;

ff = dEps_t_t2' ; ff = ff(:); f1_t= Dmat_t2_t2*ff;

dSigma_2_ = [f1_l;f1_r;f1_b;f1_t];


ff = dEps_t_l2' ; ff = ff(:); f1_l= Dmat_t3_l2*ff;

ff = dEps_t_r2' ; ff = ff(:); f1_r= Dmat_t3_r2*ff;

ff = dEps_t_b2' ; ff = ff(:); f1_b= Dmat_t3_b2*ff;

ff = dEps_t_t2' ; ff = ff(:); f1_t= Dmat_t3_t2*ff;

dSigma_3_ = [f1_l;f1_r;f1_b;f1_t];

ff = dEps_t_l2' ; ff = ff(:); f1_l= Dmat_t4_l2*ff;

ff = dEps_t_r2' ; ff = ff(:); f1_r= Dmat_t4_r2*ff;

ff = dEps_t_b2' ; ff = ff(:); f1_b= Dmat_t4_b2*ff;

ff = dEps_t_t2' ; ff = ff(:); f1_t= Dmat_t4_t2*ff;

dSigma_4_ = [f1_l;f1_r;f1_b;f1_t];

dSigma_1=dSigma_1_;dSigma_2=dSigma_2_;dSigma_3=dSigma_3_;dSigma_4=dSigma_4_;
dSigma_2 = [dSigma_1 dSigma_2 dSigma_3 dSigma_4];

%}

dd__=[dEps_t_l2;dEps_t_r2;dEps_t_b2;dEps_t_t2];
dSigma_2 = [D11.*dd__(:,1)+D12.*dd__(:,2)+D13.*dd__(:,3)+D14.*dd__(:,4) ...
    D21.*dd__(:,1)+D22.*dd__(:,2)+D23.*dd__(:,3)+D24.*dd__(:,4)...
    D31.*dd__(:,1)+D32.*dd__(:,2)+D33.*dd__(:,3)+D34.*dd__(:,4)...
    D41.*dd__(:,1)+D42.*dd__(:,2)+D43.*dd__(:,3)+D44.*dd__(:,4)];

