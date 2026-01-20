% % tic
% dSigma_1 = [diag(Dmat_t1_l(:,1:4)*dEps_t_l(:,:)');
%         diag(Dmat_t1_r(:,1:4)*dEps_t_r(:,:)');
%         diag(Dmat_t1_b(:,1:4)*dEps_t_b(:,:)');
%         diag(Dmat_t1_t(:,1:4)*dEps_t_t(:,:)');];
% dSigma_2 = [diag(Dmat_t2_l(:,1:4)*dEps_t_l(:,:)');
%         diag(Dmat_t2_r(:,1:4)*dEps_t_r(:,:)');
%         diag(Dmat_t2_b(:,1:4)*dEps_t_b(:,:)');
%         diag(Dmat_t2_t(:,1:4)*dEps_t_t(:,:)');];
% dSigma_3 = [diag(Dmat_t3_l(:,1:4)*dEps_t_l(:,:)');
%         diag(Dmat_t3_r(:,1:4)*dEps_t_r(:,:)');
%         diag(Dmat_t3_b(:,1:4)*dEps_t_b(:,:)');
%         diag(Dmat_t3_t(:,1:4)*dEps_t_t(:,:)');];
% dSigma_4 = [diag(Dmat_t4_l(:,1:4)*dEps_t_l(:,:)');
%         diag(Dmat_t4_r(:,1:4)*dEps_t_r(:,:)');
%         diag(Dmat_t4_b(:,1:4)*dEps_t_b(:,:)');
%         diag(Dmat_t4_t(:,1:4)*dEps_t_t(:,:)');];

dd__=[dEps_t_l;dEps_t_r;dEps_t_b;dEps_t_t];
dSigma_ = [D11.*dd__(:,1)+D12.*dd__(:,2)+D13.*dd__(:,3)+D14.*dd__(:,4) ...
    D21.*dd__(:,1)+D22.*dd__(:,2)+D23.*dd__(:,3)+D24.*dd__(:,4)...
    D31.*dd__(:,1)+D32.*dd__(:,2)+D33.*dd__(:,3)+D34.*dd__(:,4)...
    D41.*dd__(:,1)+D42.*dd__(:,2)+D43.*dd__(:,3)+D44.*dd__(:,4)];
% % % toc
%%
%{
% % tic
ff = dEps_t_l' ; ff = ff(:); f1_l= Dmat_t1_l2*ff;

ff = dEps_t_r' ; ff = ff(:); f1_r= Dmat_t1_r2*ff;

ff = dEps_t_b' ; ff = ff(:); f1_b= Dmat_t1_b2*ff;

ff = dEps_t_t' ; ff = ff(:); f1_t= Dmat_t1_t2*ff;

dSigma_1_ = [f1_l;f1_r;f1_b;f1_t];

ff = dEps_t_l' ; ff = ff(:); f1_l= Dmat_t2_l2*ff;

ff = dEps_t_r' ; ff = ff(:); f1_r= Dmat_t2_r2*ff;

ff = dEps_t_b' ; ff = ff(:); f1_b= Dmat_t2_b2*ff;

ff = dEps_t_t' ; ff = ff(:); f1_t= Dmat_t2_t2*ff;

dSigma_2_ = [f1_l;f1_r;f1_b;f1_t];


ff = dEps_t_l' ; ff = ff(:); f1_l= Dmat_t3_l2*ff;

ff = dEps_t_r' ; ff = ff(:); f1_r= Dmat_t3_r2*ff;

ff = dEps_t_b' ; ff = ff(:); f1_b= Dmat_t3_b2*ff;

ff = dEps_t_t' ; ff = ff(:); f1_t= Dmat_t3_t2*ff;

dSigma_3_ = [f1_l;f1_r;f1_b;f1_t];

ff = dEps_t_l' ; ff = ff(:); f1_l= Dmat_t4_l2*ff;

ff = dEps_t_r' ; ff = ff(:); f1_r= Dmat_t4_r2*ff;

ff = dEps_t_b' ; ff = ff(:); f1_b= Dmat_t4_b2*ff;

ff = dEps_t_t' ; ff = ff(:); f1_t= Dmat_t4_t2*ff;

dSigma_4_ = [f1_l;f1_r;f1_b;f1_t];

dSigma_1=dSigma_1_;
dSigma_2=dSigma_2_;
dSigma_3=dSigma_3_;
dSigma_4=dSigma_4_;

dSigma_ = [dSigma_1 dSigma_2 dSigma_3 dSigma_4];

% norm(dSigma_1-dSigma_1_)
% norm(dSigma_2-dSigma_2_)
% norm(dSigma_3-dSigma_3_)
% norm(dSigma_4-dSigma_4_)
% tic
% Mul1 = zeros(size(Dmat_t1_l,1),1);
% % Mul2 = zeros(size(Dmat_t1_r,1),1);
% % Mul3 = zeros(size(Dmat_t1_b,1),1);
% % Mul4 = zeros(size(Dmat_t1_t,1),1);
% for ii = 1 : size(Dmat_t1_l,1)
%    Mul1(ii) = Dmat_t1_l(ii,1:4)*dEps_t_l(ii,:)';
% %    Mul2(ii) = Dmat_t1_r(ii,1:4)*dEps_t_r(ij,:)';
% %    Mul3(ii) = Dmat_t1_b(ii,1:4)*dEps_t_r(ij,:)';
% %    Mul4(ii) = Dmat_t1_t(ii,1:4)*dEps_t_r(ij,:)';
% %    Mul1(ii) = Dmat_t1_l(ii,1:4)*dEps_t_r(ij,:)';
% %    Mul2(ii) = Dmat_t1_r(ii,1:4)*dEps_t_r(ij,:)';
% %    Mul3(ii) = Dmat_t1_b(ii,1:4)*dEps_t_r(ij,:)';
% %    Mul4(ii) = Dmat_t1_t(ii,1:4)*dEps_t_r(ij,:)';
% end
% norm(Mul1-f2)

% % toc
% pause
%}
%%
% diag(Dmat_t1_l(:,1:4)*dEps_t_l(:,:)')
% out1 = sum(Dmat_t1_l(:,1:4).*dEps_t_l(:,:)',2); % Elementwise

% %%
% F = rand(10);
% B = rand(10);
% out1 = sum(F*B.*F,2); % Elementwise
% out2 = diag(F*B*F');
% isequal(out1,out2)