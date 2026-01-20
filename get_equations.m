% clc
% clear all
% load test2.mat 

% x momentum equation - left
    [~,B,en,Dmat] = get_data_v2 ( x_l , node , di , form , mat) ;   % Dmat  (3,3) = Dmat  (3,3)*2;

    Rx = Dmat(1,:)*B/deltaX  ;
    K(2*ij-1,en) = K(2*ij-1,en) - Rx ;
    Sigma_t = Dmat *B* dU(en);
    Residual(2*ij-1) = Residual(2*ij-1) - Sigma_t(1)/deltaX;
%     Residural(2*ij-1) = Residual(2*ij-1) - Rx* dU(en);
% Residual(2*ij-1)
    if ij == 38 &&  plast_it == 2
       Sigma_t 
       
       pause
    end

% x momentum equation - right
    [~,B,en,Dmat] = get_data_v2 ( x_r , node , di , form , mat) ;   % Dmat  (3,3) = Dmat  (3,3)*2;

    Rx = Dmat(1,:)*B/deltaX  ;
    K(2*ij-1,en) = K(2*ij-1,en) + Rx ;
    Sigma_t = Dmat *B* dU(en);
    Residual(2*ij-1) = Residual(2*ij-1) + Sigma_t(1)/deltaX;
%     Residual(2*ij-1) = Residual(2*ij-1) + Rx* dU(en);
% Residual(2*ij-1)
% x momentum equation - bot
    [~,B,en,Dmat] = get_data_v2 ( x_b , node , di , form , mat) ;   % Dmat  (3,3) = Dmat  (3,3)*2;
% B
% Dmat
    Rx = Dmat(3,:)*B/deltaY  ;
    K(2*ij-1,en) = K(2*ij-1,en) - Rx ;
    Sigma_t = Dmat *B* dU(en);
    Residual(2*ij-1) = Residual(2*ij-1) - Sigma_t(3)/deltaY;
%     1
%     Sigma_t(3)/deltaY
%     Residual(2*ij-1) = Residual(2*ij-1) - Rx* dU(en);
% Residual(2*ij-1)
% x momentum equation - top
    [~,B,en,Dmat] = get_data_v2 ( x_t , node , di , form , mat) ;  %  Dmat  (3,3) = Dmat  (3,3)*2;

    Rx = Dmat(3,:)*B/deltaY  ;
    K(2*ij-1,en) = K(2*ij-1,en) + Rx ;
    Sigma_t = Dmat *B* dU(en);
    Residual(2*ij-1) = Residual(2*ij-1) + Sigma_t(3)/deltaY;
%     Residual(2*ij-1) = Residual(2*ij-1) + Rx* dU(en);
% Residual(2*ij-1)

% y momentum equation - left
    [~,B,en,Dmat] = get_data_v2 ( x_l , node , di , form , mat) ;  %  Dmat  (3,3) = Dmat  (3,3)*2;

    Ry = Dmat(3,:)*B/deltaX  ;
    K(2*ij  ,en) = K(2*ij  ,en) - Ry ;
    Sigma_t = Dmat *B* dU(en);
    Residual(2*ij  ) = Residual(2*ij  ) - Sigma_t(3)/deltaX;
%     Residual(2*ij  ) = Residual(2*ij  ) - Ry* dU(en);

% y momentum equation - right
    [~,B,en,Dmat] = get_data_v2 ( x_r , node , di , form , mat) ;   % Dmat  (3,3) = Dmat  (3,3)*2;

    Ry = Dmat(3,:)*B/deltaX  ;
    K(2*ij  ,en) = K(2*ij  ,en) + Ry ;
    Sigma_t = Dmat *B* dU(en);
    Residual(2*ij  ) = Residual(2*ij  ) + Sigma_t(3)/deltaX;
%     Residual(2*ij  ) = Residual(2*ij  ) + Ry* dU(en);

% y momentum equation - bot
    [~,B,en,Dmat] = get_data_v2 ( x_b , node , di , form , mat) ;   % Dmat  (3,3) = Dmat  (3,3)*2;

    Ry = Dmat(2,:)*B/deltaY  ;
    K(2*ij  ,en) = K(2*ij  ,en) - Ry ;
    Sigma_t = Dmat *B* dU(en);
    Residual(2*ij  ) = Residual(2*ij  ) - Sigma_t(2)/deltaY;
%     Residual(2*ij  ) = Residual(2*ij  ) - Ry* dU(en);

% y momentum equation - top
    [~,B,en,Dmat] = get_data_v2 ( x_t , node , di , form , mat) ;   % Dmat  (3,3) = Dmat  (3,3)*2;

    Ry = Dmat(2,:)*B/deltaY  ;
    K(2*ij  ,en) = K(2*ij  ,en) + Ry ;
    Sigma_t = Dmat *B* dU(en);
    Residual(2*ij  ) = Residual(2*ij  ) + Sigma_t(2)/deltaY;
%     Residual(2*ij  ) = Residual(2*ij  ) + Ry* dU(en);
