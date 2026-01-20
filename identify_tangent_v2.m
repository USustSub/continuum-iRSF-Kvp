    
function [Dmat_t,Dmat,G,K] = identify_tangent_v2 ( node_c , mat , dynamic , node) 


K1 = mat(1) ;
G1 = mat(2) ; 
K2 = mat(3); 
G2 = mat(4) ;
L = max(node(:,1)) ; 
D = max(node(:,2)) ;
% if (node_c(1)>=0.0 && node_c(1)<=3 && node_c(2)>=0.0 && node_c(2)<=3)
% if dynamic == 0
    if (sqrt(node_c(1)^2+node_c(2)^2)<=0.05 )
            K = K2 ; % K2
            G = G2;  % G2
    else
            K = K1 ; % K1
            G = G1;  % G1
    end
%     if (sqrt((node_c(1)-L)^2+(node_c(2)-0)^2)<=0.05 )
%             K = K2 ; % K2
%             G = G2;  % G2
%     else
%             K = K1 ; % K1
%             G = G1;  % G1
%     end
% else
%     if (abs(node_c(2)-max(node(:,2))/2)<=3.8462e+03*1. )
%             K = K2 ; % K2
%             G = G2;  % G2
%     else
%             K = K1 ; % K1
%             G = G1;  % G1
%     end
% end
Dmat_t = [
    K+4/3*G K-2/3*G 0           K-2/3*G ;
    K-2/3*G K+4/3*G 0           K-2/3*G ;
    0         0         G   0 ;
    K-2/3*G K-2/3*G 0           K+4/3*G ;  
    ];

Dmat = [
    K+4/3*G K-2/3*G 0            ;
    K-2/3*G K+4/3*G 0            ;
    0         0         G    ;
];