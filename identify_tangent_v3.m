    
function [Dmat_t,Dmat,G,K] = identify_tangent_v3 ( node_c , mat , id ) 


K1 = mat(1) ;
G1 = mat(2) ; 
K2 = mat(3); 
G2 = mat(4) ;
% if (node_c(1)>=0.0 && node_c(1)<=3 && node_c(2)>=0.0 && node_c(2)<=3)
if id == 1
        K = K2 ; % K2
        G = G2;  % G2
elseif id == 0
        K = K1 ; % K1
        G = G1;  % G1
end

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