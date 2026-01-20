% Dimension of the domain (it is simply a rectangular region L x W)
L = 150000/1 ;
D = 150000/1 ;

% MLS_params
% Inputs special to meshless methods, domain of influence
shape = 'circle' ;         % shape of domain of influence
dmax  = 1.2 ;              % radius = dmax * nodal spacing
form  = 'cubic_spline' ;   % using cubic spline weight function

% Node density defined by number of nodes along two directions
% nnx = 30 ;
nny = nnx ;

% node generation, node is of size (nnode x 2)
node = square_node_array([0 0],[L 0],[L D],[0 D],nnx,nny);
node(abs(node(:,2)-D)<0.0001,2)=D;
numnode = size(node,1);

deltaX = L/(nnx-1);
deltaY = D/(nny-1);
delta  = max(deltaX,deltaY);
di     = ones(1,numnode)*dmax*delta ;

% Meshing, Q4 elements
inc_u = 1;
inc_v = nnx;
node_pattern = [ 1 2 nnx+2 nnx+1 ];
element = make_elem(node_pattern,nnx-1,nny-1,inc_u,inc_v);



midY = max(node(:,2))/2;
dY = node(nnx+1,2)-node(1,2);
fault_nodes = find(abs(node(:,2)-midY)<deltaY);
node_id = zeros(size(node,1),1);
node_id(fault_nodes,1) = 1 ; 
figure
plot(node(:,1),node(:,2),'r.')
hold on
plot(node(node_id~=0,1),node(node_id~=0,2),'bsq')
axis equal

    element_des = [   ] ; 
    for iel = 1 : size(element,1)
        sctr = element(iel,:) ; 
        node_sctr = node(sctr,:);
        if all(ismember(sctr,fault_nodes))
            element_des = [ element_des ; iel ] ; 
        end
    end
element_id = zeros(size(element,1),1);
element_id(element_des,1) = 1 ; 
    


% material properties
eta_vp = 2.5e8;       % Kelvin element viscosity
eta0   = 2e50;        % Shear viscosity (DEACTIVATED)
coh0   = 1000000.0 ;     % Cohesion
phi0   = 30*pi/180;   % Friction angle
psi0   = 10*pi/180;   % Dilatancy angle
mu_ = sin(phi0); 
cos_phi = cos(phi0); 
sin_psi = sin(psi0);
h      = -0*0.01;     % Hardening/softening modulus
sin_phi = sin(phi0); 

K1 = 50*1e9 ; % N/m2
G1 = 30*1e9 ; % N/m2
K2 = 50*1e9; 
G2 = 30*1e9;
Gve1  = 1./(dt./eta0 + 1./G1);
VE1   = Gve1./G1;
Gve2  = 1./(dt./eta0 + 1./G2);
VE2   = Gve2./G2;
rho = 2700 ; % kg/m3

mat = [ K1 Gve1 K2 Gve2 ] ; 


% RSF inputs
mu0 = 0.2 ; 
V0 = 4e-9 ;  % m/s
a_ = 0.01 ;
L_ = 0.01 ;
b1_ = 0.001 ; 
b2_ = 0.017 ; 
theta_ = L_/V0*exp(-1);
LL1 = 32000 ;
LL2 = 108000 ; 

