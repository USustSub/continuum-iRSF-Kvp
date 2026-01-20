% Dimension of the domain (it is simply a rectangular region L x W)
L = 65 ;
D = 65 ;

% MLS_params
% Inputs special to meshless methods, domain of influence
shape = 'circle' ;         % shape of domain of influence
dmax  = 1.7 ;              % radius = dmax * nodal spacing
form  = 'cubic_spline' ;   % using cubic spline weight function


% Node density defined by number of nodes along two directions
nnx = 30 ;
nny = 30 ;

% node generation, node is of size (nnode x 2)
node = square_node_array([0 0],[L 0],[L D],[0 D],nnx,nny);
node(abs(node(:,2)-D)<0.0001,2)=D;

% read_node ;
% node(:,1) = [] ;
% element(:,1) = [] ; 
numnode = size(node,1);
%%
% deltaX = L/(nnx-1);
% deltaY = D/(nny-1);
% delta  = max(deltaX,deltaY);
% di     = ones(1,numnode)*dmax*delta ;

% Meshing, Q4 elements
inc_u = 1;
inc_v = nnx;
node_pattern = [ 1 2 nnx+2 nnx+1 ];
element = make_elem(node_pattern,nnx-1,nny-1,inc_u,inc_v);


% material properties
eta_vp = 1.0 ;       % Kelvin element viscosity
eta0   = 2e50;        % Shear viscosity (DEACTIVATED)
coh0   = 0.2 ;     % Cohesion
phi0   = 30*pi/180;   % Friction angle
psi0   = 10*pi/180;   % Dilatancy angle
sin_phi = sin(phi0); 
cos_phi = cos(phi0); 
sin_psi = sin(psi0);
h      = -0*0.01;     % Hardening/softening modulus

K1 = 10000 ;
G1 = 7500 ; 
K2 = 10000/10; 
G2 = 7500/10 ;
Gve1  = 1./(dt./eta0 + 1./G1);
VE1   = Gve1./G1;
Gve2  = 1./(dt./eta0 + 1./G2);
VE2   = Gve2./G2;


mat = [ K1 Gve1 K2 Gve2 ] ; 
