dt            = 1e+6*1;
% Dimension of the domain (it is simply a rectangular region L x W)
L = 1*150000 ;
D = 1*150000 ;
D = 1*150000/1 ;

% MLS_params
% Inputs special to meshless methods, domain of influence
shape = 'circle' ;         % shape of domain of influence
dmax  = 1.15 ;              % radius = dmax * nodal spacing
form  = 'cubic_spline' ;   % using cubic spline weight function


% Node density defined by number of nodes along two directions
% nnx = 30 ;
nny = nnx/1 ;

% node generation, node is of size (nnode x 2)
node = square_node_array([0 0],[L 0],[L D],[0 D],nnx,nny);
node(abs(node(:,2)-D)<0.0001,2)=D;
node(abs(node(:,1)-L)<0.0001,1)=L;
numnode = size(node,1);

deltaX = L/(nnx-1);
deltaY = D/(nny-1);
delta  = max(deltaX,deltaY);
di     = ones(1,numnode)*dmax*delta ;



%% non-uniform grid
%{
    densy = [ 20  9/20 ; 
            10  2/20 ;
            20  9/20 ]; 
   densx = [151 10/40;
             151 10/40;
             151 10/40;
             151 10/40];
    nnx = 1 + sum(densx(:,1)) ;
    nny = 1 + sum(densy(:,1)) ;
    pt1 = [0 0]; pt2 = [L 0];  pt3 = [L D]; pt4 = [0 D];
    node = square_node_array_ir(pt1,pt2,pt3,pt4,densx,densy);
% Nodal Support
    di_x = ones(1,nnx);
    di_y = ones(1,nny);
    deltaX = zeros(1,size(densx,1));
    deltaY = zeros(1,size(densy,1));
    for i = 1 : size(densx,1)
                deltaX(1,i) = densx(i,2)*L/densx(i,1);
    end
di_x(1,1:densx(1,1)) = ones(1,densx(1,1))*deltaX(1,1);
    for i = 2 : size(densx,1)
        di_x(1,1+sum(densx(1:i-1,1))) = max(deltaX(1,i-1),deltaX(1,i));
        di_x(1,2+sum(densx(1:i-1,1)):1+sum(densx(1:i,1))) = ...
            ones(1,densx(i,1))*deltaX(1,i);
    end

    for i = 1 : size(densy,1)
        deltaY(1,i) = densy(i,2)*D/densy(i,1);
    end

    di_y(1,1:densy(1,1)) = ones(1,densy(1,1))*deltaY(1,1);

    for i = 2 : size(densy,1)
        di_y(1,1+sum(densy(1:i-1,1))) = (deltaY(1,i-1)+deltaY(1,i))/2;
        di_y(1,2+sum(densy(1:i-1,1)):1+sum(densy(1:i,1))) = ...
            ones(1,densy(i,1))*deltaY(1,i);
    end

    di = zeros(1,nnx*nny);
    num = 0;
    for j = 1 : nny
        for i = 1 : nnx
            num = num + 1;
            di(1,num) = max( di_x(1,i) , di_y(1,j) ) * dmax;
        end
    end
figure
plot(node(:,1),node(:,2),'r.')
deltaY= min(deltaY) ;
%}
%%
% Meshing, Q4 elements
inc_u = 1;
inc_v = nnx;
node_pattern = [ 1 2 nnx+2 nnx+1 ];
element = make_elem(node_pattern,nnx-1,nny-1,inc_u,inc_v);

% find connected nodes to each element:
[conn_nodes] = find_conn_nodes (node , element ) ;

fi_ = zeros(1,length(node))+deltaY;
fi_2 = [fi_' fi_' fi_' fi_'] ; 
%% T3 mesh
%{
% load meshT3.mat 'node' 'element'
% node(:,1) = [] ; 
% node(:,end) = []; 
% element(:,1:5) = [] ; 
% numnode = size(node,1);
% node = node*L; 
% figure
% plot_mesh(node,element,'T3','k-')
% axis equal
    


% [node,element,numnode]=femTriangularMeshGenerator(L,D,nnx,nny);
% nnx = ii_ ; nny = ii_ ; 
load data_13_1.mat node element
numnode = size(node,1) ; 
% figure
% hold on
% % plot_mesh(node,element,'T3','k-')
% axis equal
plot(node(:,1),node(:,2),'rsq')
deltaX = L/(nnx-1);
deltaY = D/(nny-1);
delta  = max(deltaX,deltaY);
di     = ones(1,numnode)*dmax*delta ;

% find connected nodes to each element:
[conn_nodes] = find_conn_nodes (node , element ) ;

fi_       = zeros(1,numnode) ;
fi_2      = zeros(numnode,4)+1e120 ;
maxS = 0 ;
for ig = 1 : size(node,1)
    ik=conn_nodes(ig,2:conn_nodes(ig,1)+1) ;
    nn=element(ik,:); nn = unique(nn(:));
%     plot(node(ig,1),node(ig,2),'bsq')
%     pause(0.1)
    l_min = 1e120 ; 
    xx_ = node(ig,:);
    if (abs(node(ig,1)-0)>0.0001 && abs(node(ig,1)-L)>0.0001  && abs(node(ig,2)-D)>0.0001 && abs(node(ig,2)-0)>0.0001 )
        gh = node(nn,:)-repmat(xx_,size(nn,1),1);
        x_L = abs(min(gh(:,1)));
        x_R = abs(max(gh(:,1)));
        y_B = abs(min(gh(:,2)));
        y_T = abs(max(gh(:,2)));
        
        for ig2 = 1 : length(nn)
            if nn(ig2) ~= ig
                l_min = min(l_min , norm(node(nn(ig2),:)-xx_));
            end
        end
%         pause
    fi_2(ig,:) = [x_L x_R y_B y_T];
    end    
    
    
    fi_(ig) = l_min ; 
end
deltaY = 1*min(fi_) ; 
deltaX = 1*min(fi_) ; 
% deltaY = L/50 ;
% deltaX = L/50 ;
%}

%%

% material properties
eta_vp = 2.5e6;       % Kelvin element viscosity
eta0   = 2e50;        % Shear viscosity (DEACTIVATED)
coh0   = 1e17;     % Cohesion
phi0   = 30*pi/180;   % Friction angle
psi0   = 10*pi/180;   % Dilatancy angle
sin_phi = sin(phi0); 
cos_phi = cos(phi0); 
sin_psi = sin(psi0);
sin_psi = 0 ; 
sin_phi = 0.0008*0 ;
h      = -0*0.01;     % Hardening/softening modulus

K1 = 50*1e9 ;
G1 = 30*1e9 ; 
K2 = 50*1e9 ; 
G2 = 30*1e9 ;
Gve1  = 1./(dt./eta0 + 1./G1);
VE1   = Gve1./G1;
Gve2  = 1./(dt./eta0 + 1./G2);
VE2   = Gve2./G2;


mat = [ K1 Gve1 K2 Gve2 ] ; 

dynamic = 1 ; 
rho =2700*1; 
% incr0 = 5e-6/1; % applied displacement at each increment
incr0 = 2e-9*dt; % applied displacement at each increment

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
uy = 2e-9 ;
lambda__ = 0*0.67 ; 
Mat_ = [ mu0 V0 a_ L_ ]; 

p0 = 5e6 ; % initial pressure


fault_width = deltaY/2 ;
% fault_width = L/1000 ; 
n_faults = 1 ; 
aaa = .25*L/2 ;
alph__ = 0.1/1*0 ; 

% fault_width = L/640/2 ; 