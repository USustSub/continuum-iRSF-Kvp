%% 2D MLS based FDM 
% By Mohsen Goudarzi ( m.goudarzi@uu.nl ) 
% ----------------------------------------------------------
clear all; close all ; 
clc
tic;           % help us to see the time required for each step

% Dimension of the domain (it is simply a rectangular region L x W)
L = 4 ;
D = 4 ;

% Material properties
E  = 100 ;
nu = 0.3 ;


% Inputs special to meshless methods, domain of influence
shape = 'circle' ;         % shape of domain of influence
dmax  = 2.0 ;              % radius = dmax * nodal spacing
form  = 'cubic_spline' ;   % using cubic spline weight function


% Node density defined by number of nodes along two directions
nnx = 20 ;
nny = 20 ;

% node generation, node is of size (nnode x 2)
node = square_node_array([0 0],[L 0],[L D],[0 D],nnx,nny);

load node_interface.mat 'node' 'element'

numnode = size(node,1);

deltaX = L/(nnx-1);
deltaY = D/(nny-1);
delta  = max(deltaX,deltaY);
di     = ones(1,numnode)*dmax*delta ;



% Meshing, Q4 elements
% inc_u = 1;
% inc_v = nnx;
% node_pattern = [ 1 2 nnx+2 nnx+1 ];
% element = make_elem(node_pattern,nnx-1,nny-1,inc_u,inc_v);



figure
hold on
plot(node(:,1),node(:,2),'rsq')
for icel = 1 : size(element,1)
    sc = element(icel,:) ; 
    node_c = mean(node(sc,:));
    
    x_l = node_c - [deltaX/2 0 ] ;
    x_r = node_c + [deltaX/2 0 ] ;    
    
    y_b = node_c - [0 deltaX/2 ] ;
    y_t = node_c + [0 deltaX/2 ] ;    

end
axis equal


 
%%
Consti;

total_num_node = numnode ;
K = sparse(2*total_num_node,2*total_num_node);
F = zeros(2*total_num_node,1);

% loop over nodes
for ij = 1 : size(node,1)
    node_c = node(ij,:) ; 
    
if (node_c(1) > 0 && node_c(1) < L && node_c(2) < L && node_c(2) > 0 )
    x_l = node_c - [deltaX/2 0 ] ;
    x_r = node_c + [deltaX/2 0 ] ;    
    
    x_b = node_c - [0 deltaY/2 ] ;
    x_t = node_c + [0 deltaY/2 ] ;    

% x momentum equation - left
    [~,B,en] = get_data ( x_l , node , di , form ) ;
    Rx = Dmat(1,:)*B/deltaX  ;
    K(2*ij-1,en) = K(2*ij-1,en) - Rx ;

% x momentum equation - right
    [~,B,en] = get_data ( x_r , node , di , form) ;
    Rx = Dmat(1,:)*B/deltaX  ;
    K(2*ij-1,en) = K(2*ij-1,en) + Rx ;

% x momentum equation - bot
    [~,B,en] = get_data ( x_b , node , di , form) ;
    Rx = Dmat(3,:)*B/deltaY  ;
    K(2*ij-1,en) = K(2*ij-1,en) - Rx ;

% x momentum equation - top
    [~,B,en] = get_data ( x_t , node , di , form) ;
    Rx = Dmat(3,:)*B/deltaY  ;
    K(2*ij-1,en) = K(2*ij-1,en) + Rx ;

% y momentum equation - left
    [~,B,en] = get_data ( x_l , node , di , form) ;
    Ry = Dmat(3,:)*B/deltaX  ;
    K(2*ij  ,en) = K(2*ij  ,en) - Ry ;

% y momentum equation - right
    [~,B,en] = get_data ( x_r , node , di , form) ;
    Ry = Dmat(3,:)*B/deltaX  ;
    K(2*ij  ,en) = K(2*ij  ,en) + Ry ;

% y momentum equation - bot
    [~,B,en] = get_data ( x_b , node , di , form) ;
    Ry = Dmat(2,:)*B/deltaY  ;
    K(2*ij  ,en) = K(2*ij  ,en) - Ry ;

% y momentum equation - top
    [~,B,en] = get_data ( x_t , node , di , form) ;
    Ry = Dmat(2,:)*B/deltaY  ;
    K(2*ij  ,en) = K(2*ij  ,en) + Ry ;

end


if node_c(2) == 0 % bot edge
       plot(node_c(1),node_c(2),'bsq')
% only drichlet boundary conditions
        [phi,B,en] = get_data ( node_c , node , di , form) ;
        K(2*ij-1,en(1:2:end)) = K(2*ij-1,en(1:2:end)) + phi ;
        K(2*ij-0,en(2:2:end)) = K(2*ij-0,en(2:2:end)) + phi ;
        F(2*ij-1,1) = F(2*ij-1,1) + 0 ;
        F(2*ij-0,1) = F(2*ij-0,1) + 0 ;
elseif node_c(2) == L % top edge
       plot(node_c(1),node_c(2),'ksq')
% drichlet (ux=0.05) and neumann (t2=0) boundary conditions
        [phi,B,en] = get_data ( node_c , node , di , form) ;
        K(2*ij-0,en(2:2:end)) = K(2*ij-0,en(2:2:end)) + phi ;
        F(2*ij-0,1) = F(2*ij-0,1) + node_c(2)/L*0.05 ;
        
        K(2*ij-1,en) = K(2*ij-1,en) + Dmat(3,:)*B ;
        F(2*ij-1,1) = F(2*ij-1,1) + 0 ;
end


% Neuman boundary condition
if node_c(1) == 0 && node_c(2) >0 && node_c(2) < L && node_c(2) ~= D/2 % left edge
% neumann (t1=t2=0) boundary conditions
       plot(node_c(1),node_c(2),'gsq')
        [phi,B,en] = get_data ( node_c , node , di , form) ;
        K(2*ij-1,en) = K(2*ij-1,en) + Dmat(1,:)*B; %+ Dmat(2,:)*B;
        K(2*ij-0,en) = K(2*ij-0,en) + Dmat(3,:)*B;
        F(2*ij-1,1) = F(2*ij-1,1) + 0 ;
        F(2*ij-0,1) = F(2*ij-0,1) + 0 ;
elseif node_c(1) == L && node_c(2) >0 && node_c(2) < L && node_c(2) ~= D/2 % right edge
% neumann (t1=t2=0) boundary conditions
       plot(node_c(1),node_c(2),'gsq')
        [phi,B,en] = get_data ( node_c , node , di , form) ;
        K(2*ij-1,en) = K(2*ij-1,en) + Dmat(1,:)*B; % ; + Dmat(2,:)*B;
        K(2*ij-0,en) = K(2*ij-0,en) + Dmat(3,:)*B;
        F(2*ij-1,1) = F(2*ij-1,1) + 0 ;
        F(2*ij-0,1) = F(2*ij-0,1) + 0 ;
end

%% apply traction free boundary condition on the crack surface
if node_c(2) == 2 && node_c(1) >= 0 && node_c(1) < 120/100  % left edge
%     apply tx = ty = 0 
       plot(node_c(1),node_c(2),'bsq')
        [phi,B,en] = get_data ( node_c , node , di , form) ;
        K(2*ij-1,en) = K(2*ij-1,en) + Dmat(3,:)*B; % ; + Dmat(2,:)*B;
        K(2*ij-0,en) = K(2*ij-0,en) + Dmat(2,:)*B;
        F(2*ij-1,1) = F(2*ij-1,1) + 0 ;
        F(2*ij-0,1) = F(2*ij-0,1) + 0 ;
end

end


%%

figure
spy(K)


u = K\F ; 

fac = 1;
node_deformed = node ; 
node_deformed(:,1) = node(:,1) + fac*u(1:2:end);
node_deformed(:,2) = node(:,2) + fac*u(2:2:end);

figure
hold on
plot(node(:,1),node(:,2),'bsq')
plot(node_deformed(:,1),node_deformed(:,2),'rsq')
axis equal 



%% compute stress at nodes


% loop over nodes
    Stress = zeros(numnode,3); 
    Disp = zeros(numnode,2); 
for ij = 1 : size(node,1)
    node_c = node(ij,:) ; 
    [phi,B,en] = get_data ( node_c , node , di , form) ;
    Stress(ij,:) = Dmat * B * u(en) ;
     
    Disp(ij,:) = [phi * u(en(1:2:end)) phi * u(en(2:2:end))] ;

    
end



tri = delaunay(node(:,1),node(:,2));
% tri = tricheck(node(:,1:2),tri) ; 

VTKPostProcess(node,tri,1,'Tri3','stress',Stress,Disp)
!paraview stress.vtu&




