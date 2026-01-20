%% 2D MLS based FDM 
% By Mohsen Goudarzi ( m.goudarzi@uu.nl ) 
% ----------------------------------------------------------
clear all; %close all ; 
clc
tic;           % help us to see the time required for each step

% Dimension of the domain (it is simply a rectangular region L x W)
L = 4 ;
D = 4 ;

% Material properties
E  = 100 ;
nu = 0.3 ;
BC = 2 ; 
incr0 = 5e-6; 

% Inputs special to meshless methods, domain of influence
shape = 'circle' ;         % shape of domain of influence
dmax  = 1.3 ;              % radius = dmax * nodal spacing
form  = 'cubic_spline' ;   % using cubic spline weight function


% Node density defined by number of nodes along two directions
nnx = 155 ;
nny = 155 ;

% node generation, node is of size (nnode x 2)
node = square_node_array([0 0],[L 0],[L D],[0 D],nnx,nny);
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

%

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


%%
Consti;

K = sparse(2*numnode,2*numnode);
F = zeros(2*numnode,1) ; 

% profile on
tic
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
    [~,B,en] = get_data ( x_r , node , di , form ) ;
    Rx = Dmat(1,:)*B/deltaX  ;
    K(2*ij-1,en) = K(2*ij-1,en) + Rx ;

% x momentum equation - bot
    [~,B,en] = get_data ( x_b , node , di , form ) ;
    Rx = Dmat(3,:)*B/deltaY  ;
    K(2*ij-1,en) = K(2*ij-1,en) - Rx ;
        
% x momentum equation - top
    [~,B,en] = get_data ( x_t , node , di , form ) ;
    Rx = Dmat(3,:)*B/deltaY  ;
    K(2*ij-1,en) = K(2*ij-1,en) + Rx ;

    
% y momentum equation - left
    [~,B,en] = get_data ( x_l , node , di , form ) ;
    Ry = Dmat(3,:)*B/deltaX  ;
    K(2*ij  ,en) = K(2*ij  ,en) - Ry ;

% y momentum equation - right
    [~,B,en] = get_data ( x_r , node , di , form ) ;
    Ry = Dmat(3,:)*B/deltaX  ;
    K(2*ij  ,en) = K(2*ij  ,en) + Ry ;
    
% y momentum equation - bot
    [~,B,en] = get_data ( x_b , node , di , form ) ;
    Ry = Dmat(2,:)*B/deltaY  ;
    K(2*ij  ,en) = K(2*ij  ,en) - Ry ;
        
% y momentum equation - top
    [~,B,en] = get_data ( x_t , node , di , form ) ;
    Ry = Dmat(2,:)*B/deltaY  ;
    K(2*ij  ,en) = K(2*ij  ,en) + Ry ;

end

if BC == 1

    if node_c(1) == 0 % left edge
    % only drichlet boundary conditions
            [phi,B,en] = get_data ( node_c , node , di , form ) ;
            K(2*ij-1,en(1:2:end)) = K(2*ij-1,en(1:2:end)) + phi ;
            K(2*ij-0,en(2:2:end)) = K(2*ij-0,en(2:2:end)) + phi ;
            F(2*ij-1,1) = F(2*ij-1,1) + 0 ;
            F(2*ij-0,1) = F(2*ij-0,1) + 0 ;
    elseif node_c(1) == L % right edge
    % drichlet (ux=0.05) and neumann (t2=0) boundary conditions
            [phi,B,en] = get_data ( node_c , node , di , form ) ;
            K(2*ij-1,en(1:2:end)) = K(2*ij-1,en(1:2:end)) + phi ;
            F(2*ij-1,1) = F(2*ij-1,1) + node_c(1)/L*0.05 ;

            K(2*ij-0,en) = K(2*ij-0,en) + Dmat(3,:)*B ;
            F(2*ij-0,1) = F(2*ij-0,1) + 0 ;
    end
    
% Neuman boundary condition
    if node_c(2) == 0 & ~( node_c(1) ==0 | node_c(1) == L ) % bot edge
    % neumann (t1=t2=0) boundary conditions
    %        plot(x(1),x(2),'gsq')
            [phi,B,en] = get_data ( node_c , node , di , form ) ;
            K(2*ij-1,en) = K(2*ij-1,en) + Dmat(3,:)*B; %+ Dmat(2,:)*B;
            K(2*ij-0,en) = K(2*ij-0,en) + Dmat(2,:)*B;
            F(2*ij-1,1) = F(2*ij-1,1) + 0 ;
            F(2*ij-0,1) = F(2*ij-0,1) + 0 ;
    elseif node_c(2) == L & ~( node_c(1) ==0 | node_c(1) == L ) % top edge
    % neumann (t1=t2=0) boundary conditions
            [phi,B,en] = get_data ( node_c , node , di , form ) ;
            K(2*ij-1,en) = K(2*ij-1,en) + Dmat(3,:)*B; % ; + Dmat(2,:)*B;
            K(2*ij-0,en) = K(2*ij-0,en) + Dmat(2,:)*B;
            F(2*ij-1,1) = F(2*ij-1,1) + 0 ;
            F(2*ij-0,1) = F(2*ij-0,1) + 0 ;
    end

elseif BC == 2 
            [phi,B,en] = get_data ( node_c , node , di , form ) ;
    if node_c(1) == 0 % left edge
% set ux = 0 
            K(2*ij-1,en(1:2:end)) = K(2*ij-1,en(1:2:end)) + phi ;
            F(2*ij-1,1) = F(2*ij-1,1) + 0 ;
% set uy = 0 
%             K(2*ij-0,en(2:2:end)) = K(2*ij-0,en(2:2:end)) + phi ;
%             F(2*ij-0,1) = F(2*ij-0,1) + 0 ;
% set duydx = 0 
            dphidx = B(1,1:2:end) ;
            K(2*ij-0,en(2:2:end)) = K(2*ij-0,en(2:2:end)) + dphidx ;
            F(2*ij-0,1) = F(2*ij-0,1) + 0 ;
    elseif node_c(1) == L % right edge
% set ux = -x*incr0 
            K(2*ij-1,en(1:2:end)) = K(2*ij-1,en(1:2:end)) + phi ;
            F(2*ij-1,1) = F(2*ij-1,1) - node_c(1)/L*0.05 ;
% set uy = y*incr0 
            K(2*ij-0,en(2:2:end)) = K(2*ij-0,en(2:2:end)) + phi ;
            F(2*ij-0,1) = F(2*ij-0,1) + node_c(2)/L*0.05 ;
    end

    if node_c(2) == 0% & ~( node_c(1) ==0 | node_c(1) == L ) % bot edge
% set uy = 0 
            K(2*ij-0,en(2:2:end)) = K(2*ij-0,en(2:2:end)) + phi ;
            F(2*ij-0,1) = F(2*ij-0,1) + 0 ;
% set duxdy = 0 
            dphidy = B(2,2:2:end) ;
            K(2*ij-1,en(1:2:end)) = K(2*ij-1,en(1:2:end)) + dphidy ;
            F(2*ij-1,1) = F(2*ij-1,1) + 0 ;
    elseif node_c(2) == L% & ~( node_c(1) ==0 | node_c(1) == L ) % top edge
% set ux = -x*incr0 
            K(2*ij-1,en(1:2:end)) = K(2*ij-1,en(1:2:end)) + phi ;
            F(2*ij-1,1) = F(2*ij-1,1) - node_c(1)/L*0.05 ;
% set uy = y*incr0 
            K(2*ij-0,en(2:2:end)) = K(2*ij-0,en(2:2:end)) + phi ;
            F(2*ij-0,1) = F(2*ij-0,1) + node_c(2)/L*0.05 ;
    end

end
    
end
toc
u = K\F ; 

fac = 20;
node_deformed = node ; 
node_deformed(:,1) = node(:,1) + fac*u(1:2:end);
node_deformed(:,2) = node(:,2) + fac*u(2:2:end);

figure
hold on
plot(node(:,1),node(:,2),'bsq')
plot(node_deformed(:,1),node_deformed(:,2),'rsq')
axis equal 

% profile viewer
% profsave


%% compute stress at nodes


% loop over nodes
    Stress = zeros(numnode,3); 
    Disp = zeros(numnode,2); 
for ij = 1 : size(node,1)
    node_c = node(ij,:) ; 
    [phi,B,en] = get_data ( node_c , node , di , form ) ;
    Stress(ij,:) = Dmat * B * u(en) ;
     
    Disp(ij,:) = [phi * u(en(1:2:end)) phi * u(en(2:2:end))] ;

    
end



% tri = delaunay(node(:,1),node(:,2));
% tri = tricheck(node(:,1:2),tri) ; 
tri = element ; 

VTKPostProcess(node,tri,1,'Quad4','stress',Stress,Disp)
!paraview stress.vtu&

return

%% Node density defined by number of nodes along two directions
nnx = 200 ;
nny = 200 ;

% node generation, node is of size (nnode x 2)
node_ = square_node_array([0 0],[L 0],[L D],[0 D],nnx,nny);
    SH = zeros(length(node_),4) ; 
for ij = 1 : length(node_)
    ij/length(node_)
    x = node_(ij,:); 
    [index] = define_support(node,x,di);
    if any(index==13)
        [phi,B,en] = get_data ( x , node , di , form ) ;
        sx = find(index==13);
        SH(ij,:) = [x phi(sx) B(1,2*sx-1)];
    end
end

SH(:,1:2) = node_ ; 
tri = delaunay(SH(:,1),SH(:,2));
% tri = tricheck(node(:,1:2),tri) ; 

VTKPostProcess(SH(:,1:2),tri,1,'Tri3','mls',SH(:,2:4),SH(:,3:4))
!paraview mls.vtu&



