
clear all ; 
close all ;
! rm -r stress
! mkdir stress
! rm -r out
! mkdir out
! rm -r out/paraview
! mkdir out/paraview

dt            = 1e+5;
inc       = 0; % number of increment
incr0 = 2e-9*dt; % applied displacement at each increment

ninc          = 400;     % Number of increments
gitmax        = 50;     % Max. number of global iterations
tol_plast     = 1e-15;  % Tolerance of global iterations
nnx = 40 ;
Input_herr ; % input parameters and grid generation
LD = [] ; 

% initialize stresses at all nodes 
TSigma0 = 0 ; 
P0 = 0 ; 
C0 = coh0;
HH = zeros(size(node,1),(4+1+1)*4); % history parameter for L(1), B(2), R(3), T(4)
HH(:,[6 12 18 24]) = C0;  
HH_ = 0*HH ; 
HH_2 = zeros(size(node,1),4*4); % history parameter
HH_2_ = 0*HH_2 ; 

% get initial elastic solution matrix
stype = 'elastic'; dU_e = zeros(2*numnode,1); 


%% Newmark constants:
beta = 2;
gama = 1.5 ; 
theta = 1;
if (beta>=(0.25*(0.5+gama)^2))&&(theta>=0.5)&&(gama>=0.5)
    disp('Newmark Constants are OK!')
else
    error('Newmark Constants are not OK!')
end
a0 = 1/(beta*dt^2);
a1 = gama / (beta*dt);
a2 = 1 / (beta*dt) ; 
a3 = (gama/beta)-1 ;
a4 = (1/2/beta)-1 ; 
a5 = dt*((gama/2/beta)-1);

dynamic = 0 ; 
%%
mm_ = 0 ; 
for ij = 1 : numnode % loop over nodes 
    node_c = node(ij,:) ; 
    if (node_c(1) > 0 && node_c(1) < L && node_c(2) < D && node_c(2) > 0 )
        for ii = 1 : 4
        mm_ = mm_  + 1 ;
        end
    end
end

%%
dU_e = 0*dU_e ;
dU = dU_e;
dU0 = 0*dU_e ; 
plastic = 0 ; % set initial flag 
u = sparse(2*numnode,1);
au = u ;vu = u ;u0 = u ;au0 = u ;vu0 = u ;
h1_ = figure ; 
h2_ = figure ; 
while inc < 50000  % time stepping 
    inc = inc + 1 ; 

    disp(['step ---> ' num2str(inc)])
% Global Newton iterations : 
% Initialize the Nodal variables increment
%     dU( : ) = 0.0;

    for plast_it = 1 : gitmax  % until convergence
        SS = zeros(mm_,7) ; 
        DD =[ ] ; 
        cc_ij = 1 ;
        cc_ij2 = 1 ; 
%         SS = [] ; 
        K = sparse(2*numnode,2*numnode);
        Residual = zeros(2*numnode,1) ;  maxF = 0 ;  
        F = zeros(2*numnode,1) ;

        stepF = 1000 ;
        Cell = cell(stepF ,1);
        Cell2 = cell(stepF ,1);
        c_ = 1 ; 
        for ij = 1 : numnode % loop over nodes 
            node_c = node(ij,:) ; 
            [~,Dmat,~,~] = identify_tangent_v2 ( node_c , mat, dynamic , node) ;

            % only for nodes inside domain
            if (node_c(1) > 0 && node_c(1) < L && node_c(2) < D && node_c(2) > 0 )
                x_l = node_c - [deltaX/2 0 ] ;
                x_r = node_c + [deltaX/2 0 ] ;    

                x_b = node_c - [0 deltaY/2 ] ;
                x_t = node_c + [0 deltaY/2 ] ;    
%%
% knowing P0, TSigma0, C0 ---> P, TSigma, C, Dvp,
                side = 'L'; get_equations_vp ;
                side = 'R'; get_equations_vp ;
                side = 'B'; get_equations_vp ;
                side = 'T'; get_equations_vp ;
            end

            get_BCs_v3 () ;
        if rem(ij,stepF)==0
            IJV = cell2mat( Cell );
            A = sparse(IJV(:,1),IJV(:,2),IJV(:,3),size(K,1),size(K,1));
            K = K + A;
            Cell = cell(stepF,1);

            IJV = cell2mat( Cell2 );
            A = sparse(IJV(:,1),IJV(:,2),IJV(:,3),size(K,1),size(K,1));
            K = K + A;
            Cell2 = cell(stepF,1);            
            
            c_ = 1 ; 
        end
        end % on nodes 
        IJV = cell2mat( Cell );
        K = K + sparse(IJV(:,1),IJV(:,2),IJV(:,3),size(K,1),size(K,1));
        IJV = cell2mat( Cell2 );
        K = K + sparse(IJV(:,1),IJV(:,2),IJV(:,3),size(K,1),size(K,1));

%         pause
%         maxF
        if plastic == 1
%             disp('plasticity') 
        end
        
        nR = norm(F-Residual,1)/length(Residual); % residual norm
        disp(['iter: ' num2str(plast_it) '  norm: ' num2str(nR) ])

%         Residual = K*dU ;
%         norm(K*ddU-F)
        % check convergence 
%         if nR < tol_plast % converged 
%             break;
%         else % update displacements
            ddU    = K\(F-Residual) ;
            if any(isnan(ddU))
                error('ri')
            end
            dU     = dU + ddU;
%         end

    if nR < tol_plast % converged 
        break;
    end
    
au = a0*(dU - u0) - (a2*vu0) - (a4*au0) ;
vu = a1*(dU - u0) - (a3*vu0) - (a5*au0) ;

    end % on global Newton iterations 
    
    
u0 = dU ;
au0 = au ;
vu0 = vu ;
    
    LD = [LD ; inc*incr0 mean(SS(:,end))];
%             dU(find(node(:,1)==L)*2-1)
if rem(inc,1)==0
    fac = 100000000;
    node_deformed = node ; 
    node_deformed(:,1) = node(:,1) + fac*dU(1:2:end);
    node_deformed(:,2) = node(:,2) + fac*dU(2:2:end);
%     node_deformed(:,2) = node_deformed(:,2)-0.35; 
% if rem(inc,5)==0
    figure(h1_)
    clf
    hold on
    plot(node(:,1),node(:,2),'bsq')
    plot(node_deformed(:,1),node_deformed(:,2),'ksq')
    axis equal 
% end
% !paraview stress126.vtu&
    figure(h2_)
    clf
    hold on
    load LD_ref.mat 'LD_ref'
    plot(LD_ref(:,1),LD_ref(:,2),'rsq')
    plot(LD(:,1),LD(:,2),'bsq-')

    if exist('tri')==0
    tri = delaunay(SS(:,1:2));
    end
    VTKPostProcess(SS(:,1:2),tri,1,'Tri3',['stress/stress' num2str(inc)],SS(:,3:5),[SS(:,6) SS(:,7)])
    OutPut ;
    
    pause%(0.1)
end

    % update history parameters:
%     P0 = P; 
%     TSigma0 = TSigma ; 
%     C0        = C;
    HH = HH_;
    HH_2 = HH_2_ ; 
    dU0 = dU ;

end % on time stepping

% profile viewer
% profsave

