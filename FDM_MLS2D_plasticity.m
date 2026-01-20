%% 2D MLS based FDM 
% By Mohsen Goudarzi ( m.goudarzi@uu.nl ) 
% ----------------------------------------------------------
clear all; 
close all ; 
clc
tic;           % help us to see the time required for each step

% Dimension of the domain (it is simply a rectangular region L x W)
L = 1 ;
D = 0.7 ;

% Input params
% Physics
eta0   = 2e50;        % Shear viscosity (DEACTIVATED)
etai   = 2e50;        % Inclusion shear viscosity (DEACTIVATED)
K0     = 2;           % Bulk modulus
G0     = 1;           % Shear modulus
Gi     = 1/4;         % Inclusion shear modulus
rad    = 0.05;        % Inclusion radius
incr0  = 5e-6;        % Incremental deformation
coh0   = 1.75e-4;     % Cohesion
phi0   = 30*pi/180;   % Friction angle
psi0   = 10*pi/180;   % Dilatancy angle
rho0   = 2700;        % Density
pconf  = 0;           % Confining pressure     
gy     = -0*9.81;     % Vertical gravity component
vm     = 1;           % factor 3 in von Mises - set to 1 for Drucker-Prager
eta_vp = 2.5e0;       % Kelvin element viscosity
h      = -0*0.01;     % Hardening/softening modulus
Cmin   = coh0/2;      % Minimum cohesion

% Input numerics
% Numerics
dt            = 1e4;
ninc          = 40;     % Number of increments
gitmax        = 50;     % Max. number of global iterations
tol_plast     = 1e-11;  % Tolerance of global iterations
nitPicNewt    = 0;      % Number of Picard steps before Newton
LineSearch    = 1;      % Activates line search
alpha         = 1.0;    % Default correction step
alpha_min     = 0.25;   % Minimum correction step
increase_step = 1;      % Allow for increasing increment
from_step     = 1160;   % Then set from_step = last recorded step
strain_max    = 4*5e-4; % Maximum strain before code stops running
imp_hard      = 1;      % Implicit formulation of hardening/softening
nout          = 1;      % visualises/outputs each nout step

% MLS_params
% Inputs special to meshless methods, domain of influence
shape = 'circle' ;         % shape of domain of influence
dmax  = 1.1 ;              % radius = dmax * nodal spacing
form  = 'cubic_spline' ;   % using cubic spline weight function


% Node density defined by number of nodes along two directions
nnx = 60 ;
nny = 60 ;

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


% Initialize
strain = 0; 
sin_phi = sin(phi0); 
cos_phi = cos(phi0); 
sin_psi = sin(psi0);
Kv   =   K0*ones(nnx ,nny ); % Elastic bulk modulus
Gv   =   G0*ones(nnx ,nny ); % Elastic shear modulus
etav = eta0*ones(nnx ,nny ); % Viscous shear modulus
Cv   = coh0*ones(nnx ,nny ); % cohesion
phic = phi0*ones(nnx,nny); % phi
phiv = phi0*ones(nnx ,nny ); % phi
psic = psi0*ones(nnx,nny); % psi
psiv = psi0*ones(nnx ,nny ); % psi
Kc   =   K0*ones(nnx,nny); % Elastic bulk modulus
Gc   =   G0*ones(nnx,nny); % Elastic shear modulus
etac = eta0*ones(nnx,nny); % Elastic shear modulus
Cc   = coh0*ones(nnx,nny); % cohesion
Sxxc =     zeros(nnx,nny); % Normal stress
Syyc =     zeros(nnx,nny); % Normal stress
Txyc =     zeros(nnx,nny); % Shear stress
Szzc =     zeros(nnx,nny); % Normal stress
Sxxv =     zeros(nnx ,nny ); % Normal stress
Syyv =     zeros(nnx ,nny ); % Normal stress
Txyv =     zeros(nnx ,nny ); % Shear stress
Szzv =     zeros(nnx ,nny ); % Normal stress
Exxc =     zeros(nnx,nny); % Normal strain
Eyyc =     zeros(nnx,nny); % Normal strain
Exyc =     zeros(nnx,nny); % Normal strain
Ezzc =     zeros(nnx,nny); % Normal strain
Exxv =     zeros(nnx ,nny ); % Normal strain
Eyyv =     zeros(nnx ,nny ); % Normal strain
Ezzv =     zeros(nnx ,nny ); % Normal strain
Ekkc =     zeros(nnx,nny); % Normal strain
Exyv =     zeros(nnx ,nny ); % Shear strain
Ux   =     zeros(nnx ,nny); % Horizontal displacement
Uy   =     zeros(nnx,nny ); % Vertical displacement
plc  =     zeros(nnx,nny); % Plastic flag centroids
plv  =     zeros(nnx ,nny ); % Plastic flag vertices
plc2 =     zeros(nnx,nny); % Plastic flag centroids
plv2 =     zeros(nnx ,nny ); % Plastic flag vertices
dgc  =     zeros(nnx,nny); % Plastic increment centroids
dgv  =     zeros(nnx ,nny ); % Plastic increment vertices
Eii  =     zeros(nnx,nny); % accumulated strain
Ep_acc_c = zeros(nnx,nny);  % accumulated strain
Ep_acc_v = zeros(nnx ,nny );  % accumulated strain
hc       = h*ones(nnx,nny); % hardening/softening modulus
hv       = h*ones(nnx ,nny ); % hardening/softening modulus
h1c      = h*ones(nnx,nny); % hardening/softening modulus
h1v      = h*ones(nnx ,nny ); % hardening/softening modulus
% For post-processing
Txxc   =    zeros(nnx,nny);
Tyyc   =    zeros(nnx,nny);
Tzzc   =    zeros(nnx,nny);
Txxc0  =    zeros(nnx,nny);
Tyyc0  =    zeros(nnx,nny);
Txyc0  =    zeros(nnx,nny);
Tzzc0  =    zeros(nnx,nny);
Txxv0  =    zeros(nnx ,nny );
Tyyv0  =    zeros(nnx ,nny );
Txyv0  =    zeros(nnx ,nny );
Tzzv0  =    zeros(nnx ,nny );
Txyv   =    zeros(nnx ,nny );
Pc0    =    zeros(nnx,nny);
Pv0    =    zeros(nnx ,nny );
% For monitoring
increment = zeros(ninc,1);
timevec   = zeros(ninc,1);
Pvec      = zeros(ninc,1);
Tiivec    = zeros(ninc,1);
epvec     = zeros(ninc,1);
strvec    = zeros(ninc,1);
nitervec  = zeros(ninc,1);
rxvec_rel = zeros(gitmax,1);    rxvec_abs = zeros(gitmax,1);
rvec_abs  = zeros(ninc,gitmax); rvec_rel  = zeros(ninc,gitmax);


% Set inclusion (shear modulus)
Gv((node(:,1)).^2+(node(:,2)-D).^2<rad^2) = Gi;
Gc = 0.25*(Gv(2:end,2:end) + Gv(1:end-1,2:end) + Gv(2:end,1:end-1) + Gv(1:end-1,1:end-1));
K  = 1e-6;
dt_diff   = dx^2/K/4.1/5;
diff_time = 15;
nt_diff   = ceil(diff_time/dt_diff);
dt_diff   = diff_time/nt_diff;
for i=1:nt_diff
    Gc_x = [Gc(1,:); Gc; Gc(end,:)];
    Gc_y = [Gc(:,1), Gc, Gc(:,end)];
    Gc   =  Gc + K*dt_diff/dx^2 * diff(Gc_x,2,1) + K*dt_diff/dy^2 * diff(Gc_y,2,2);
end
fprintf('diffusion of jump over a duration of %2.2f and using %02d steps\n', nt_diff*dt_diff, nt_diff);
[ Gv ] = M2Di_EP9_centroids2vertices( Gc );
% Elasto-viscous moduli
Gvev  = 1./(dt./etav + 1./Gv);
Gvec  = 1./(dt./etac + 1./Gc);
VEv   = Gvev./Gv;
VEc   = Gvec./Gc;
% Initial old stresses
time  = 0;
Exxc0 = Exxc; Exxv0 = Exxv;
Eyyc0 = Eyyc; Eyyv0 = Eyyv;
Exyv0 = Exyv; Exyc0 = Exyc;
Ezzc0 = Ezzc; Ezzv0 = Ezzv;
Ux0   = Ux;
Uy0   = Uy;

%%

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


% do plasticity
Get_stiffness_


% Initialize stiffness matrix for elasticity problem
% Assemble elastic matrix and RHS
[Ke, BcK] = M2Di_EP13_AssembleStiffness(plc,plv,plc2,plv2,Kc,Gvec,VEc,Cc,phic,psic,Kv,Gvev,VEv,phiv,psiv,dgc,dgv,Sxxc,Syyc,Txyc,Szzc,Sxxv,Syyv,Txyv,Szzv,dx,dy,nnx,nny,nnx,nny,NumUx,NumUy,NumUyG,Ux,Uy,0,vm, BC,sin_phi,sin_psi, cos_phi, coh0,symmetry,rho0,gy , -Pc0+Txxc0, -Pc0+Tyyc0, Txyv0, h, eta_vp, dt, 0, 0, 0, 0, 0 );

% get elastic solution
dU = Ke\BcK;

dUx_e    = dU(NumUx);
dUy_e    = dU(NumUyG);
ddU      = 0*dU;       % Initialise correction vector

%% ------------------------ Incremental loop ------------------------ %
inc       = 0;
plastic   = 0;
success   = 0;
finish    = 0;
incr0_ini = incr0;

while inc<ninc && finish == 0
    inc     = inc + 1;
    fail    = 0;
    new_BC  = 0;
    plc     =  zeros(nnx,nny); % Plastic flag centroids
    plv     =  zeros(nnx ,nny ); % Plastic flag vertices
    dgc     =  zeros(nnx,nny); % Plastic correction centroids
    dgv     =  zeros(nnx ,nny ); % Plastic correction vertices
    
%     initiaize dQds
    dQdsxxc = zeros(nnx,nny);
    dQdsyyc = zeros(nnx,nny);
    dQdszzc = zeros(nnx,nny);
    dQdsxyc = zeros(nnx,nny);
    dQdsxxv = zeros(nnx,nny);
    dQdsyyv = zeros(nnx,nny);
    dQdszzv = zeros(nnx,nny);
    dQdsxyv = zeros(nnx,nny);
    
%     initiaize gamdot
    gamdotc = zeros(nnx,nny);
    gamdotv = zeros(nnx,nny);
%     initiaize dep
    depc    = zeros(nnx,nny);
    depv    = zeros(nnx,nny);
  
    % Accumulated strain
    Ep_acc_c0  = Ep_acc_c;
    Ep_acc_v0  = Ep_acc_v;
    Cc0        = Cc;
    Cv0        = Cv;
    
    % Total Strain increment - the trial stress is elastic
    if plastic == 0
        fprintf('ELASTIC TRIAL\n');
        dUxc_t = dUx_e;
        dUyc_t = dUy_e;
    end
%
    fprintf('\n#######################################\n');
    fprintf(  '########## Loading step %04d ##########\n', inc);
    fprintf(  '#######################################\n');
%
    dU        = [dUxc_t(:); dUyc_t(:)];
    rxvec_rel = zeros(gitmax,1);    
    rxvec_abs = zeros(gitmax,1);    
    corvec    = zeros(gitmax,1);    
    alpvec    = zeros(gitmax,1);
    % ------------------------ Global equilibrium: Newton iterations ------------------------ %

    LSsuccess = 1; LIsuccessc=1; LIsuccessv = 1;
    for plast_it = 1:gitmax

        % Initial guess or iterative solution for strain increments
    
        
%  Start loop over nodes (TO DO)    
        
% get strain increments using dUxc_t (TO DO)
        dExxc_t  = diff(dUxc_t,1,1)/dx;
        dEyyc_t  = diff(dUyc_t,1,2)/dy;
        dEzzc_t  = 0; 
        dEzzv_t  = 0;
        dExyv_t  = 0.5*( dUxdy + dUydx );
        % Extrapolate trial strain components 
        dExyc_t     = 0.25*(dExyv_t(1:end-1,1:end-1) + dExyv_t(2:end,1:end-1) + dExyv_t(1:end-1,2:end) + dExyv_t(2:end,2:end));
        [ dExxv_t ] = M2Di_EP9_centroids2vertices( dExxc_t );
        [ dEyyv_t ] = M2Di_EP9_centroids2vertices( dEyyc_t );
        
        % Updat total strains (TO DO)
        Exxc_t = Exxc + dExxc_t;
        Eyyc_t = Eyyc + dEyyc_t;
        Exyv_t = Exyv + dExyv_t;
        Ezzc_t = Ezzc + 0;
        % Save strains 
        dExxc_0 = dExxc_t; dExxv_0 = dExxv_t; dEyyc_0 = dEyyc_t;  
        dEyyv_0 = dEyyv_t; dExyc_0 = dExyc_t; dExyv_0 = dExyv_t;
        dEzzc_0 = dEzzc_t; dEzzv_0 = dEzzv_t;
        % Total stress increments from total strain increment:
        dSxxc = (Kc + 4/3*Gvec).*dExxc_t + (Kc - 2/3*Gvec).*dEyyc_t + (Kc - 2/3*Gvec).*dEzzc_t;
        dSyyc = (Kc + 4/3*Gvec).*dEyyc_t + (Kc - 2/3*Gvec).*dExxc_t + (Kc - 2/3*Gvec).*dEzzc_t;
        dSzzc = (Kc + 4/3*Gvec).*dEzzc_t + (Kc - 2/3*Gvec).*dExxc_t + (Kc - 2/3*Gvec).*dEyyc_t;
        dTxyv = 2*Gvev.*dExyv_t ;
        dSxxv = (Kv + 4/3*Gvev).*dExxv_t + (Kv - 2/3*Gvev).*dEyyv_t + (Kv - 2/3*Gvev).*dEzzv_t;
        dSyyv = (Kv + 4/3*Gvev).*dEyyv_t + (Kv - 2/3*Gvev).*dExxv_t + (Kv - 2/3*Gvev).*dEzzv_t;
        dSzzv = (Kv + 4/3*Gvev).*dEzzv_t + (Kv - 2/3*Gvev).*dExxv_t + (Kv - 2/3*Gvev).*dEyyv_t;
        dTxyc = 2*Gvec.*dExyc_t;
        
        % Total stresses
        dPc    =-1/3*(dSxxc + dSyyc + dSzzc);
        dPv    =-1/3*(dSxxv + dSyyv + dSzzv);
        Sxxc_t = VEc .* Txxc0 - Pc0 + dSxxc;
        Syyc_t = VEc .* Tyyc0 - Pc0 + dSyyc;
        Txyv_t = VEv .* Txyv0       + dTxyv;
        Szzc_t = VEc .* Tzzc0 - Pc0 + dSzzc;
        Sxxv_t = VEv .* Txxv0 - Pv0 + dSxxv;
        Syyv_t = VEv .* Tyyv0 - Pv0 + dSyyv;
        Txyc_t = VEc .* Txyc0       + dTxyc;
        Szzv_t = VEv .* Tzzv0 - Pv0 + dSzzv;
        
        % Compute pressure and deviatoric components
        Pc   = Pc0 + dPc;
        Pv   = Pv0 + dPv;
        Txxc = Pc + Sxxc_t; 
        Tyyc = Pc + Syyc_t; 
        Tzzc = Pc + Szzc_t; 
        Txyc = Txyc_t;
        Txxv = Pv + Sxxv_t; 
        Tyyv = Pv + Syyv_t; 
        Tzzv = Pv + Szzv_t; 
        Txyv = Txyv_t;
        
        % Check yield function
        J2c    = 1/2*(Txxc.^2 + Tyyc.^2 + Tzzc.^2) + Txyc.^2;
        J2v    = 1/2*(Txxv.^2 + Tyyv.^2 + Tzzv.^2) + Txyv.^2;
        Cc     = Cc0; 
        Cv     = Cv0;
        syc    = Cc.*cos_phi + (Pc+pconf).*sin_phi;
        syv    = Cv.*cos_phi + (Pv+pconf).*sin_phi;
        Fc     =  sqrt(vm*J2c) - syc;
        Fv     =  sqrt(vm*J2v) - syv;
        fprintf('max. trial     Fc = %2.2e\n', max(Fc(:)))
        fprintf('max. trial     Fv = %2.2e\n', max(Fv(:)))
        % Plastic corrections
        if max(Fc(:))>0 || max(Fv(:))>0
            % Plastic flags
            plc = Fc>0;
            plv = Fv>0;
            % dQds centers
            txx = Pc+Sxxc_t; 
            tyy = Pc+Syyc_t;
            txy = Txyc_t; 
            tzz = Pc+Szzc_t; 
            J2 = J2c;
            dQdsxxc = vm * txx .* (vm * J2) .^ (-0.1e1 / 0.2e1) / 0.2e1 + sin_psi / 0.3e1;
            dQdsyyc = vm * tyy .* (vm * J2) .^ (-0.1e1 / 0.2e1) / 0.2e1 + sin_psi / 0.3e1;
            dQdsxyc = vm * txy .* (vm * J2) .^ (-0.1e1 / 0.2e1);
            dQdszzc = vm * tzz .* (vm * J2) .^ (-0.1e1 / 0.2e1) / 0.2e1 + sin_psi / 0.3e1;
            % dQds vertices
            txx = Pv+Sxxv_t; tyy = Pv+Syyv_t; txy = Txyv_t; tzz = Pv+Szzv_t; J2 = J2v;
            dQdsxxv = vm * txx .* (vm * J2) .^ (-0.1e1 / 0.2e1) / 0.2e1 + sin_psi / 0.3e1;
            dQdsyyv = vm * tyy .* (vm * J2) .^ (-0.1e1 / 0.2e1) / 0.2e1 + sin_psi / 0.3e1;
            dQdsxyv = 0.1000000000e1 * vm * txy .* (vm * J2) .^ (-0.1e1 / 0.2e1);
            dQdszzv = vm * tzz .* (vm * J2) .^ (-0.1e1 / 0.2e1) / 0.2e1 + sin_psi / 0.3e1;
            % dgamma - increment of plastic multiplier
            h1c = cos_phi.*hc.*sqrt(2/3).*sqrt((dQdsxxc.^2 + 2*dQdsxyc.^2 + dQdsyyc.^2 + dQdszzc.^2));
            h1v = cos_phi.*hv.*sqrt(2/3).*sqrt((dQdsxxv.^2 + 2*dQdsxyv.^2 + dQdsyyv.^2 + dQdszzv.^2));
            dgc = plc.*(Fc./(vm.*Gvec + Kc.*sin_phi*sin_psi + eta_vp/dt + h1c ));
            dgv = plv.*(Fv./(vm.*Gvev + Kv.*sin_phi*sin_psi + eta_vp/dt + h1v ));
            % Corrected strain on centers
            dExxc_t = dExxc_0 - plc.*dgc.*dQdsxxc;
            dEyyc_t = dEyyc_0 - plc.*dgc.*dQdsyyc;
            dExyc_t = dExyc_0 - plc.*dgc.*dQdsxyc/2;
            dEzzc_t = dEzzc_0 - plc.*dgc.*dQdszzc;
            % Corrected strain on vertices
            dExxv_t = dExxv_0 - plv.*dgv.*dQdsxxv;
            dEyyv_t = dEyyv_0 - plv.*dgv.*dQdsyyv;
            dExyv_t = dExyv_0 - plv.*dgv.*dQdsxyv/2;
            dEzzv_t = dEzzv_0 - plv.*dgv.*dQdszzv;
            % Total stress increments
            dSxxc = (Kc + 4/3*Gvec).*dExxc_t + (Kc - 2/3*Gvec).*dEyyc_t + (Kc - 2/3*Gvec).*dEzzc_t;
            dSyyc = (Kc + 4/3*Gvec).*dEyyc_t + (Kc - 2/3*Gvec).*dExxc_t + (Kc - 2/3*Gvec).*dEzzc_t;
            dSzzc = (Kc + 4/3*Gvec).*dEzzc_t + (Kc - 2/3*Gvec).*dExxc_t + (Kc - 2/3*Gvec).*dEyyc_t;
            dTxyv = 2*Gvev.*dExyv_t ;
            dSxxv = (Kv + 4/3*Gvev).*dExxv_t + (Kv - 2/3*Gvev).*dEyyv_t + (Kv - 2/3*Gvev).*dEzzv_t;
            dSyyv = (Kv + 4/3*Gvev).*dEyyv_t + (Kv - 2/3*Gvev).*dExxv_t + (Kv - 2/3*Gvev).*dEzzv_t;
            dSzzv = (Kv + 4/3*Gvev).*dEzzv_t + (Kv - 2/3*Gvev).*dExxv_t + (Kv - 2/3*Gvev).*dEyyv_t;
            dTxyc = 2*Gvec.*dExyc_t;
            % Total stresses
            dPc    =-1/3*(dSxxc + dSyyc + dSzzc);
            dPv    =-1/3*(dSxxv + dSyyv + dSzzv);
            Sxxc_t = VEc .* Txxc0 - Pc0 + dSxxc;
            Syyc_t = VEc .* Tyyc0 - Pc0 + dSyyc;
            Txyv_t = VEv .* Txyv0       + dTxyv;
            Szzc_t = VEc .* Tzzc0 - Pc0 + dSzzc;
            Sxxv_t = VEv .* Txxv0 - Pv0 + dSxxv;
            Syyv_t = VEv .* Tyyv0 - Pv0 + dSyyv;
            Txyc_t = VEc .* Txyc0       + dTxyc;
            Szzv_t = VEv .* Tzzv0 - Pv0 + dSzzv;
            % Compute pressure and deviatoric components
            Pc     = Pc0 + dPc;
            Pv     = Pv0 + dPv;
            Txxc   = Pc + Sxxc_t; Tyyc = Pc + Syyc_t; Tzzc = Pc + Szzc_t; Txyc = Txyc_t;
            Txxv   = Pv + Sxxv_t; Tyyv = Pv + Syyv_t; Tzzv = Pv + Szzv_t; Txyv = Txyv_t;
            % Check yield function
            J2c    = 1/2*(Txxc.^2 + Tyyc.^2 + Tzzc.^2) + Txyc.^2;
            J2v    = 1/2*(Txxv.^2 + Tyyv.^2 + Tzzv.^2) + Txyv.^2;
            depc   = sqrt(2/3).*sqrt((dgc.*dQdsxxc).^2+(dgc.*dQdsyyc).^2+(dgc.*dQdszzc).^2+2*(dgc.*dQdsxyc).^2);
            depv   = sqrt(2/3).*sqrt((dgv.*dQdsxxv).^2+(dgv.*dQdsyyv).^2+(dgv.*dQdszzv).^2+2*(dgv.*dQdsxyv).^2);
            Cc     = Cc0 + imp_hard*hc.*depc;
            Cv     = Cv0 + imp_hard*hv.*depv;
            syc    = Cc.*cos_phi + (Pc+pconf).*sin_phi;
            syv    = Cv.*cos_phi + (Pv+pconf).*sin_phi;
            Fc     =  sqrt(vm*J2c) - syc;                       % Backbone plastic yield function
            Fv     =  sqrt(vm*J2v) - syv;
            gamdotc = plc.*(sqrt(vm*J2c) - syc) ./ eta_vp; % plastic strain rate
            gamdotv = plv.*(sqrt(vm*J2v) - syv) ./ eta_vp; % plastic strain rate
            %             Fc_vp   =  sqrt(vm*J2c) - syc - gamdotc.*eta_vp; % equivalent to line below
            %             Fv_vp   =  sqrt(vm*J2v) - syv - gamdotv.*eta_vp;
            Fc_vp   =  sqrt(vm*J2c) - syc - (dgc./dt).*eta_vp;  % Consistency viscoplasticity model
            Fv_vp   =  sqrt(vm*J2v) - syv - (dgv./dt).*eta_vp;
            fprintf('    Backbone EP max. Fc = %2.2e \n', max(Fc(:)))
            fprintf('    Backbone EP max. Fv = %2.2e \n', max(Fv(:)))
            fprintf('    Viscoplast. max. Fc = %2.2e \n', max(Fc_vp(:)))
            fprintf('    Viscoplast. max. Fv = %2.2e \n', max(Fv_vp(:)))
        end
        % Evaluate global residual
        Rx = [zeros(1,nny); (diff(Sxxc_t,1,1)/dx + diff(Txyv_t(2:end-1,:),1,2)/dy)          ; zeros(1,nny)];
        Ry = [zeros(nnx,1), (diff(Syyc_t,1,2)/dy + diff(Txyv_t(:,2:end-1),1,1)/dx) + rho0*gy, zeros(nnx,1)];
        R  = [Rx(:); Ry(:)];
        nR = norm(R,1)/length(R);
        if plast_it==1, nR0 = nR; end
        fprintf('Iteration %03d ||F|| = %2.8e - ||F0|| = %2.8e - ||F/F0|| = %2.8e - Newton = %d\n', plast_it, nR, nR0, nR/nR0, Newton);
        rxvec_abs(plast_it) = nR;
        rxvec_rel(plast_it) = nR/nR0;
        % Exit if return mapping failed
        if LIsuccessc==0 || LIsuccessv==0 || nR/nR0>1e3, break; end
        % Exit if convergence criteria is met
        if nR<tol_plast
            break;
        else
            % Identify plastic nodes
            plastic = 1;
            % Assemble elasto-plastic tangent matrix
                [Kep,~] = M2Di_EP13_AssembleStiffness(plc,plv,plc2,plv2,Kc,Gvec,VEc,Cc,phic,psic,Kv,Gvev,VEv,phiv,psiv,dgc,dgv,Sxxc_t,Syyc_t,Txyc_t,Szzc_t,Sxxv_t,Syyv_t,Txyv_t,Szzv_t,dx,dy,nnx,nny,nnx,nny,NumUx,NumUy,NumUyG,Ux,Uy,Newton,vm, BC,sin_phi,sin_psi, cos_phi, coh0,symmetry,rho0,gy , -Pc0+Txxc0, -Pc0+Tyyc0, Txyv0, h, eta_vp, dt, syc, syv, 0, h1c, h1v );
                ddU       = Kep\R;
                fprintf('BACKSLASH solve took = %2.2e s\n', toc)
            % Line search algorithm
            % Update incremental displcaments 
            if LSsuccess == 1
                dU     = dU + alpha*ddU;
                dUxc_t = dU(NumUx);
                dUyc_t = dU(NumUyG);
                nC     = norm(ddU)/length(ddU);
            else
                % If line search failed: exit
                plast_it = gitmax;
                break;
            end
        end
        alpvec(plast_it) = alpha;
        corvec(plast_it) = nC;
    end
    if  nR<tol_plast
        % Update total strains and displacements
        Exxc = Exxc_t;
        Eyyc = Eyyc_t;
        Ezzc = Ezzc_t;
        Exyv = Exyv_t;
        Ekkc = Exxc_t + Eyyc_t + Ezzc_t;
        dExxdc_0 = dExxc_0 - 1/3*(dExxc_0+dEyyc_0);
        dEyydc_0 = dEyyc_0 - 1/3*(dExxc_0+dEyyc_0);
        dEzzdc_0 = dEzzc_0 - 1/3*(dExxc_0+dEyyc_0);
        dExydv_0 = dExyv_0;
        Ux   = Ux0 + dUxc_t;
        Uy   = Uy0 + dUyc_t;
        Ux0  = Ux;
        Uy0  = Uy;
        Ep_acc_c = Ep_acc_c0 + depc;
        Ep_acc_v = Ep_acc_v0 + depv;
        % Softening (explicit form)
        if imp_hard==0
            Cc     = Cc0 + hc.*depc ;
            Cv     = Cv0 + hv.*depv ;
        end
        % limit softening
        hc(Cc<=Cmin) = 0;
        hv(Cv<=Cmin) = 0;
        Cc(Cc<=Cmin) = Cmin;
        Cv(Cv<=Cmin) = Cmin;
        strain          = strain + incr0;
        time            = time + dt;
        timevec(inc)    = time;
        increment(inc)  = incr0;
        Tiivec(inc)     = mean(sqrt(vm*J2c(:)));
        strvec(inc)     = strain;
        rvec_rel(inc,:) = rxvec_rel(:);
        rvec_abs(inc,:) = rxvec_abs(:);
        nitervec(inc)   = plast_it;
        Pvec(inc)       = Pc(fix(nnx/2), fix(nny/2));
        epvec(inc)      = norm( Ep_acc_c );
        % Post-processing
        % Deviatoric strain invariant
        Exxdc = Exxc - 1/3*Ekkc;
        Eyydc = Eyyc - 1/3*Ekkc;
        Ezzdc = Ezzc - 1/3*Ekkc;
        Eii   = M2Di_EP10_InvariantOnCentroids(Exxdc   , Eyydc,   Ezzdc,   Exyv   );
        % Deviatoric strain rate invariant
        Exxcr = (Exxc-Exxc0)/dt; Eyycr = (Eyyc-Eyyc0)/dt; Ezzcr = (Ezzc-Ezzc0)/dt; Exyvr = (Exyv-Exyv0)/dt;
        Ekkcr = Exxcr + Eyycr + Ezzcr;
        Exxdr = Exxcr - 1/3*Ekkcr;
        Eyydr = Eyycr - 1/3*Ekkcr;
        Ezzdr = Ezzcr - 1/3*Ekkcr;
        Eiir  = M2Di_EP10_InvariantOnCentroids(Exxdr   , Eyydr,   Ezzdr,   Exyvr   );
        % Deviatoric strain rates decomposed
        Ekkp    = dgc.*(dQdsxxc + dQdsyyc + dQdszzc);
        Exxcr   = dExxdc_0/dt;                     Eyycr   = dEyydc_0/dt;                     Ezzcr   = dEzzdc_0/dt;                     Exyvr   = dExyv_0/dt;
        Exxcr_e = (Txxc-Txxc0)./2./Gc/dt;          Eyycr_e = (Tyyc-Tyyc0)./2./Gc/dt;          Ezzcr_e = (Tzzc-Tzzc0)./2./Gc/dt;          Exyvr_e = (Txyv-Txyv0)./2./Gv/dt;
        Exxcr_v = Txxc./2./etac;                   Eyycr_v = Tyyc./2./etac;                   Ezzcr_v = Tzzc./2./etac;                   Exyvr_v = Txyv./2./etav;
        Exxcr_p = 1/dt*(dgc.*dQdsxxc - 1/3.*Ekkp); Eyycr_p = 1/dt*(dgc.*dQdsyyc - 1/3.*Ekkp); Ezzcr_p = 1/dt*(dgc.*dQdszzc - 1/3.*Ekkp); Exyvr_p = 1/dt*dgv.*dQdsxyv / 2;
        Eiicr   = M2Di_EP10_InvariantOnCentroids(Exxcr   , Eyycr,   Ezzcr,   Exyvr   );
        Eiicr_e = M2Di_EP10_InvariantOnCentroids(Exxcr_e , Eyycr_e, Ezzcr_e, Exyvr_e );
        Eiicr_v = M2Di_EP10_InvariantOnCentroids(Exxcr_v , Eyycr_v, Ezzcr_v, Exyvr_v );
        Eiicr_p = M2Di_EP10_InvariantOnCentroids(Exxcr_p , Eyycr_p, Ezzcr_p, Exyvr_p );
        Exx_net = Exxcr-Exxcr_e-Exxcr_v-Exxcr_p;
        Eyy_net = Eyycr-Eyycr_e-Eyycr_v-Eyycr_p;
        Ezz_net = Ezzcr-Ezzcr_e-Ezzcr_v-Ezzcr_p;
        Exy_net = Exyvr-Exyvr_e-Exyvr_v-Exyvr_p;
        Eii_net = M2Di_EP10_InvariantOnCentroids(Exx_net   , Eyy_net,   Ezz_net,   Exy_net   );
        % Update old stresses
        Pc0 = Pc; Pv0 = Pv;
        Sxxc0 = Sxxc_t; Sxxv0 = Sxxv_t;
        Syyc0 = Syyc_t; Syyv0 = Syyv_t;
        Txyv0 = Txyv_t; Txyc0 = Txyc_t;
        Szzc0 = Szzc_t; Szzv0 = Szzv_t;
        Txxc0 = Txxc  ; Tyyc0 = Tyyc  ; Tzzc0 = Tzzc;
        Txxv0 = Txxv  ; Tyyv0 = Tyyv  ; Tzzv0 = Tzzv;
        Exxdc0 = Exxdc;Eyydc0 = Eyydc;Ezzdc0 = Ezzdc; Exyv0 = Exyv;
        last_num_its = plast_it;
        success = 1;
    else
        % Then try with smaller increment
        success = 0;
        if (nR/nR0)>100; success=0; end
        inc = inc-1;
    end
    % Visualize (only after elasto-plastic iterations)
% %     if plast_it>1 && success==1
% %         if mod(inc,nout)==0
% %             figure(3), clf
% %             subplot(311), hold on
% %             imagesc(xc, yc, log10(abs(Eii'))), colorbar, axis image, set(gca, 'ydir', 'normal'), title(['min. Eii = ', num2str(min(Eii(:))), ' max. Eii = ', num2str(max(Eii(:)))]), colormap('jet')
% %             subplot(312), hold on             
% %             plot(strvec(1:inc),Tiivec(1:inc), 'dr'), title('Stress vs. strain')
% %             subplot(313), hold on
% %             plot(1:plast_it, log10(rxvec_abs(1:plast_it)),'ok')
% %             title('Global residual vs. iterations')
% %             drawnow
% % %             pause
% %             % Save
% %             if saveRunData==1 && mod(inc,nout)==0,
% %                 print([path,'/TimeEvol',num2str(inc,'%04d')], '-r300','-dpng');
% %             end 
% %         end
% %         % Exit increments if strain_max is exceeded
% %         if  strain> strain_max, finish = 1; end
% %     end
end



%% old mfree 
Consti;

K = sparse(2*numnode,2*numnode);
F = zeros(2*numnode,1) ; 

% loop over nodes and assemble stiffness matrix
for ij = 1 : size(node,1)
    node_c = node(ij,:) ; 
    
if (node_c(1) > 0 && node_c(1) < L && node_c(2) < D && node_c(2) > 0 )
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
elseif node_c(2) == D & ~( node_c(1) ==0 | node_c(1) == L ) % top edge
% neumann (t1=t2=0) boundary conditions
        [phi,B,en] = get_data ( node_c , node , di , form ) ;
        K(2*ij-1,en) = K(2*ij-1,en) + Dmat(3,:)*B; % ; + Dmat(2,:)*B;
        K(2*ij-0,en) = K(2*ij-0,en) + Dmat(2,:)*B;
        F(2*ij-1,1) = F(2*ij-1,1) + 0 ;
        F(2*ij-0,1) = F(2*ij-0,1) + 0 ;
end

    
end

u = K\F ; 

fac = 10;
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
    [phi,B,en] = get_data ( node_c , node , di , form ) ;
    Stress(ij,:) = Dmat * B * u(en) ;
     
    Disp(ij,:) = [phi * u(en(1:2:end)) phi * u(en(2:2:end))] ;

    
end



% tri = delaunay(node(:,1),node(:,2));
% tri = tricheck(node(:,1:2),tri) ; 
tri = element ; 

VTKPostProcess(node,tri,1,'Quad4','stress',Stress,Disp)
!paraview stress.vtu&




