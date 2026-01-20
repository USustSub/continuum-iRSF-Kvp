
clear all ; 
close all ;
! rm -r stress
! mkdir stress

! rm -r out
! mkdir out
! rm -r out/slip
! mkdir out/slip
! rm -r out/mu
! mkdir out/mu
! rm -r out/deformed
! mkdir out/deformed
! rm -r out/paraview
! mkdir out/paraview

print_step1 = 10 ; 


% dt            = 0.0005;
dt            = 1e+7;
inc       = 0; % number of increment
incr0 = 2e-9*dt; % applied displacement at each increment
uy = 2e-9 ; 
ninc          = 10000;     % Number of increments
gitmax        = 50;     % Max. number of global iterations
tol_plast     = 1e-12;  % Tolerance of global iterations
nnx = 20 ;
dynamic = 1; 
Time = 0 ;

% Input_ ; % input parameters and grid generation
Input_herr
% profile on

% Gauss quadrature
order = 2 ;
[W,Q] = quadrature(order,'GAUSS',2);
nq = size(Q,1) ; 
tic

LD=[];

rightNodes = find ( abs ( node(:,1) - max(node(:,1)) ) < 0.0001 )  ;
leftNodes = find ( abs ( node(:,1) - min(node(:,1)) ) < 0.0001 )  ;
topNodes = find ( abs ( node(:,2) - max(node(:,2)) ) < 0.0001 )  ;
botNodes = find ( abs ( node(:,2) - min(node(:,2)) ) < 0.0001 )  ;

figure
hold on
plot(node(:,1),node(:,2),'bsq')
plot(node(fault_nodes,1),node(fault_nodes,2),'rsq')
plot_mesh(node,element(element_des,:),'Q4','k-');

plot(node(rightNodes,1),node(rightNodes,2),'sq')
plot(node(leftNodes,1),node(leftNodes,2),'sq')
plot(node(topNodes,1),node(topNodes,2),'sq')
plot(node(botNodes,1),node(botNodes,2),'sq')
axis equal

conn_nodes = find_conn_nodes (node, element);
%%
% initialize stresses at all nodes 
TSigma0 = 0 ; 
P0 = 1*5e6 ; 
C0 = coh0;
HH = zeros(size(element,1)*size(Q,1),4+1+1+4+2+4); % history parameter for L(1), B(2), R(3), T(4)
HH(:,6) = C0;
HH(:,5) = P0;
HH(:,11) = mu0;
HH(:,12) = theta_; 
HH_ = 0*HH ;     

%% Newmark constants:
% dt = 0.005 ;
beta = 2;
gama = 1.5 ; 
theta = 2;
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

if dynamic == 0
    a0 = 0 ; 
end
%%
dU_e = zeros(2*numnode,1); 
dU = dU_e;
dU0 = 0*dU_e ;
au0 = 0*dU_e ;
vu0 = 0*dU_e ;
plastic = 0 ; % set initial flag 
while inc < ninc  % time stepping 
    inc = inc + 1 ; 
    Time = Time + dt ; 
    disp(['Time is --> ' num2str(Time/60/60/24/365) ' years'])
     
    ap1 =  1*( uy + a3*vu0(topNodes(1)*2-1) + a5*au0(topNodes(1)*2-1) ) / a1 ;
    ap2 =  1*(-uy + a3*vu0(botNodes(1)*2-1) + a5*au0(botNodes(1)*2-1) ) / a1 ;
    incr0 = ap1 ; 
    
% Global Newton iterations : 
for plast_it = 1 : gitmax  % until convergence
    SS = zeros(size(element,1)*size(Q,1),13) ; 
    
    K = sparse(2*numnode,2*numnode);
    Mass = sparse(2*numnode,2*numnode);
    Residual = zeros(2*numnode,1) ;  maxF = -1e100 ;  
    volume_ = 0 ;  
    Cell = cell(size(element,1) ,1);
    Cell2 = cell(size(element,1) ,1);
    stepF = 1000 ; 
    maxVP = 0 ;
for iel = 1 : size(element,1) 
    sctr = element(iel,:); % element connectivity
    nn   = length(sctr);   % number of nodes per element

    % scatter vector for element assembly
    sctrB = zeros(1,length(sctr)*2) ;
    sctrB(1:2:end) = 2*sctr-1 ; 
    sctrB(2:2:end) = 2*sctr ;
    
    % ---------------------
    B = zeros(3,2*nn);
    Nu = zeros(2,2*nn);
    KK = zeros(8,8) ;
    NN = zeros(8,8) ;
    for kk = 1 : size(W,1)
        plastic = 0 ; % set initial flag 
        pt = Q(kk,:);                             % quadrature point
        % Shape functions and its derivatives
        [N,dNdxi] = lagrange_basis('Q4',pt);  % element shape functions
        J0 = node(sctr,:)'*dNdxi(:,:);                 % element Jacobian matrix
        xx = N(1:4)'*node(sctr,:) ;
        invJ0 = inv(J0);
        dNdx  = dNdxi*invJ0;                      % derivatives of N w.r.t XY
        
        % B matrix
        B(1,1:2:2*nn)      = dNdx(:,1)';
        B(2,2:2:2*nn) = dNdx(:,2)';
        B(3,1:2:2*nn)      = dNdx(:,2)';
        B(3,2:2:2*nn) = dNdx(:,1)';
        Nu(1,1:2:2*nn) = N; 
        Nu(2,2:2:2*nn) = N; 
        en = sctrB ; 
        uu = [N(1:4)'*dU(sctrB(1:2:end)) N(1:4)'*dU(sctrB(2:2:end)) 0];

%% get D,...
    TSigma0 = HH(iel*nq-(nq)+kk,1:4)';
    P0 = HH(iel*nq-(nq)+kk,5);
    P0 = 5e6 ;
    C0 = HH(iel*nq-(nq)+kk,6); 
    Eps = HH(iel*nq-(nq)+kk,7:10);
    mu_ = HH(iel*nq-(nq)+kk,11) ;
    theta_ = HH(iel*nq-(nq)+kk,12) ; 
    depsp_0 = HH(iel*nq-(nq)+kk,13:16) ;

    if element_id(iel)==0
        C0 = 10000000000; % elastic bulk 
        mu_ = 0.2 ; 
    else
        C0 = 0.0; % elastic bulk 
        mu_ = 0.05 ;         
    end
    % get strain
    dEps_t = B*(dU(en)-dU0(en)); 
    dEps_t = [dEps_t;0]; % add z strain
    dEps_0 = dEps_t; % save total strain  

    % get (visco)elastic tangent 
    [Dmat_t,Dmat,Gve,K_] = identify_tangent_v3 ( xx , mat , element_id(iel)) ;
    dSigma = Dmat_t * dEps_t ;% total stress increment 
    dP =-1/3*(sum(dSigma([1 2 4]))); % pressure increment
%       -1/3*(sum(dSigma([1 2 4])))
    % Updat total strains
    Eps_t = Eps' + [dEps_t(1) dEps_t(2) 0.5*dEps_t(3) dEps_t(4)]';
    
    % Deviatoric strain invariant
    Ekkc = Eps_t(1) + Eps_t(2) + Eps_t(4) ;
    Epsdc = Eps_t ; 
    Epsdc([1 2 4]) = Eps_t([1 2 4]) - 1/3*Ekkc;
    Eii  = sqrt(1/2*(Epsdc(1)).^2 + 1/2*(Epsdc(2)).^2 + 1/2*(Epsdc(4)).^2 + Epsdc(3).^2);

    Sigma_t = VE1* TSigma0 - [P0;P0;0;P0] + dSigma; % compute total trial stress (eq. 1)
    P   = P0 + dP; % get pressure
%     max(abs(dP))
%     dP
    TSigma = [P;P;0;P] + Sigma_t; % get deviatoric components
    J2 = 1/2*(TSigma(1)^2+TSigma(2)^2+TSigma(4)^2)+TSigma(3)^2;
    C     = C0; 
    F_     =  sqrt(J2) - C.*cos_phi - P.*mu_; % trial yield function
    % Plastic corrections
    dg = 0 ; dQdsigma = [0;0;0;0] ; 
    if F_>0 % if trial stresses exceed yield 
%         mu_
%         mu_
        plastic = 1;
        dQdsigma = get_dQdSigma ( TSigma , sin_psi , J2 );
        h1 = cos_phi.*h.*sqrt(2/3).*sqrt((dQdsigma(1).^2 + 2*dQdsigma(3).^2 + dQdsigma(2).^2 + dQdsigma(4).^2));
        dg = (F_./(Gve + K_.*mu_*sin_psi + eta_vp/dt + h1 ));% increment of plastic multiplier
        
        dEps_t = [dEps_0(1); dEps_0(2) ; dEps_0(3); dEps_0(4)]  - ...
            dg*[dQdsigma(1);dQdsigma(2);dQdsigma(3)/1;dQdsigma(4)]; % compute corrected strain
        dSigma = Dmat_t * dEps_t ;% total stress increment 
        dP    =-1/3*(sum(dSigma([1 2 4]))); % pressure increment
        Sigma_t = VE1 .* TSigma0 - [P0;P0;0;P0] + dSigma; % compute total stresses
        P   = P0 + dP; % get pressure
        TSigma = [P;P;0;P] + Sigma_t; % get deviatoric components

        % Check yield function (F_vp=0)
        J2 = 1/2*(TSigma(1).^2+TSigma(2).^2+TSigma(4).^2)+TSigma(3).^2;
        dep=sqrt(2/3).*sqrt((dg.*dQdsigma(1)).^2+(dg.*dQdsigma(2)).^2+(dg.*dQdsigma(4)).^2+2*(dg.*dQdsigma(3)).^2); % eq. 

        C     = C0 + h.*dep; 
        F_     =  sqrt(J2) - (C.*cos_phi + P.*mu_) ;                       % Backbone plastic yield function
        F_vp   =  sqrt(J2) - (C.*cos_phi + P.*mu_) - (dg./dt).*eta_vp;  % Consistency viscoplasticity model
%         fprintf('    Backbone EP max. Fc = %12.12e \n', max(F_(:)))
%         fprintf('    Viscoplast. max. Fc = %12.12e \n', max(F_vp(:)))

        % now get the tangent operator:
        dQdsigma = get_dQdSigma ( TSigma , sin_psi , J2 )';
        dFdsigma = get_dQdSigma ( TSigma , mu_ , J2 )';
        d2Qdsigma2 = get_dQ2dSigma2 ( TSigma  , J2 );
        [Dmat_t,~,Gve,K_] = identify_tangent_v3 ( xx , mat , element_id(iel)) ;

        I = [ 1 0 0 0 ; 0 1 0 0 ; 0 0 1 0 ; 0 0 0 1] ;    
        M = I + dg * Dmat_t * d2Qdsigma2 ;
        Mi = inv(M) ; 
        Dvp = Mi * Dmat_t -  ( Mi * Dmat_t * dQdsigma * dFdsigma' * Mi * Dmat_t ) / (  eta_vp/dt + h1 + dFdsigma' * Mi * Dmat_t * dQdsigma ) ; 
%                     
        Dmat = Dvp(1:3,1:3); % use viscoplastic tangent
%         maxFvp = max(abs(F_vp),maxFvp);
    end % F_ > 0
%}
        maxF = max((F_),maxF);
    
% %         Eps_p = Eps' + dg*dQdsigma;
    
    % Deviatoric strain invariant
% %         Ekkc = Eps_p(1) + Eps_p(2) + Eps_p(4) ;
% %         Epsdc = Eps_p ; 
% %         Epsdc([1 2 4]) = Eps_p([1 2 4]) - 1/3*Ekkc;
% %         Eii  = sqrt(1/2*(Epsdc(1)).^2 + 1/2*(Epsdc(2)).^2 + 1/2*(Epsdc(4)).^2 + Epsdc(3).^2);
% %         Eii = Eps_p(1) ; 
%         J2_all  = [J2_all ; J2 ] ; 

        % Stiffness matrix = B^T * D * B * detJ * w
        KK = KK + B'*Dmat*B*W(kk)*det(J0);
        NN = NN + Nu'*rho*Nu*W(kk)*det(J0);

        Residual(sctrB) = Residual(sctrB) + B'*[Sigma_t(1:3)]*W(kk)*det(J0);
        volume_ = volume_ +  W(kk)*det(J0) ; 

%% update friction coefficient mu_        
%         mu_ = mu0 ; 
%{
        depsp = 1*depsp_0' + dg*dQdsigma;
        Ekk = depsp(1) + depsp(2) + depsp(4) ;
%         if Ekk ~= 0
%             pause
%         end
        depsp_d = depsp ; 
        depsp_d([1 2 4]) = depsp_d([1 2 4]) - 1/3*Ekk;
        Eii_p  = sqrt(1/2*(depsp_d(1)).^2 + 1/2*(depsp_d(2)).^2 + 1/2*(depsp_d(4)).^2 + depsp_d(3).^2);
        Vp = 2*Eii_p/dt*deltaY ;

        if Vp~=0
            maxVP = max(abs(Vp),maxVP);
        end
        
        if xx(1) < LL1 || xx(1) > LL2
            b_ = b1_; 
        else
            b_ = b2_ ; 
        end
        
% update theta 
        S1 = 1+theta_/dt;
        S2 = (1/dt + Vp/L_);
        theta_n = S1/S2;
% update friction coefficient
        param = Vp/2/V0*exp((mu0+b_*log(theta_n*V0/L_))/a_);
        mu_gp = a_*asinh(param) ;
        if element_id(iel) ~= 0 && plast_it > 1 
            if mu_gp == 0 
                pause
            end
        end
%}
        mu_gp = mu0 ; Eii_p = 0 ; Vp = 0 ; theta_n = theta_ ;
%% store date         
%         SS = [SS ; xx Sigma_t' J2 log10(abs(Eii)) uu];
        SS(iel*nq-(nq)+kk,:) = [xx Sigma_t' mu_gp Eii_p uu Vp J2];
        HH_(iel*nq-(nq)+kk,1:4) = TSigma' ;
        HH_(iel*nq-(nq)+kk,5) =  P ;
        HH_(iel*nq-(nq)+kk,6) = C ; 
        HH_(iel*nq-(nq)+kk,7:10) = Eps_t; 
        HH_(iel*nq-(nq)+kk,11) = mu_gp; 
        HH_(iel*nq-(nq)+kk,12) = theta_n; 
%         HH_(iel*nq-(nq)+kk,13:16) = depsp'; 
        
    end                  % end of looping on GPs
%         theta_n
        [i j s] = find(KK); 
        Cell{iel} = [sctrB(i)' ,sctrB(j)', s];
    
        [i2 j2 s2] = find(NN); 
        Cell2{iel} = [sctrB(i2)' ,sctrB(j2)', s2];
        
        if rem(iel,stepF)==0
            IJV = cell2mat( Cell );
            A = sparse(IJV(:,1),IJV(:,2),IJV(:,3),size(K,1),size(K,1));
            K = K + A;
            Cell = cell(stepF,1);
            c = 0 ; 
        end
        
        if rem(iel,stepF)==0
            IJV = cell2mat( Cell2 );
            A = sparse(IJV(:,1),IJV(:,2),IJV(:,3),size(Mass,1),size(Mass,1));
            Mass = Mass + A;
            Cell2 = cell(stepF,1);
        end
end

IJV = cell2mat( Cell );
K = K + sparse(IJV(:,1),IJV(:,2),IJV(:,3),size(K,1),size(K,1));

IJV = cell2mat( Cell2 );
Mass = Mass + sparse(IJV(:,1),IJV(:,2),IJV(:,3),size(Mass,1),size(Mass,1));

%     Fint = Residual ;
% ************************

%     K = K + a0*Mass ;
    Converged_last_step = Mass*a0*dU0 + Mass*a2*vu0 + Mass*a4*au0 ;
    Residual = [Residual-dynamic*Converged_last_step] ;
    Residual = Residual + Mass*a0*dU; 

%     ApplyBCs ; 
    ApplyBCs_RSF_2 ;
%     nR = norm(Residual,1)/length(Residual); % residual norm
     
%     fprintf('inc %3d  itr %2d  Resisdual. norm = %2.2e \n', inc , plast_it ,nR)

%     maxVP

%{
% % % %     apply free slip boundary condition at left and right edges
% % %         ww = 1e6 ; 
% % %         Kp = sparse(length(K),length(K));
% % %         Rp = sparse(length(K),1);
% % %     for  ui = 1 : length(leftNodes) 
% % %         cur_node1 = leftNodes ( ui )  ;
% % %         iel = conn_nodes(cur_node1,2);
% % %         sctr =  element(iel,:) ;
% % %         sc_ = 2*sctr-1 ; 
% % %         local = element_coordinate ( node(cur_node1,:) , sctr , node );
% % %         [N,dNdxi] = lagrange_basis('Q4',local);  % element shape functions
% % %         J0 = node(sctr,:)'*dNdxi(:,:);                 % element Jacobian matrix
% % %         xx = N(1:4)'*node(sctr,:) ;
% % %         invJ0 = inv(J0);
% % %         dNdx  = dNdxi*invJ0;                      % derivatives of N w.r.t XY
% % %         scc = dNdx(:,1)' ; 
% % %         Kp(sc_,sc_) = Kp(sc_,sc_) + ww*scc'*scc;
% % %         Rp(sc_) = Rp(sc_) - ww*scc'*scc*dU(2*sctr-1) ;
% % %     end
% % %     for  ui = 1 : length(rightNodes) 
% % %         cur_node1 = rightNodes ( ui )  ;
% % %         iel = conn_nodes(cur_node1,2);
% % %         sctr =  element(iel,:) ;
% % %         sc_ = 2*sctr-1 ; 
% % %         local = element_coordinate ( node(cur_node1,:) , sctr , node );
% % %         [N,dNdxi] = lagrange_basis('Q4',local);  % element shape functions
% % %         J0 = node(sctr,:)'*dNdxi(:,:);                 % element Jacobian matrix
% % %         xx = N(1:4)'*node(sctr,:) ;
% % %         invJ0 = inv(J0);
% % %         dNdx  = dNdxi*invJ0;                      % derivatives of N w.r.t XY
% % %         scc = dNdx(:,1)' ; 
% % %         Kp(sc_,sc_) = Kp(sc_,sc_) + ww*scc'*scc;
% % %         Rp(sc_) = Rp(sc_) - ww*scc'*scc*dU(2*sctr-1) ;
% % %     end
% % %     if plast_it==1
% % %         K = K + Kp ; 
% % %         Residual = Residual + Rp ; 
% % %     end
%}
    
    
        % check convergence 
%     ddU    = K\(Residual) ;
%     dU     = dU + ddU;
        
    if nR < tol_plast % converged 
        break;
    end
end % on global Newton iterations 


    LD = [ LD ; Time/60/60/24/365 mean(sqrt(SS(:,13))) ];% abs(sum(Fint(2*botNodes-1)))  ] ;

    dlmwrite('LD.out',LD,'delimiter',' ')


%% calculate acceleration, velocity, ...
    au = a0*(dU - dU0) - (a2*vu0) - (a4*au0) ;
    vu = a1*(dU - dU0) - (a3*vu0) - (a5*au0) ;

        
% update history parameters:
    HH = HH_;
    dU0 = dU ;
    au0 = au ;
    vu0 = vu ;
%     max(HH(:,5)/5e6)
% postprocess    
    if rem(inc,1)==0
%         figure(h1_)
%         clf
%         load /home/mohsen/Downloads/m2di/M2Di_EP_KelvinViscoplastic/LD_TD.mat 'ff'
% %         plot(ff(:,1),ff(:,2),'bsq')
%         hold on
%         plot(LD(:,1),LD(:,2),'r-sq')
    
    end
    
if rem(inc , print_step1 ) == 0    

    
    fac = 10000 ; 
    node_deformed = node ; 
    node_deformed(:,1) = node(:,1) + fac*dU(1:2:end);
    node_deformed(:,2) = node(:,2) + fac*dU(2:2:end);
    figure('visible','off')
    hold on
    plot(node(:,1),node(:,2),'bsq')
    plot(node_deformed(:,1),node_deformed(:,2),'rsq')
    axis equal 
    %             plot_mesh(node_deformed,element(:,:),'Q4','k-');
    eval(['print -djpeg out/deformed/jj' num2str(inc) , '.jpeg'])

    
%     figure('visible','off')
%     plot(slip_(:,1),slip_(:,2),'rsq')
%     ylabel('slip')
%     eval(['print -djpeg out/slip/jj' num2str(itime) , '.jpeg'])

%     figure('visible','off')
%     plot(slip_rate(:,1),slip_rate(:,2),'rsq')
%     ylabel('slip_rate')
%     eval(['print -djpeg out/slip_rate/jj' num2str(itime) , '.jpeg'])

%     figure('visible','off')
%     plot(mu_(:,1),mu_(:,2),'rsq')
%     ylabel('mu')
%     eval(['print -djpeg out/mu/jj' num2str(itime) , '.jpeg'])
% 
%     figure('visible','off')
%     plot(theta_p(:,1),theta_p(:,2),'rsq')
%     ylabel('mu')
%     eval(['print -djpeg out/theta/jj' num2str(itime) , '.jpeg'])

    figure('visible','off')
    plot(LD(:,1),LD(:,2)/1e6,'r-sq','MarkerSize',1)
    ylabel('LD')
    hold on
    load /home/mohsen/Desktop/PhD/MyMatlabToolkit/Herr1.mat
    plot(Herr1(:,1),Herr1(:,2),'b-')
    ylabel('Stress[MPa]')
    eval(['print -djpeg out/LD.jpeg'])

    
    if exist('tri')==0    tri = delaunay(SS(:,1:2));    end
    VTKPostProcess(SS(:,1:2),tri,1,'Tri3',['out/paraview/jj' num2str(inc)],SS(:,3:5),[SS(:,7) SS(:,13)])
    
    
%     fac = 1000 ; 
%     DISP = u_elastic;
%     figure('visible','off')
%     hold on
%     DISP_ =  [DISP(1:2:end) DISP(2:2:end)];
%     node_deformed = node + fac * [DISP(1:2:end) DISP(2:2:end)] ; 
%     plot(node_deformed(:,1),node_deformed(:,2),'r.')
%     axis equal
%     eval(['print -djpeg out/deformed.jpeg'])
close all 

end
    
end % on time stepping
% profile viewer
% profsave

! paraview &

figure(h1_)
clf
load /home/mohsen/Downloads/m2di/M2Di_EP_KelvinViscoplastic/LD_TD.mat 'ff'
% plot(ff(:,1),ff(:,2),'bsq')
hold on
plot(LD(:,1),LD(:,2),'r-sq')

toc


%% RSF inputs
clear all ; close all ; 
mu0 = 0.2 ; 
V0 = 4e-9 ;  % m/s
a_ = 0.01 ;
L_ = 0.01 ;
b1_ = 0.001 ; 
b2_ = 0.017 ; 
theta_ = L_/V0*exp(1);
LL1 = 32000 ;
LL2 = 108000 ; 
dt  = 1e+4;
%
xx = [0:500:150000];
b_ = 0*xx ;
b_(:) = b2_;
b_(xx<LL1 | xx>LL2) = b1_;
h1 = figure ; 
Vp_ = [1e-14:(1e-9-1e-14)/10000:1e-9] ; 
for inc = 1 : 10000
    inc/10000
    if inc > 1 
        theta_ = theta_n;
    end
    Vp= Vp_(inc) ;
%     Vp = 1e-15 ; 
% (theta_n-theta_)/dt = 1 - (Vp*theta_n)/L_;
    S1 = 1+theta_/dt;
    S2 = (1/dt + Vp/L_);
    theta_n = S1/S2;
%     theta_n = L_/V0*exp(1) ; 
% mu = a_*asinh(Vp/2/V0*exp((mu0+b_.*log(theta_n*V0/L_))/a_));
    param = Vp/2/V0*exp((mu0+b_.*log(theta_n*V0/L_))/a_);
    mu_gp = a_*asinh(param) ;
    figure(h1)
    plot(xx,mu_gp,'r.')
%     ylim([0 1])
    pause(0.001)
end