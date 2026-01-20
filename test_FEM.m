
clear all ; 
close all ;
! rm -r stress
! mkdir stress
Ex = 'Duretz' ; 
addpath '../Test2_'
% dt            = 0.0005;
dt            = 10^4;
inc       = 0; % number of increment
incr0 = 0.000005; % applied displacement at each increment
ninc          = 30;     % Number of increments
gitmax        = 50;     % Max. number of global iterations
tol_plast     = 1e-13;  % Tolerance of global iterations
nnx = 20 ;

Input_ ; % input parameters and grid generation
% figure
% plot_mesh(node,element(:,:),'Q4','k-');

% profile on

% Gauss quadrature
order = 2 ;
[W,Q] = quadrature(order,'GAUSS',2);
nq = size(Q,1) ; 
tic

LD=[];
%%
% initialize stresses at all nodes 
TSigma0 = 0 ; 
P0 = 0 ; 
C0 = coh0;
HH = zeros(size(element,1)*size(Q,1),4+1+1+4); % history parameter for L(1), B(2), R(3), T(4)
HH(:,6) = C0;  
HH_ = 0*HH ; 

%%
h1_ = figure ; 
dU_e = zeros(2*numnode,1); 
dU = dU_e;
dU0 = 0*dU_e ; 
plastic = 0 ; % set initial flag 
while inc < ninc  % time stepping 
    inc = inc + 1 ; 

% Global Newton iterations : 
for plast_it = 1 : gitmax  % until convergence
    SS = zeros(size(element,1)*size(Q,1),11) ; 
    
    K = sparse(2*numnode,2*numnode);
    Residual = zeros(2*numnode,1) ;  maxF = 0 ;  
    volume_ = 0 ;  
    Cell = cell(size(element,1) ,1);
    stepF = 1000 ;
for iel = 1 : size(element,1) 
    sctr = element(iel,:); % element connectivity
    nn   = length(sctr);   % number of nodes per element

    % scatter vector for element assembly
    sctrB = zeros(1,length(sctr)*2) ;
    sctrB(1:2:end) = 2*sctr-1 ; 
    sctrB(2:2:end) = 2*sctr ;
    
    % ---------------------
    B = zeros(3,2*nn);
    KK = zeros(8,8) ;
    for kk = 1 : size(W,1)
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
        en = sctrB ; 
        uu = [N(1:4)'*dU(sctrB(1:2:end)) N(1:4)'*dU(sctrB(2:2:end)) 0];

%% get D,...
    TSigma0 = HH(iel*nq-(nq)+kk,1:4)';
    P0 = HH(iel*nq-(nq)+kk,5);
    C0 = HH(iel*nq-(nq)+kk,6); 
    Eps = HH(iel*nq-(nq)+kk,7:10); 
    
    % get strain
    dEps_t = B*(dU(en)-dU0(en)); 
    dEps_t = [dEps_t;0]; % add z strain
    dEps_0 = dEps_t; % save total strain  

    % get (visco)elastic tangent 
    [Dmat_t,Dmat,Gve,K_] = identify_tangent_v2 ( xx , mat) ;
    dSigma = Dmat_t * dEps_t ;% total stress increment 
    dP =-1/3*(sum(dSigma([1 2 4]))); % pressure increment
    if plast_it > 1
%       dP
%       pause
    end
%    dSigma
%       pause(0.1)
    % Updat total strains
    Eps_t = Eps' + [dEps_t(1) dEps_t(2) 0.5*dEps_t(3) dEps_t(4)]';
    
    % Deviatoric strain invariant
    Ekkc = Eps_t(1) + Eps_t(2) + Eps_t(4) ;
    Epsdc = Eps_t ; 
    Epsdc([1 2 4]) = Eps_t([1 2 4]) - 1/3*Ekkc;
    Eii  = sqrt(1/2*(Epsdc(1)).^2 + 1/2*(Epsdc(2)).^2 + 1/2*(Epsdc(4)).^2 + Epsdc(3).^2);

    Sigma_t = VE1* TSigma0 - [P0;P0;0;P0] + dSigma; % compute total trial stress (eq. 1)
    P   = P0 + dP; % get pressure
    TSigma = [P;P;0;P] + Sigma_t; % get deviatoric components
    J2 = 1/2*(TSigma(1)^2+TSigma(2)^2+TSigma(4)^2)+TSigma(3)^2;
    C     = C0; 
    F_     =  sqrt(J2) - C.*cos_phi - P.*sin_phi; % trial yield function
    % Plastic corrections
    dg = 0 ; dQdsigma = [0;0;0;0] ; 
    if F_>0 % if trial stresses exceed yield 
        plastic = 1;
        dQdsigma = get_dQdSigma ( TSigma , sin_psi , J2 );
        h1 = cos_phi.*h.*sqrt(2/3).*sqrt((dQdsigma(1).^2 + 2*dQdsigma(3).^2 + dQdsigma(2).^2 + dQdsigma(4).^2));
        dg = (F_./(Gve + K_.*sin_phi*sin_psi + eta_vp/dt + h1 ));% increment of plastic multiplier
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
        F_     =  sqrt(J2) - (C.*cos_phi + P.*sin_phi) ;                       % Backbone plastic yield function
        F_vp   =  sqrt(J2) - (C.*cos_phi + P.*sin_phi) - (dg./dt).*eta_vp;  % Consistency viscoplasticity model
% %         fprintf('    Backbone EP max. Fc = %2.2e \n', max(F_(:)))
% %         fprintf('    Viscoplast. max. Fc = %2.2e \n', max(F_vp(:)))


        % now get the tangent operator:
        dQdsigma = get_dQdSigma ( TSigma , sin_psi , J2 )';
        dFdsigma = get_dQdSigma ( TSigma , sin_phi , J2 )';
        d2Qdsigma2 = get_dQ2dSigma2 ( TSigma  , J2 );
        [Dmat_t,~,Gve,K_] = identify_tangent_v2 ( xx , mat) ;

        I = [ 1 0 0 0 ; 0 1 0 0 ; 0 0 1 0 ; 0 0 0 1] ;    
        M = I + dg * Dmat_t * d2Qdsigma2 ;
        Mi = inv(M) ; 
        Dvp = Mi * Dmat_t -  ( Mi * Dmat_t * dQdsigma * dFdsigma' * Mi * Dmat_t ) / (  eta_vp/dt + h1 + dFdsigma' * Mi * Dmat_t * dQdsigma ) ; 
%                     
        Dmat = Dvp(1:3,1:3); % use viscoplastic tangent
%         maxFvp = max(abs(F_vp),maxFvp);
    end % F_ > 0
%}
        maxF = max(abs(F_),maxF);
    
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
        Residual(sctrB) = Residual(sctrB) - B'*[Sigma_t(1:3)]*W(kk)*det(J0);
        volume_ = volume_ +  W(kk)*det(J0) ; 
%         SS = [SS ; xx Sigma_t' J2 log10(abs(Eii)) uu];
        SS(iel*nq-(nq)+kk,:) = [xx Sigma_t' J2 log10(abs(Eii)) uu];
        HH_(iel*nq-(nq)+kk,1:4) = TSigma' ;
        HH_(iel*nq-(nq)+kk,5) =  P ;
        HH_(iel*nq-(nq)+kk,6) = C ; 
        HH_(iel*nq-(nq)+kk,7:10) = Eps_t; 

    end                  % end of looping on GPs
        
        [i j s] = find(KK); 
        Cell{iel} = [sctrB(i)' ,sctrB(j)', s];
    
        if rem(iel,stepF)==0
            IJV = cell2mat( Cell );
            A = sparse(IJV(:,1),IJV(:,2),IJV(:,3),size(K,1),size(K,1));
            K = K + A;
            Cell = cell(stepF,1);
            c = 0 ; 
        end
end

IJV = cell2mat( Cell );
K = K + sparse(IJV(:,1),IJV(:,2),IJV(:,3),size(K,1),size(K,1));

Fint = Residual ;
% volume_
% ************************
    rightNodes = find ( abs ( node(:,1) - max(node(:,1)) ) < 0.0001 )  ;
    leftNodes = find ( abs ( node(:,1) - min(node(:,1)) ) < 0.0001 )  ;
    topNodes = find ( abs ( node(:,2) - max(node(:,2)) ) < 0.0001 )  ;
    botNodes = find ( abs ( node(:,2) - min(node(:,2)) ) < 0.0001 )  ;

    K0 = K; 
    for  ui = 1 : length(rightNodes) 
        cur_node1 = rightNodes ( ui )  ;
        xx = node(cur_node1,1);%/max(node(:,1)); 
        yy = node(cur_node1,2);%/max(node(:,2)); 
        [K,Residual] = boundary_1point(K,Residual,2*cur_node1,(yy-0.7)*(plast_it==1)*incr0);
        [K,Residual] = boundary_1point(K,Residual,2*cur_node1-1,-xx*(plast_it==1)*incr0);  
    end

    for  ui = 1 : length(botNodes) 
        cur_node1 = botNodes ( ui )  ;
        xx = node(cur_node1,1);%/max(node(:,1)); 
        yy = node(cur_node1,2);%/max(node(:,2)); 
        [K,Residual] = boundary_1point(K,Residual,2*cur_node1,(yy-0.7)*(plast_it==1)*incr0);
        [K,Residual] = boundary_1point(K,Residual,2*cur_node1-1,-xx*(plast_it==1)*incr0);  
    end

    for  ui = 1 : length(topNodes) 
        cur_node1 = topNodes ( ui )  ;
        xx = node(cur_node1,1);%/max(node(:,1)); 
        yy = node(cur_node1,2);%/max(node(:,2)); 
        [K,Residual] = boundary_1point(K,Residual,2*cur_node1,(yy-0.7)*(plast_it==1)*incr0);
%         [K,Residual] = boundary_1point(K,Residual,2*cur_node1-1,-xx*(plast_it==1)*incr0);  
    end

    for  ui = 1 : length(leftNodes) 
        cur_node1 = leftNodes ( ui )  ;
        xx = node(cur_node1,1);%/max(node(:,1)); 
        yy = node(cur_node1,2);%/max(node(:,2)); 
        [K,Residual] = boundary_1point(K,Residual,2*cur_node1,(yy-0.7)*(plast_it==1)*incr0);
        [K,Residual] = boundary_1point(K,Residual,2*cur_node1-1,-xx*(plast_it==1)*incr0);  
    end



%     apply free slip boundary condition at left and right edges
% % % %         ww = 1e8 ; 
% % % %         Kp = sparse(length(K),length(K));
% % % %         Rp = sparse(length(K),1);
% % % %     for  ui = 1 : length(leftNodes) 
% % % %         cur_node1 = leftNodes ( ui )  ;
% % % %         iel = conn_nodes(cur_node1,2);
% % % %         sctr =  element(iel,:) ;
% % % %         sc_ = 2*sctr-1 ; 
% % % %         local = element_coordinate ( node(cur_node1,:) , sctr , node );
% % % %         [N,dNdxi] = lagrange_basis('Q4',local);  % element shape functions
% % % %         J0 = node(sctr,:)'*dNdxi(:,:);                 % element Jacobian matrix
% % % %         xx = N(1:4)'*node(sctr,:) ;
% % % %         invJ0 = inv(J0);
% % % %         dNdx  = dNdxi*invJ0;                      % derivatives of N w.r.t XY
% % % %         scc = dNdx(:,1)' ; 
% % % %         Kp(sc_,sc_) = Kp(sc_,sc_) + ww*scc'*scc;
% % % %         Rp(sc_) = Rp(sc_) - ww*scc'*scc*dU(2*sctr-1) ;
% % % %     end
% % % %     for  ui = 1 : length(rightNodes) 
% % % %         cur_node1 = rightNodes ( ui )  ;
% % % %         iel = conn_nodes(cur_node1,2);
% % % %         sctr =  element(iel,:) ;
% % % %         sc_ = 2*sctr-1 ; 
% % % %         local = element_coordinate ( node(cur_node1,:) , sctr , node );
% % % %         [N,dNdxi] = lagrange_basis('Q4',local);  % element shape functions
% % % %         J0 = node(sctr,:)'*dNdxi(:,:);                 % element Jacobian matrix
% % % %         xx = N(1:4)'*node(sctr,:) ;
% % % %         invJ0 = inv(J0);
% % % %         dNdx  = dNdxi*invJ0;                      % derivatives of N w.r.t XY
% % % %         scc = dNdx(:,1)' ; 
% % % %         Kp(sc_,sc_) = Kp(sc_,sc_) + ww*scc'*scc;
% % % %         Rp(sc_) = Rp(sc_) - ww*scc'*scc*dU(2*sctr-1) ;
% % % %     end
% % % %     if plast_it==1
% % % %         K = K + Kp ; 
% % % %         Residual = Residual + Rp ; 
% % % %     end

    nR = norm(Residual,1)/length(Residual); % residual norm

    fprintf('inc %3d  itr %2d  Resisdual. norm = %2.2e \n', inc , plast_it ,nR)

        % check convergence 
    ddU    = K\(Residual) ;
    dU     = dU + ddU;
        
    if nR < tol_plast % converged 
        break;
    end
    
%     Fint = K*dU ;

end % on global Newton iterations 
%             dU(find(node(:,1)==L)*2-1)
% norm(dU)
% pause
    if exist('tri')==0    tri = delaunay(SS(:,1:2));    end
    VTKPostProcess(SS(:,1:2),tri,1,'Tri3',['stress/stress' num2str(inc)],SS(:,9:11),[SS(:,8) SS(:,7)])

    LD = [ LD ; inc*incr0 mean(sqrt(SS(:,7))) ];% abs(sum(Fint(2*botNodes-1)))  ] ;

    dlmwrite('LD.out',LD,'delimiter',' ')
    % update history parameters:
    HH = HH_;
    dU0 = dU ;

    if rem(inc,1)==0
        figure(h1_)
        clf
        load /home/mohsen/Downloads/m2di/M2Di_EP_KelvinViscoplastic/LD_TD.mat 'ff'
        plot(ff(:,1),ff(:,2),'bsq')
        hold on
        plot(LD(:,1),LD(:,2),'r-sq')
        
%         if plastic == 1 
            fac = 100000 ; 
            node_deformed = node ; 
            node_deformed(:,1) = node(:,1) + fac*dU(1:2:end);
            node_deformed(:,2) = node(:,2) + fac*dU(2:2:end);

%             figure
%             hold on
%             plot(node(:,1),node(:,2),'bsq')
%             plot(node_deformed(:,1),node_deformed(:,2),'rsq')
%             axis equal 
%             plot_mesh(node_deformed,element(:,:),'Q4','k-');
%             pause%(0.1)
%         end       
    end
    
end % on time stepping
% profile viewer
% profsave


figure(h1_)
clf
load /home/mohsen/Downloads/m2di/M2Di_EP_KelvinViscoplastic/LD_TD.mat 'ff'
plot(ff(:,1),ff(:,2),'bsq')
hold on
plot(LD(:,1),LD(:,2),'r-sq')

toc
