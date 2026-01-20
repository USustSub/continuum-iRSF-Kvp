% initialize stresses at all nodes 
TSigma0 = 0 ; P0 = 0 ; 

HH_ = 0*HH ;
HH_TT = 0*HH ;
HH_2 = zeros(size(node,1),4*4); % history parameter
HH_2_ = 0*HH_2 ;
HH_3 = zeros(size(node,1),4*1); % history parameteri
HH_3_ = 0*HH_3+sin_phi ; 
sin_phi_ = [HH_3_(:,1);HH_3_(:,2);HH_3_(:,3);HH_3_(:,4)];

HH_4 = zeros(size(node,1),4*1); % history parameter
HH_4_ = 0*HH_4+theta_ ; 

% initial C
C0 = [HH(:,1*6-0);HH(:,2*6-0);HH(:,3*6-0);HH(:,4*6-0)];
           
% initial P
P0_ = p0 + [HH(:,1*6-1);HH(:,2*6-1);HH(:,3*6-1);HH(:,4*6-1)];
HH(:,1*6-1) = p0 ; 
HH(:,2*6-1) = p0 ; 
HH(:,3*6-1) = p0 ; 
HH(:,4*6-1) = p0 ; 


MM = sparse(2*numnode,2*numnode);

sin_phi_ = 0*sin_phi_+0.000000;  
mu_gp = sin_phi_ ; 

% a and b params
if strcmp(Ex,'Herr') || strcmp(Ex,'Preuss') || strcmp(Ex,'Preuss2')
    b_ = zeros(length(SS0),1)+b2_;
    b_((SS0(:,1)<LL1 | SS0(:,1)>LL2))=b1_;
%     id1_= find(SS0(:,1)<LL1 & SS0(:,1)>LL1-4000)   ; 
%     b_(id1_)=b1_+(SS0(id1_,1)-28000)* (b2_-b1_)/4000;
%     id2_= find(SS0(:,1)>LL2 & SS0(:,1)<LL2+4000)   ;
%     b_(id2_)=b2_+(SS0(id2_,1)-LL2)* (b1_-b2_)/4000;
%     b_((C0~=0)) = 0.005 ; 
%     b_((SS0(:,1)<2000 | SS0(:,1)>L-2000))=0.000001;
    if strcmp(Ex,'Herr') 
        b_((C0~=0)) = b2_ ; 
    end
end
if strcmp(Ex,'Junction') 
    b_ = zeros(length(SS0),1)+b2_;
%     b_((SS0(:,1)<LL1 | SS0(:,1)>LL2))=b1_;
%     b_((SS0(:,2)<LL1 | SS0(:,2)>LL2))=b1_;
    b_((SS0(:,1)<LL1 | SS0(:,1)>LL2))=b1_;
    b_((C0~=0)) = b2_ ; 
    
    if any(m_type(:,1))~=0 % for multiple faults
        b_ = 0*b_+b2_ ; 
        b_(m_type(:,1)==1) = b2_ ; % velocity strengthening
        b_(m_type(:,1)==2) = b1_ ; % velocity weakening
    end
end

theta_n = theta_ ; 
Vp = 0 ; LD = [] ;

% get initial elastic solution matrix
stype = 'elastic'; 

dU_e = zeros(2*numnode,1); 
dU_e = 0*dU_e ;
dU = dU_e;
dU0 = 0*dU_e ; 
plastic = 0 ; % set initial flag
u = sparse(2*numnode,1);
au = u ;
vu = u ;
u0 = u ;
au0 = u ;
vu0 = u ;
time = 0 ; 
C_BC = [] ; 
K_BC_ = sparse(2*numnode,2*numnode)  ; 
topNodes = find(node(:,2)==D); 
botNodes = find(node(:,2)==0); 
AP_0 = 0 ; S0 = 0 ; 
ittt = 0 ; 


% figure
% plot(DD(:,1),DD(:,2),'ksq')
% hold on
% plot(node(:,1),node(:,2),'rsq')
% pause


if strcmp(Ex,'Preuss') || strcmp(Ex,'Preuss2')
    theta_n((C0~=0)) = L_/V0*exp(5) ;
    theta_((C0~=0)) = L_/V0*exp(5) ;
    theta_n((C0==0)) = L_/V0*exp(-1) ;
    theta_((C0==0)) = L_/V0*exp(-1) ;
    theta_n = theta_n'; 
    theta_ = theta_'; 
elseif strcmp(Ex,'Herr')
    theta_n((C0~=0)) = L_/V0*exp(40) ;
    theta_((C0~=0)) = L_/V0*exp(40) ;
    theta_n((C0==0)) = L_/V0*exp(-1) ;
    theta_((C0==0)) = L_/V0*exp(-1) ;
    theta_n = theta_n'; 
    theta_ = theta_'; 
elseif strcmp(Ex,'Junction')
    theta_n((C0~=0)) = L_/V0*exp(40) ;
    theta_((C0~=0)) = L_/V0*exp(40) ;
    theta_n((C0==0)) = L_/V0*exp(-1) ;
    theta_((C0==0)) = L_/V0*exp(-1) ;
    theta_n = theta_n'; 
    theta_ = theta_'; 
    
    
%     theta_n = T_type; 
%     theta_ = T_type; 
    
    
end


%     if any(m_type(:,2))~=0 % for multiple faults
%         m_type(m_type(:,2)<1e-2,2) = 0;
%         theta_n (:,1) = m_type(:,2).*L_/V0*exp(40)+L_/V0*exp(-1);
%         theta_ (:,1) = m_type(:,2).*L_/V0*exp(40)+L_/V0*exp(-1);
%     end

if exist('tri')==0
    tri = delaunay(full(SS(:,1:2)));
    tri2 = delaunay(full(node(:,1:2)));
end
    vtfile =['out/paraview/theta_0' num2str(inc) '.vtu'] ;
    id__=find(SS0(:,1)~=0 & SS0(:,2)~=0);
    VTU_Fiber ([SS0(id__,1:2) SS0(id__,1)*0],tri,1*theta_n(id__),vtfile);
        
    vtfile =['out/paraview/m_type' num2str(inc) '.vtu'] ;
    id__=find(SS0(:,1)~=0 & SS0(:,2)~=0);
    VTU_Fiber ([SS0(id__,1:2) SS0(id__,1)*0],tri,1*m_type(id__,2),vtfile);

    vtfile =['out/paraview/b_' num2str(inc) '.vtu'] ;
    id__=find(SS0(:,1)~=0 & SS0(:,2)~=0);
    VTU_Fiber ([SS0(id__,1:2) SS0(id__,1)*0],tri,1*b_(id__),vtfile);


    

C0 = 0*C0 ; 
C0((SS0(:,1)==0 & SS0(:,2)==0))=1e16;



if ifprofile == 1
    ninc = 100 ; 
end


profile off 
if ifprofile == 1
profile on
end


