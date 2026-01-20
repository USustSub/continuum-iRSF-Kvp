    fac = 2000;
    node_deformed = node ; 
    node_deformed(:,1) = node(:,1) + fac*dU(1:2:end);
    node_deformed(:,2) = node(:,2) + fac*dU(2:2:end);

    ixg = find(abs(SS0(:,2)-0.5*D)<1) ;

if rem(inc,pstep2)==0
    
    
% %     if strcmp(Ex,'Herr')
% %         SS_t_deformed = SS_t + fac*u_t ;
% %         SS_b_deformed = SS_b + fac*u_b ;
% %     else
% %         SS_t_deformed = [] ;
% %         SS_b_deformed = [] ; 
% %     end
    
    
    
    h1 = figure('visible','off') ;
%     clf
    hold on
%     plot(node(:,1),node(:,2),'bsq')
    plot(node_deformed(:,1),node_deformed(:,2),'k.','MarkerSize',3)
    axis equal 
%     plot(SS_t_deformed(:,1),SS_t_deformed(:,2),'r.')
%     plot(SS_b_deformed(:,1),SS_b_deformed(:,2),'b.')
    eval(['print -djpeg out/deformed/jj' num2str(inc) , '.jpeg'])
    close(h1)
    
end

% %{
if rem(inc,pstep2)==0
	h1 = figure('visible','off');
    hold on
    if strcmp(Ex,'Duretz')
        load /home/mohsen/Downloads/m2di/M2Di_EP_KelvinViscoplastic/LD.mat  'LD_REF'
        plot(LD_REF(:,1),LD_REF(:,2),'rsq')
        plot(LD(:,1),LD(:,2),'bsq-')
    elseif strcmp(Ex,'Slip')
        load LD_r.mat 'LD_r' 'LD_elas'
        plot(LD_r(:,1),LD_r(:,2),'rsq')
        plot(LD_elas(:,1),LD_elas(:,2),'ksq')
        plot(LD(:,1),LD(:,2),'bsq-')
    elseif strcmp(Ex,'Herr') ||  strcmp(Ex,'Junction')
%         load LD_r.mat 'LD_r' 'LD_elas'
%         plot(LD_r(:,1),LD_r(:,2),'rsq')
%         plot(LD_elas(:,1),LD_elas(:,2),'ksq')
        plot(LD(:,3)/60/60/24/365,LD(:,2)/1e6,'bsq-','MarkerSize',1)
        load Herr1.mat
        plot(Herr1(:,1),Herr1(:,2),'r-')
        ylabel('Stress[MPa]')
    elseif strcmp(Ex,'Preuss') ||  strcmp(Ex,'Preuss2')
        plot(LD(:,3)/60/60/24/365,LD(:,2)/1e6,'bsq-','MarkerSize',1)
        ylabel('Stress[MPa]')
    end
    eval(['print -djpeg out/LD' , '.jpeg'])
    close(h1)
    
    
    h1 = figure('visible','off');
    semilogy(LD(:,3)/60/60/24/365,LD(:,end),'bsq-','MarkerSize',1)
    eval(['print -djpeg out/maxVp' , '.jpeg'])
    close(h1)
    
% %{
    if exist('tri')==0
        tri = delaunay(full(SS(:,1:2)));
        tri2 = delaunay(full(node(:,1:2)));
    end
    
    
% %     SS_DU = SS*0 ;
% % for ii__ = 1  :size(SS,1)
% %    [phi,B,en] = get_data ( xx , node , di , form ) ;
% %     phi = phi_SS{ii__} ; 
% %     en = index_SS{ii__} ;
% %     SS_DU (ii__,: )  = [ phi*dU(en(1:2:end)) phi*dU(en(2:2:end)) ]; 
% % end
    
    
%     VTKPostProcess(SS(:,1:2),tri,1,'Tri3',['stress/stress' num2str(inc)],SS(:,3:5),[SS(:,6) SS(:,7)])
% %     vtfile =['out/paraview/Eii' num2str(inc) '.vtu'] ;
% %     EE = (Eii_); EE(EE==0) = []; EE = log10(abs(EE)) ; 
% % %     VTU_Fiber ([SS(:,1:2) SS(:,1)*0],tri,EE,vtfile);
% %     VTU_Fiber ([SS(:,1:2) SS(:,1)*0],tri,EE,vtfile);

%     vtfile =['out/paraview/SQJ' num2str(inc) '.vtu'] ;
%     VTU_Fiber ([SS(:,1:2) SS(:,1)*0],tri,SS(:,7),vtfile); 
  %}    


if (strcmp(Ex,'Herr') || strcmp(Ex,'Preuss') || strcmp(Ex,'Preuss2') ||  strcmp(Ex,'Junction'))&& plastic==1

        
    figure('visible','off')
    plot(SS0(:,1),Vp(:,1),'rsq')
%     ylim([0 1])
    hold on
    plot(SS0(ixg,1),Vp(ixg,1),'bsq')
    ylabel('slip_rate [m/s]')
    ylim([0 inf])
    eval(['print -djpeg out/slip_rate/jj' num2str(inc) , '.jpeg'])

    
    dlmwrite(['out/slip_rate/jj' num2str(inc) '.log' ],[SS0(inx,1) Vp(inx,1)],'Delimiter',' ')
    
    
    figure('visible','off')
    plot(SS0(:,1),Sn(:,1),'rsq')
    hold on
    plot(SS0(ixg,1),Sn(ixg,1),'bsq')
%     ylim([0 1])
    ylabel('slip [m]')
    ylim([0 inf])
    eval(['print -djpeg out/slip/jj' num2str(inc) , '.jpeg'])
    dlmwrite(['out/slip/jj' num2str(inc) '.log' ],[SS0(inx,1) Sn(inx,1)],'Delimiter',' ')

    
    
    figure('visible','off')
    plot(SS0(:,1),dg_(:,1),'rsq')
    hold on
%     ylim([0 1])
    ylabel('dg ')
%     ylim([0 inf])
    eval(['print -djpeg out/dlambda/jj' num2str(inc) , '.jpeg'])
%     dlmwrite(['out/dlambda/jj' num2str(inc) '.log' ],[SS0(inx,1) dg_(inx,1)],'Delimiter',' ')

    
    
    
    figure('visible','off')
%     plot(SS0(:,1),Sn(:,1),'rsq')
    hold on
    plot(SS0(ixg,1),Sn(ixg,1),'bsq')
%     ylim([0 1])
    ylabel('slip [m]')
    ylim([0 inf])
    eval(['print -djpeg out/slip_center/jj' num2str(inc) , '.jpeg'])
%     dlmwrite(['out/slip_center/jj' num2str(inc) '.log' ],[SS0(ixg,1) Sn(ixg,1)],'Delimiter',' ')

    
    figure('visible','off')
    plot(SS0(:,1),sin_phi_(:,1),'rsq')
    hold on
    plot(SS0(ixg,1),sin_phi_(ixg,1),'bsq')
%     ylim([0 1])
    ylabel('mu')
    ylim([0 1])
    eval(['print -djpeg out/mu/jj' num2str(inc) , '.jpeg'])
     dlmwrite(['out/mu/jj' num2str(inc) '.log' ],[SS0(inx,1) sin_phi_(inx,1)],'Delimiter',' ')
    
    figure('visible','off')
    semilogy(SS0(:,1),theta_(:,1),'rsq')
    hold on
%     if inc  > 11
%         plot(SS0(ixg,1),theta_(ixg,1),'bsq')
%     end
%     ylim([0 1])
    ylabel('theta')
    eval(['print -djpeg out/theta/jj' num2str(inc) , '.jpeg'])
    dlmwrite(['out/theta/jj' num2str(inc) '.log' ],[SS0(inx,1) theta_(inx,1)],'Delimiter',' ')

%     figure('visible','off')
%     plot(DD(:,1),slip(:,1),'rsq')
%     ylabel('slip')
%     eval(['print -djpeg out/slip/jj' num2str(inc) , '.jpeg'])

%     figure('visible','off')
%     plot(Slip(:,1),Slip(:,?),'rsq')
%     ylabel('slip')
%     eval(['print -djpeg out/slip/jj' num2str(inc) , '.jpeg'])

% %     vtfile =['out/paraview/Eii_p' num2str(inc) '.vtu'] ;
% %     id__=find(SS0(:,1)~=0 & SS0(:,2)~=0);
% %     Eii_p(isnan(Eii_p))= 0 ; 
% %     VTU_Fiber ([SS0(id__,1:2) SS0(id__,1)*0],tri,1e15*Eii_p(id__),vtfile);

% %     vtfile =['out/paraview/Vp' num2str(inc) '.vtu'] ;
% %     id__=find(SS0(:,1)~=0 & SS0(:,2)~=0);
% %     Vp(isnan(Vp))= 0 ; 
% %     VTU_Fiber ([SS0(id__,1:2) SS0(id__,1)*0],tri,1*Vp(id__),vtfile);


end

end

if rem(inc,pstep1)==0
    iii = find(abs(SS0(:,1)-max(SS0(:,1)))<1e-5) ;
    SS0(iii,1) = max(SS0(:,1)) ;
    iii = find(abs(SS0(:,2)-max(SS0(:,2)))<1e-5) ;
    SS0(iii,2) = max(SS0(:,2)) ;
    iii = find(abs(SS0(:,2)-min(SS0(:,2)))<1e-5) ;
    SS0(iii,2) = min(SS0(:,2)) ;

    vtfile =['out/paraview/dggg/dggg' num2str(inc) '.vtu'] ;
    id__=find(SS0(:,1)~=0 & SS0(:,2)~=0);
    dggg = dg_2; 
    dggg(isnan(dggg))= 0 ; 
    VTU_Fiber ([SS0(id__,1:2) SS0(id__,1)*0],tri,1*dggg(id__),vtfile);

    vtfile =['out/paraview/J/J' num2str(inc) '.vtu'] ;
    id__=find(SS0(:,1)~=0 & SS0(:,2)~=0);
    JJJ = sqrt(J2__);
%     JJJ(isnan(JJJ))= 0 ; 
    VTU_Fiber ([SS0(id__,1:2) SS0(id__,1)*0],tri,1*JJJ(id__),vtfile);


    vtfile =['out/paraview/P/P' num2str(inc) '.vtu'] ;
    JJJ = P__;
%     JJJ(isnan(JJJ))= 0 ; 
    VTU_Fiber ([SS0(id__,1:2) SS0(id__,1)*0],tri,1*JJJ(id__),vtfile);
    
    
    vtfile =['out/paraview/Eii/Eii' num2str(inc) '.vtu'] ;
    JJJ = Eii_;
%     JJJ(isnan(JJJ))= 0 ; 
    VTU_Fiber ([SS0(id__,1:2) SS0(id__,1)*0],tri,1*JJJ(id__),vtfile);
    
    
    vtfile =['out/paraview/vx/vx' num2str(inc) '.vtu'] ;
    JJJ = vu(1:2:end);
    VTU_Fiber ([node(:,1:2) node(:,1)*0],tri2,1*JJJ(:),vtfile);
    
    vtfile =['out/paraview/ax/ax' num2str(inc) '.vtu'] ;
    JJJ = au(1:2:end);
    VTU_Fiber ([node(:,1:2) node(:,1)*0],tri2,1*JJJ(:),vtfile);

% %     vtfile =['out/paraview/vy' num2str(inc) '.vtu'] ;
% %     JJJ = vu(2:2:end);
% %     VTU_Fiber ([node(:,1:2) node(:,1)*0],tri2,1*JJJ(:),vtfile);

    vtfile =['out/paraview/sxx/sxx' num2str(inc) '.vtu'] ;
    JJJ = Sigma_t__(:,1);
%     JJJ(isnan(JJJ))= 0 ; 
    VTU_Fiber ([SS0(id__,1:2) SS0(id__,1)*0],tri,1*JJJ(id__),vtfile);
    
    
    vtfile =['out/paraview/syy/syy' num2str(inc) '.vtu'] ;
    JJJ = Sigma_t__(:,2);
%     JJJ(isnan(JJJ))= 0 ; 
    VTU_Fiber ([SS0(id__,1:2) SS0(id__,1)*0],tri,1*JJJ(id__),vtfile);
    
    
    vtfile =['out/paraview/sxy/sxy' num2str(inc) '.vtu'] ;
    JJJ = Sigma_t__(:,3);
%     JJJ(isnan(JJJ))= 0 ; 
    VTU_Fiber ([SS0(id__,1:2) SS0(id__,1)*0],tri,1*JJJ(id__),vtfile);
    
    
    vtfile =['out/paraview/szz/szz' num2str(inc) '.vtu'] ;
    JJJ = Sigma_t__(:,4);
%     JJJ(isnan(JJJ))= 0 ; 
    VTU_Fiber ([SS0(id__,1:2) SS0(id__,1)*0],tri,1*JJJ(id__),vtfile);

    vtfile =['out/paraview/mu/mu' num2str(inc) '.vtu'] ;
    id__=find(SS0(:,1)~=0 & SS0(:,2)~=0);
%     sin_phi_(isnan(Vp))= 0 ; 
    VTU_Fiber ([SS0(id__,1:2) SS0(id__,1)*0],tri,1*sin_phi_(id__),vtfile);

    vtfile =['out/paraview/theta/theta' num2str(inc) '.vtu'] ;
    id__=find(SS0(:,1)~=0 & SS0(:,2)~=0);
%     sin_phi_(isnan(Vp))= 0 ; 
    VTU_Fiber ([SS0(id__,1:2) SS0(id__,1)*0],tri,1*theta_n(id__),vtfile);

    
    
end


if (inc == 1 && (  strcmp(Ex,'Slip') ||  strcmp(Ex,'Junction') ||  strcmp(Ex,'Herr') ||  strcmp(Ex,'Preuss') ||  strcmp(Ex,'Preuss2')) ) 
        figure('visible','off')
        plot(node(:,1),node(:,2),'b.','MarkerSize',1)
        hold on
        plot(DD(:,1),DD(:,2),'rsq','MarkerSize',1) 
        axis equal
        plot(SS0(ixg,1),SS0(ixg,2),'b.')
        eval(['print -djpeg out/DD' , '.jpeg'])
%     pause
end

%}







    fclose('all'); 
    close all ;
    
%%
% % % pause
% % % ff = square_node_array([0 0],[L 0],[L D],[0 D],200,200);
% % ff = [0.:10:D];
% % ff = [ ff'*0+L/2 ff'];
% % U_ = zeros(length(ff),2);
% % for ii = 1 : size(ff,1)
% %     xx = ff(ii,:) ; 
% %     [index] = define_support(node,xx,di);
% %     [phi,~,~] = MLS_ShapeFunction(xx,index,node,di,form);
% %     u__ = [phi*dU(index*2-1) phi*dU(index*2)] ;  
% %     U_(ii,:) = u__ ; 
% % end
% % % n_d = ff + fac*[U_(:,1) U_(:,2)];
% % % figure
% % % plot(n_d(:,1),n_d(:,2),'bsq')
% % % hold on
% % % plot(node_deformed(:,1),node_deformed(:,2),'ksq')
% % figure
% % plot(ff(:,2),U_(:,1),'r-sq')