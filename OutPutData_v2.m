    fac = 2000;
    node_deformed = node ; 
    node_deformed(:,1) = node(:,1) + fac*dU(1:2:end);
    node_deformed(:,2) = node(:,2) + fac*dU(2:2:end);

    ixg = find(abs(SS0(:,2)-0.5*D)<1) ;

    
if rem(inc,pstep2)==0
    h1 = figure('visible','off') ;
    hold on
    plot(node_deformed(:,1),node_deformed(:,2),'ksq','MarkerSize',3)
    axis equal 
    eval(['print -djpeg out/deformed/jj' num2str(inc) , '.jpeg'])
    close(h1)
end

if rem(inc,pstep2)==0 % LD curve 
	h1 = figure('visible','off');
    subplot(1,2,1)
    hold on
    plot(LD(:,3)/60/60/24/365,LD(:,2)/1e6,'bsq-','MarkerSize',1)
    load Herr1.mat
    plot(Herr1(:,1),Herr1(:,2),'r-')
    ylabel('Stress[MPa]')
    subplot(1,2,2)
    hold on
    plot(LD(:,1),LD(:,2)/1e6,'bsq-','MarkerSize',1)
    ylabel('Stress[MPa]')
    eval(['print -djpeg out/LD' , '.jpeg'])
    close(h1)
    
% maxVp curve 
    h1 = figure('visible','off');
    subplot(1,2,1)
    semilogy(LD(:,3)/60/60/24/365,LD(:,end),'bsq-','MarkerSize',1)
    subplot(1,2,2)
    hold on
    plot(LD(:,1),LD(:,end),'bsq-','MarkerSize',1)
    eval(['print -djpeg out/maxVp' , '.jpeg'])
    close(h1)

    if exist('tri')==0
        tri = delaunay(full(SS(:,1:2)));
        tri2 = delaunay(full(node(:,1:2)));
    end


if (strcmp(Ex,'Herr') || strcmp(Ex,'Preuss') || strcmp(Ex,'Preuss2') ||  strcmp(Ex,'Junction'))&& plastic==1        
    figure('visible','off')
    plot(SS0(:,1),Vp(:,1),'rsq')
    hold on
    plot(SS0(ixg,1),Vp(ixg,1),'bsq')
    ylabel('slip_rate [m/s]')
    ylim([0 inf])
    eval(['print -djpeg out/slip_rate/jj' num2str(inc) , '.jpeg'])
%     dlmwrite(['out/slip_rate/data/jj' num2str(inc) '.log' ],[SS0(inx,1) Vp(inx,1)],'Delimiter',' ')
    save_mat ( 'slip_rate' , [SS0(inx,1) Vp(inx,1)] , inc )     ; 

    figure('visible','off')
    plot(SS0(:,1),Sn(:,1),'rsq')
    hold on
    plot(SS0(ixg,1),Sn(ixg,1),'bsq')
    ylabel('slip [m]')
    ylim([0 inf])
    eval(['print -djpeg out/slip/jj' num2str(inc) , '.jpeg'])
%     dlmwrite(['out/slip/data/jj' num2str(inc) '.log' ],[SS0(inx,1) Sn(inx,1)],'Delimiter',' ')
    save_mat ( 'slip' , [SS0(inx,1) Sn(inx,1)] , inc )     ; 

    

    figure('visible','off')
    plot(SS0(:,1),dg_(:,1),'rsq')
    hold on
    ylabel('dg ')
    eval(['print -djpeg out/dlambda/jj' num2str(inc) , '.jpeg'])

    figure('visible','off')
    hold on
    plot(SS0(ixg,1),Sn(ixg,1),'bsq')
    ylabel('slip [m]')
    ylim([0 inf])
    eval(['print -djpeg out/slip_center/jj' num2str(inc) , '.jpeg'])
    save_mat ( 'slip_center' , [SS0(ixg,1) Sn(ixg,1)] , inc )     ; 


    figure('visible','off')
    plot(SS0(:,1),sin_phi_(:,1),'rsq')
    hold on
    plot(SS0(ixg,1),sin_phi_(ixg,1),'bsq')
    ylabel('mu')
    ylim([0 1])
    eval(['print -djpeg out/mu/jj' num2str(inc) , '.jpeg'])
    save_mat ( 'mu' , [SS0(inx,1) sin_phi_(inx,1)] , inc )     ; 
    
    figure('visible','off')
    semilogy(SS0(:,1),theta_(:,1),'rsq')
    hold on
    ylabel('theta')
    eval(['print -djpeg out/theta/jj' num2str(inc) , '.jpeg'])
    save_mat ( 'theta' , [SS0(inx,1) theta_(inx,1)] , inc )     ; 


end

end

if rem(inc,pstep1)==0
    iii = find(abs(SS0(:,1)-max(SS0(:,1)))<1e-5) ;
    SS0(iii,1) = max(SS0(:,1)) ;
    iii = find(abs(SS0(:,2)-max(SS0(:,2)))<1e-5) ;
    SS0(iii,2) = max(SS0(:,2)) ;
    iii = find(abs(SS0(:,2)-min(SS0(:,2)))<1e-5) ;
    SS0(iii,2) = min(SS0(:,2)) ;


    vtfile =['out/paraview/J/J' num2str(inc) '.vtu'] ;
    id__=find(SS0(:,1)~=0 & SS0(:,2)~=0);
    JJJ = sqrt(J2__);
%     JJJ(isnan(JJJ))= 0 ; 
    VTU_Fiber ([SS0(id__,1:2) SS0(id__,1)*0],tri,1*JJJ(id__),vtfile);


    vtfile =['out/paraview/P/P' num2str(inc) '.vtu'] ;
    JJJ = P__;
%     JJJ(isnan(JJJ))= 0 ; 
% %     VTU_Fiber ([SS0(id__,1:2) SS0(id__,1)*0],tri,1*JJJ(id__),vtfile);
    
    
    vtfile =['out/paraview/Eii/Eii' num2str(inc) '.vtu'] ;
    JJJ = Eii_;
%     JJJ(isnan(JJJ))= 0 ; 
% %     VTU_Fiber ([SS0(id__,1:2) SS0(id__,1)*0],tri,1*JJJ(id__),vtfile);
    
    vtfile =['out/paraview/dggg/dg_' num2str(inc) '.vtu'] ;
    JJJ = dg_;
%     JJJ(isnan(JJJ))= 0 ; 
% %     VTU_Fiber ([SS0(id__,1:2) SS0(id__,1)*0],tri,1*JJJ(id__),vtfile);
    
    vtfile =['out/paraview/vx/vx' num2str(inc) '.vtu'] ;
    JJJ = vu(1:2:end);
% %     VTU_Fiber ([node(:,1:2) node(:,1)*0],tri2,1*JJJ(:),vtfile);
    
    vtfile =['out/paraview/ax/ax' num2str(inc) '.vtu'] ;
    JJJ = au(1:2:end);
    VTU_Fiber ([node(:,1:2) node(:,1)*0],tri2,1*JJJ(:),vtfile);

% %     vtfile =['out/paraview/vy' num2str(inc) '.vtu'] ;
% %     JJJ = vu(2:2:end);
% %     VTU_Fiber ([node(:,1:2) node(:,1)*0],tri2,1*JJJ(:),vtfile);

    vtfile =['out/paraview/sxx/sxx' num2str(inc) '.vtu'] ;
    JJJ = Sigma_t__(:,1);
%     JJJ(isnan(JJJ))= 0 ; 
% %     VTU_Fiber ([SS0(id__,1:2) SS0(id__,1)*0],tri,1*JJJ(id__),vtfile);

    vtfile =['out/paraview/mu/mu' num2str(inc) '.vtu'] ;
    id__=find(SS0(:,1)~=0 & SS0(:,2)~=0);
%     sin_phi_(isnan(Vp))= 0 ; 
% %     VTU_Fiber ([SS0(id__,1:2) SS0(id__,1)*0],tri,1*sin_phi_(id__),vtfile);

    vtfile =['out/paraview/theta/theta' num2str(inc) '.vtu'] ;
    id__=find(SS0(:,1)~=0 & SS0(:,2)~=0);
%     sin_phi_(isnan(Vp))= 0 ; 
% %     VTU_Fiber ([SS0(id__,1:2) SS0(id__,1)*0],tri,1*theta_n(id__),vtfile);

end

% output{inc}{1} = dU;
% output{inc}{2} = Vp;
% output{inc}{3} = Sn;
% output{inc}{4} = dg_;
% output{inc}{5} = sin_phi_;
% output{inc}{6} = theta_;
% output{inc}{7} = J2__;
% output{inc}{8} = vu;
% output{inc}{9} = au;
if inc == 1
    save out/history/data/output0.mat node SS0 SS  D inx
end
% if rem(inc,pstep2)==0
% save out/history/data/output.mat  output LD 
% end
%     eval(['save out/history/data/output' num2str(inc) '.mat LD dU Vp  Sn dg_ sin_phi_ theta_ J2__ vu au'])
%     save(['out/history/data/output' num2str(inc) '.mat'], 'LD','dU','Vp','Sn','dg_','sin_phi_','theta_','J2__','vu','au','-v7.3')
%     save(['out/history/data/output' num2str(inc) '.mat'], 'LD','dU','vu','au','-v7.3')
% ByteSize(output)

    fclose('all'); 
    close all ;
