   
clear all ; 
close all ; 
    load out/history/data/output0.mat node SS0 SS  D inx
%     load out/history/data/output.mat
    L = 150000 ; 
    D = 150000 ; 
    
    ii__ = size(find(node(:,2)==0),1) ; 
    nnx = ii__  ;
    nny = nnx ; 
    inc_u = 1;
    inc_v = nnx;
    node_pattern = [ 1 2 nnx+2 nnx+1 ];
    element = make_elem(node_pattern,nnx-1,nny-1,inc_u,inc_v);

    top_id = zeros(length(element),1); 
    bot_id = zeros(length(element),1); 
    for ii = 1 : size(element,1)
        sctr = element(ii,:) ;
        n_sctr = node(sctr,:) ; 
        if all(n_sctr(:,2)>D/2)
            top_id(ii) = 1 ; 
        elseif all(n_sctr(:,2)<D/2)
            bot_id(ii) = 1 ;         
        end
    end

    top_node_id = zeros(length(node),1); 
    bot_node_id = zeros(length(node),1); 
    for ii = 1 : size(node,1)
        n_sctr = node(ii,:) ; 
        if all(n_sctr(:,2)>D/2)
            top_node_id(ii) = 1 ; 
        elseif all(n_sctr(:,2)<D/2)
            bot_node_id(ii) = 1 ;         
        end
    end
    
        node_top = node(top_node_id==1,:) ;
        node_bot = node(bot_node_id==1,:) ;
        BC_top = find(abs(node_top(:,2)-min(node_top(:,2)))<0.0001);
        BC_bot = find(abs(node_bot(:,2)-max(node_bot(:,2)))<0.0001);

        BC_top2 = find(abs(node_top(:,2)-max(node_top(:,2)))<0.0001);
        BC_top3 = find(abs(node_top(:,1)-max(node_top(:,1)))<0.0001);
        BC_top4 = find(abs(node_top(:,1)-min(node_top(:,1)))<0.0001);

        BC_bot2 = find(abs(node_bot(:,2)-min(node_bot(:,2)))<0.0001);
        BC_bot3 = find(abs(node_bot(:,1)-max(node_bot(:,1)))<0.0001);
        BC_bot4 = find(abs(node_bot(:,1)-min(node_bot(:,1)))<0.0001);



        
%%
        h1 = figure('visible','on') ;

for inc =  [ 1:20000 ]
    inc
        eval(['load out/history/data/output' num2str(inc) '.mat '])
%     vu = output{inc}{8} ;
%     au = output{inc}{9} ;
%     dU = output{inc}{1} ;

    if exist('tri')==0
        tri = delaunay(full(SS(:,1:2)));
        tri2 = delaunay(full(node(:,1:2)));
    end

        
        fac = 5000;
        node_deformed = node ; 
        node_deformed(:,1) = node(:,1) + fac*dU(1:2:end);
        node_deformed(:,2) = node(:,2) + fac*dU(2:2:end);
        au_ = [ au(1:2:end,1) au(2:2:end,1)] ;
        au_top = au_(top_node_id==1,:) ;
        au_bot = au_(bot_node_id==1,:) ;

        
        node_top = node_deformed(top_node_id==1,:) ;
        node_bot = node_deformed(bot_node_id==1,:) ;
        node_top(:,2) = node_top(:,2) - 0 ; 
        node_bot(:,2) = node_bot(:,2) + 0 ; 
        
        
        ixg = find(abs(SS0(:,2)-0.5*D)<1) ;
        
        
%         figure(h1)
%         plot(node_deformed(top_node_id==1,1),node_deformed(top_node_id==1,2),'r.','MarkerSize',2)
        hold on
%         plot(node_deformed(bot_node_id==1,1),node_deformed(bot_node_id==1,2),'b.','MarkerSize',2)
        axis equal 
        xlim([-L/5 L+L/5 ])
        ylim([-D/10 D+D/10 ])
        axis off
        mTextBox = uicontrol('style','text') ;
        set(mTextBox,'String',['i=' num2str(inc) ]) ;
        plot(node_top(BC_top,1),node_top(BC_top,2),'k-','LineWidth',1.5 )
        plot(node_top(BC_top2,1),node_top(BC_top2,2),'k-','LineWidth',1.5)
        plot(node_top(BC_top3,1),node_top(BC_top3,2),'k-','LineWidth',1.5)
        plot(node_top(BC_top4,1),node_top(BC_top4,2),'k-','LineWidth',1.5)
        plot(node_bot(BC_bot,1),node_bot(BC_bot,2),'k-','LineWidth',1.5)
        plot(node_bot(BC_bot2,1),node_bot(BC_bot2,2),'k-','LineWidth',1.5)
        plot(node_bot(BC_bot3,1),node_bot(BC_bot3,2),'k-','LineWidth',1.5)
        plot(node_bot(BC_bot4,1),node_bot(BC_bot4,2),'k-','LineWidth',1.5)
        
        
        X2 = reshape(node_top(:,1),nnx,nny/2);
        Y2 = reshape(node_top(:,2),nnx,nny/2);
        Z2 = reshape(au_top(:,1),nnx,nny/2);
        contourf(X2,Y2,Z2,'edgecolor','none')
        colorbar
        X2 = reshape(node_bot(:,1),nnx,nny/2);
        Y2 = reshape(node_bot(:,2),nnx,nny/2);
        Z2 = reshape(au_bot(:,1),nnx,nny/2);
        contourf(X2,Y2,Z2,'edgecolor','none')
        colorbar
        colormap(hot)
       caxis([-0.1 0.1]/1)
%         colormap(jet(50))
%         pause
        
        fpat = 'out' ;
        fpat2 = 'deformed/data/' ;
        filename = [fpat,filesep,fpat2,filesep,'d' num2str(inc) ];
        print(h1,'-djpeg','-r250',filename)
        clf

%     vtfile =['out/paraview/ax/ax_top' num2str(inc) '.vtu'] ;
%     JJJ = au(1:2:end);
%     VTU_Fiber ([node_deformed(:,1:2) node_deformed(:,1)*0],element(top_id==1,:),1*JJJ(:),vtfile);
%     vtfile =['out/paraview/ax/ax_bot' num2str(inc) '.vtu'] ;
%     JJJ = au(1:2:end);
%     VTU_Fiber ([node_deformed(:,1:2) node_deformed(:,1)*0],element(bot_id==1,:),1*JJJ(:),vtfile);
%     vtfile =['out/paraview/ax/ax' num2str(inc) '.vtu'] ;
%     JJJ = au(1:2:end);
%     VTU_Fiber ([node_deformed(:,1:2) node_deformed(:,1)*0],element(:,:),1*JJJ(:),vtfile);

    
%     pause(0.01)
%         clf
end

%%
figure
[X,Y] = meshgrid(1:0.5:10,1:20);
Z = sin(X) + cos(Y);
surface(X,Y,Z)
%%
au_ = [ au(1:2:end,1) au(2:2:end,1)] ;
X = zeros(nnx,nny);
Y = zeros(nnx,nny);
Z = zeros(nnx,nny);
figure
% for ii = 1 : nnx
%     for jj = 1 : nny
%         X(ii,jj) = node( nny*(jj-1)+ii , 1) ;
%         Y(ii,jj) = node( nny*(jj-1)+ii , 2) ;
%         Z(ii,jj) = au_( nny*(jj-1)+ii , 1) ;
%         
%     end
% end
X2 = reshape(node(:,1),nnx,nny);
Y2 = reshape(node(:,2),nnx,nny);
Z2 = reshape(au(1:2:end,1),nnx,nny);

contourf(X2,Y2,Z2)
% contourslice(X,Y,Z*0,Z,[],[],0:5:60)

