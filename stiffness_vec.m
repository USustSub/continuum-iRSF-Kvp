        D11_ = D11 ;         D12_ = D12 ;          D13_ = D13 ;         D14_ = D14 ; 
        D21_ = D21 ;         D22_ = D22 ;          D23_ = D23 ;         D24_ = D24 ; 
        D31_ = D31 ;         D32_ = D32 ;          D33_ = D33 ;         D34_ = D34 ; 
        D41_ = D41 ;         D42_ = D42 ;          D43_ = D43 ;         D44_ = D44 ; 
        
        inx = find(F__>0) ;
        if isempty(inx)==0
            D11_(inx) = Dvp_(inx,1) ;
            D12_(inx) = Dvp_(inx,2) ;
            D13_(inx) = Dvp_(inx,3) ;
            D14_(inx) = Dvp_(inx,4) ;

            D21_(inx) = Dvp_(inx,5) ;
            D22_(inx) = Dvp_(inx,6) ;
            D23_(inx) = Dvp_(inx,7) ;
            D24_(inx) = Dvp_(inx,8) ;

            D31_(inx) = Dvp_(inx,9) ;
            D32_(inx) = Dvp_(inx,10) ;
            D33_(inx) = Dvp_(inx,11) ;
            D34_(inx) = Dvp_(inx,12) ;

            D41_(inx) = Dvp_(inx,13) ;
            D42_(inx) = Dvp_(inx,14) ;
            D43_(inx) = Dvp_(inx,15) ;
            D44_(inx) = Dvp_(inx,16) ;
        end
            FX = zeros(size(D11_,1),size(Bx,2)) ;
            FY = zeros(size(D11_,1),size(Bx,2)) ;
            FX2 = zeros(size(D11_,1),size(Bx,2)) ;
            FY2 = zeros(size(D11_,1),size(Bx,2)) ;
            for hhh_ = 1 : size(Bx,2)
                f1__= (D11_.*Bx(:,hhh_)+D13_.*By(:,hhh_))/deltaX;
                f2__= (D12_.*By(:,hhh_)+D13_.*Bx(:,hhh_))/deltaX;
                FX(:,hhh_*2-1:hhh_*2) = [ f1__ f2__];
                
                f1__= (D31_.*Bx(:,hhh_)+D33_.*By(:,hhh_))/deltaX;
                f2__= (D32_.*By(:,hhh_)+D33_.*Bx(:,hhh_))/deltaX;
                FY(:,hhh_*2-1:hhh_*2) = [ f1__ f2__];
                
                f1__= (D31_.*Bx(:,hhh_)+D33_.*By(:,hhh_))/deltaY;
                f2__= (D32_.*By(:,hhh_)+D33_.*Bx(:,hhh_))/deltaY;
                FX2(:,hhh_*2-1:hhh_*2) = [ f1__ f2__];
                
                f1__= (D21_.*Bx(:,hhh_)+D23_.*By(:,hhh_))/deltaY;
                f2__= (D22_.*By(:,hhh_)+D23_.*Bx(:,hhh_))/deltaY;
                FY2(:,hhh_*2-1:hhh_*2) = [ f1__ f2__];
            end

%{     
%             pause
%             GX = sparse(4*numnode,2*numnode) ;
%             GY = sparse(4*numnode,2*numnode) ;
%             GX2 = sparse(4*numnode,2*numnode) ;
%             GY2 = sparse(4*numnode,2*numnode) ;

%             for hhh_ = 1 : size(Bx,2)
%                 f1__= (D11_.*Bx(:,hhh_)+D13_.*By(:,hhh_))/deltaX;
%                 f2__= (D12_.*By(:,hhh_)+D13_.*Bx(:,hhh_))/deltaX;
%                 en1__ = En(1:numnode,hhh_*2-1); en1__(en1__==0) = 1;
%                 en2__ = En(1:numnode,hhh_*2); en2__(en2__==0) = 1;
%                 GX(:,en1__) = GX(:,en1__) + f1__ ;
%                 GX(:,en2__) = GX(:,en2__) + f2__ ;
%                 
%                 f1__= (D31_.*Bx(:,hhh_)+D33_.*By(:,hhh_))/deltaX;
%                 f2__= (D32_.*By(:,hhh_)+D33_.*Bx(:,hhh_))/deltaX;
%                 en1__ = En(1:numnode,hhh_*2-1); en1__(en1__==0) = 1;
%                 en2__ = En(1:numnode,hhh_*2); en2__(en2__==0) = 1;
%                 GY(:,en1__) = GY(:,en1__) + f1__ ;
%                 GY(:,en2__) = GY(:,en2__) + f2__ ;
%                 
%                 f1__= (D31_.*Bx(:,hhh_)+D33_.*By(:,hhh_))/deltaY;
%                 f2__= (D32_.*By(:,hhh_)+D33_.*Bx(:,hhh_))/deltaY;
%                 en1__ = En(1:numnode,hhh_*2-1); en1__(en1__==0) = 1;
%                 en2__ = En(1:numnode,hhh_*2); en2__(en2__==0) = 1;
%                 GX2(:,en1__) = GX2(:,en1__) + f1__ ;
%                 GX2(:,en2__) = GX2(:,en2__) + f2__ ;
% 
%                 f1__= (D21_.*Bx(:,hhh_)+D23_.*By(:,hhh_))/deltaY;
%                 f2__= (D22_.*By(:,hhh_)+D23_.*Bx(:,hhh_))/deltaY;
%                 en1__ = En(1:numnode,hhh_*2-1); en1__(en1__==0) = 1;
%                 en2__ = En(1:numnode,hhh_*2); en2__(en2__==0) = 1;
%                 GY2(:,en1__) = GY2(:,en1__) + f1__ ;
%                 GY2(:,en2__) = GY2(:,en2__) + f2__ ;
% 
%             end
            


            
%             K(1:2:2*numnode,:) = K(1:2:2*numnode,:) - GX(1:numnode,:);
%             K(2:2:2*numnode,:) = K(1:2:2*numnode,:) - GY(1:numnode,:);
%             K(1:2:2*numnode,:) = K(1:2:2*numnode,:) + GX(numnode+1:2*numnode,:);
%             K(2:2:2*numnode,:) = K(1:2:2*numnode,:) + GY(numnode+1:2*numnode,:);
%             K(1:2:2*numnode,:) = K(1:2:2*numnode,:) - GX2(2*numnode+1:3*numnode,:);
%             K(2:2:2*numnode,:) = K(1:2:2*numnode,:) - GY2(2*numnode+1:3*numnode,:);
%             K(1:2:2*numnode,:) = K(1:2:2*numnode,:) + GX2(3*numnode+1:4*numnode,:);
%             K(2:2:2*numnode,:) = K(1:2:2*numnode,:) + GY2(3*numnode+1:4*numnode,:);

% K = 0*K; 
% K2 = K  ; 

K3 = K  ; 

        stepF = 10 ;
        Cell = cell(stepF ,1);  Cell2 = cell(stepF ,1); c_ = 1 ;
        C1= [] ; 
            for ij = 1 : numnode
                node_c = node(ij,:) ; 
                if (node_c(1) > 0 && node_c(1) < L && node_c(2) < D && node_c(2) > 0 )

                    idx = 1; 
                    id1 = ij+idx*numnode-numnode ; 
                    en = En(id1,:); en = en(en~=0); en = en(en>0);
                    kx = FX(id1,1:length(en));
                    ky = FY(id1,1:length(en));
                    K(2*ij-1,en) = K(2*ij-1,en) - kx ; 
                    K(2*ij-0,en) = K(2*ij-0,en) - ky ;

                    i=[1:1:length(kx)]; j = 1:length(kx) ; s = -kx ;  
                    C1 = [ C1;  [en(i)'*0+2*ij-1 ,en(j)', s'] ] ;
%                     pause
                    Cell{c_} = [en(i)'*0+2*ij-1 ,en(j)', s'];
                    i=[1:1:length(ky)]; j = 1:length(ky) ; s = -ky ;  
                    Cell2{c_} = [en(i)'*0+2*ij ,en(j)', s'];
                    c_ = c_ + 1 ; 
                    C1 = [ C1;  [en(i)'*0+2*ij-0 ,en(j)', s'] ] ;

                    idx = 2; 
                    id1 = ij+idx*numnode-numnode ; 
                    en = En(id1,:); en = en(en~=0); en = en(en>0);
                    kx = FX(id1,1:length(en));
                    ky = FY(id1,1:length(en));
                    K(2*ij-1,en) = K(2*ij-1,en) + kx ; 
                    K(2*ij-0,en) = K(2*ij-0,en) + ky ; 

                    i=[1:1:length(kx)]; j = 1:length(kx) ; s = kx ;  
                    Cell{c_} = [en(i)'*0+2*ij-1 ,en(j)', s'];
                    C1 = [ C1;  [en(i)'*0+2*ij-1 ,en(j)', s'] ] ;
                    i=[1:1:length(ky)]; j = 1:length(ky) ; s = ky ;  
                    Cell2{c_} = [en(i)'*0+2*ij ,en(j)', s'];
                    C1 = [ C1;  [en(i)'*0+2*ij-0 ,en(j)', s'] ] ;
                    c_ = c_ + 1 ; 

                    idx = 3; 
                    id1 = ij+idx*numnode-numnode ; 
                    en = En(id1,:); en = en(en~=0); en = en(en>0);
                    kx = FX2(id1,1:length(en));
                    ky = FY2(id1,1:length(en));
                    K(2*ij-1,en) = K(2*ij-1,en) - kx ; 
                    K(2*ij-0,en) = K(2*ij-0,en) - ky ; 

                    i=[1:1:length(kx)]; j = 1:length(kx) ; s = -kx ;  
                    Cell{c_} = [en(i)'*0+2*ij-1 ,en(j)', s'];
                    C1 = [ C1;  [en(i)'*0+2*ij-1 ,en(j)', s'] ] ;
                    i=[1:1:length(ky)]; j = 1:length(ky) ; s = -ky ;  
                    Cell2{c_} = [en(i)'*0+2*ij ,en(j)', s'];
                    C1 = [ C1;  [en(i)'*0+2*ij-0 ,en(j)', s'] ] ;
                    c_ = c_ + 1 ; 

                    idx = 4; 
                    id1 = ij+idx*numnode-numnode ; 
                    en = En(id1,:); en = en(en~=0); en = en(en>0);
                    kx = FX2(id1,1:length(en));
                    ky = FY2(id1,1:length(en));
                    K(2*ij-1,en) = K(2*ij-1,en) + kx ; 
                    K(2*ij-0,en) = K(2*ij-0,en) + ky ; 

                    i=[1:1:length(kx)]; j = 1:length(kx) ; s = kx ;  
                    Cell{c_} = [en(i)'*0+2*ij-1 ,en(j)', s'];
                    C1 = [ C1;  [en(i)'*0+2*ij-1 ,en(j)', s'] ] ;
                    i=[1:1:length(ky)]; j = 1:length(ky) ; s = ky ;  
                    Cell2{c_} = [en(i)'*0+2*ij ,en(j)', s'];
                    C1 = [ C1;  [en(i)'*0+2*ij-0 ,en(j)', s'] ] ;
                    c_ = c_ + 1 ; 
                end
            end
            
% %         IJV = cell2mat( Cell );
% %         K2 = K2 + sparse(IJV(:,1),IJV(:,2),IJV(:,3),size(K,1),size(K,1));
% %         IJV = cell2mat( Cell2 );
% %         K2 = K2 + sparse(IJV(:,1),IJV(:,2),IJV(:,3),size(K,1),size(K,1));
        
        K3 = K3 + sparse(C1(:,1),C1(:,2),C1(:,3),size(K,1),size(K,1));

%         max(max(abs(K-K2)))
        max(max(abs(K-K3)))
        K = K3 ; 
%         pause
%         %%
%         for ij = 1 : size(K,1) 
%             if norm(abs(K(ij,:)-K2(ij,:)))>1e-11
%                 ij
%                 norm(abs(K(ij,:))-abs(K2(ij,:)))
%                 pause
%             end
%         end
% %         K =  K2  ; 
        pause
        %}
%%
K = 0*K; 
K3 = K  ; 
stepF = 10 ;
% C1 = [] ; C2 = [] ; 
% for ij = 1 : numnode
%                     idx = 1; 
%                     id1 = ij+idx*numnode-numnode ; 
%                     en = En(id1,:); en = en(en~=0); en = en(en>0);
%                     kx = FX(id1,1:length(en));
%                     ky = FY(id1,1:length(en));
%                     K(2*ij-1,en) = K(2*ij-1,en) - kx ; 
%                     K(2*ij-0,en) = K(2*ij-0,en) - ky ;
% 
%                     en = En(id1,:);
%                     kx = FX(id1,1:length(en));
%                     ky = FY(id1,1:length(en));
%                     i=[1:1:length(kx)]; j = 1:length(kx) ; s = -kx ;  
%                     C1 = [ C1;  [en(i)'*0+2*ij-1 ,en(j)', s'] ] ;
%                     i=[1:1:length(ky)]; j = 1:length(ky) ; s = -ky ;  
%                     c_ = c_ + 1 ; 
%                     C2 = [ C2;  [en(i)'*0+2*ij-0 ,en(j)', s'] ] ;
% end
% C1
    I1 = 1:1:1*numnode;    I1=I1(ones(1,size(En,2)),:);    I1 = I1(:);
    En_ = En(1:numnode,:)';
    I2 = (En_(:)); %I2(I2==0) =1 ;
    idx_ = 1:numnode ; 
    FX_ = FX(idx_,:)';    I3 = -(FX_(:));
    C1_1 = [2*I1-1 I2 I3]; 
%     norm(C1_1-C1)
    FY_ = FY(idx_,:)';    I3 = -(FY_(:));
    C2_1 = [2*I1-0 I2 I3]; 
%     norm(C2_1-C2)
%
% C1 = [] ; C2 = [] ; 
% for ij = 1 : numnode
% 
%                     idx = 2; 
%                     id1 = ij+idx*numnode-numnode ; 
%                     en = En(id1,:); en = en(en~=0); en = en(en>0);
%                     kx = FX(id1,1:length(en));
%                     ky = FY(id1,1:length(en));
%                     K(2*ij-1,en) = K(2*ij-1,en) + kx ; 
%                     K(2*ij-0,en) = K(2*ij-0,en) + ky ; 
% 
%                     en = En(id1,:); 
%                     kx = FX(id1,1:length(en));
%                     ky = FY(id1,1:length(en));
%                     i=[1:1:length(kx)]; j = 1:length(kx) ; s = kx ;  
%                     C1 = [ C1;  [en(i)'*0+2*ij-1 ,en(j)', s'] ] ;
%                     i=[1:1:length(ky)]; j = 1:length(ky) ; s = ky ;  
%                     C2 = [ C2;  [en(i)'*0+2*ij-0 ,en(j)', s'] ] ;
%                     c_ = c_ + 1 ; 
% end
    idx_ = numnode+1:2*numnode ; 
    En_ = En(idx_,:)';
    I2 = (En_(:)); %I2(I2==0) =1 ;
    FX_ = FX(idx_,:)';    I3 = (FX_(:));
    C1_2 = [2*I1-1 I2 I3]; 
%     norm(C1_2-C1)
    FY_ = FY(idx_,:)';    I3 = (FY_(:));
    C2_2 = [2*I1-0 I2 I3]; 
%     norm(C2_2-C2)
%
% C1 = [] ; C2 = [] ; 
% for ij = 1 : numnode
%                     idx = 3; 
%                     id1 = ij+idx*numnode-numnode ; 
%                     en = En(id1,:); en = en(en~=0); en = en(en>0);
%                     kx = FX2(id1,1:length(en));
%                     ky = FY2(id1,1:length(en));
%                     K(2*ij-1,en) = K(2*ij-1,en) - kx ; 
%                     K(2*ij-0,en) = K(2*ij-0,en) - ky ; 
% 
%                     en = En(id1,:); 
%                     kx = FX2(id1,1:length(en));
%                     ky = FY2(id1,1:length(en));
%                     i=[1:1:length(kx)]; j = 1:length(kx) ; s = -kx ;  
%                     C1 = [ C1;  [en(i)'*0+2*ij-1 ,en(j)', s'] ] ;
%                     i=[1:1:length(ky)]; j = 1:length(ky) ; s = -ky ;  
%                     C2 = [ C2;  [en(i)'*0+2*ij-0 ,en(j)', s'] ] ;
%                     c_ = c_ + 1 ; 
% end
    idx_ = 2*numnode+1:3*numnode ; 
    En_ = En(idx_,:)';
    I2 = (En_(:)); %I2(I2==0) =1 ;
    FX2_ = FX2(idx_,:)';    I3 = -(FX2_(:));
    C1_3 = [2*I1-1 I2 I3]; 
%     norm(C1_3-C1)
    FY2_ = FY2(idx_,:)';    I3 = -(FY2_(:));
    C2_3 = [2*I1-0 I2 I3]; 
%     norm(C2_3-C2)
%     
%
% C1 = [] ; C2 = [] ; 
% for ij = 1 : numnode
%                     idx = 4; 
%                     id1 = ij+idx*numnode-numnode ; 
%                     en = En(id1,:); en = en(en~=0); en = en(en>0);
%                     kx = FX2(id1,1:length(en));
%                     ky = FY2(id1,1:length(en));
%                     K(2*ij-1,en) = K(2*ij-1,en) + kx ; 
%                     K(2*ij-0,en) = K(2*ij-0,en) + ky ; 
% 
%                     en = En(id1,:); 
%                     kx = FX2(id1,1:length(en));
%                     ky = FY2(id1,1:length(en));
%                     i=[1:1:length(kx)]; j = 1:length(kx) ; s = kx ;  
%                     C1 = [ C1;  [en(i)'*0+2*ij-1 ,en(j)', s'] ] ;
%                     i=[1:1:length(ky)]; j = 1:length(ky) ; s = ky ;  
%                     C2 = [ C2;  [en(i)'*0+2*ij-0 ,en(j)', s'] ] ;
%                     c_ = c_ + 1 ; 
% end
    idx_ = 3*numnode+1:4*numnode ; 
    En_ = En(idx_,:)';
    I2 = (En_(:)); %I2(I2==0) =1 ;
    FX2_ = FX2(idx_,:)';    I3 = (FX2_(:));
    C1_4 = [2*I1-1 I2 I3]; 
%     norm(C1_4-C1)
    FY2_ = FY2(idx_,:)';    I3 = (FY2_(:));
    C2_4 = [2*I1-0 I2 I3]; 
%     norm(C2_4-C2)
    
    C1 = [C1_1;C1_2;C1_3;C1_4];
    C2 = [C2_1;C2_2;C2_3;C2_4];
    
 %   
 C1(C1(:,2)<=0,:) = [] ;
 C2(C2(:,2)<=0,:) = [] ;
 C1(C1(:,1)<=0,:) = [] ;
 C2(C2(:,1)<=0,:) = [] ;

K3 = K3 + sparse(C1(:,1),C1(:,2),C1(:,3),size(K,1),size(K,1));
K3 = K3 + sparse(C2(:,1),C2(:,2),C2(:,3),size(K,1),size(K,1));
%         max(max(abs(K-K3)))
        K = K3 ; 
%%
% % tic
% % % % %     I1 = 1:1:1*numnode;    I1=I1(ones(1,size(En,2)),:);    I1 = I1(:);
% % % % %     
% % % % %     En_ = En(1:numnode,:)';
% % % % %     I2 = (En_(:)); I2(I2==0) =1 ; 
% % % % %     
% % % % %     idx_ = 1:numnode ; 
% % % % %     FX_ = FX(idx_,:)';    I3 = -(FX_(:));
% % % % %     A1 = sparse(2*I1-1,I2,I3,size(K,1),size(K,1));
% % % % % 
% % % % %     FY_ = FY(idx_,:)';    I3 = -(FY_(:));
% % % % %     A2 = sparse(2*I1-0,I2,I3,size(K,1),size(K,1));
% % % % % 
% % % % %     idx_ = numnode+1:2*numnode ; 
% % % % %     FX_ = FX(idx_,:)';    I3 = (FX_(:));
% % % % %     A3 = sparse(2*I1-1,I2,I3,size(K,1),size(K,1));
% % % % % 
% % % % %     FY_ = FY(idx_,:)';    I3 = (FY_(:));
% % % % %     A4 = sparse(2*I1-0,I2,I3,size(K,1),size(K,1));
% % % % % 
% % % % % 
% % % % %     idx_ = 2*numnode+1:3*numnode ; 
% % % % %     FX_ = FX2(idx_,:)';    I3 = -(FX_(:));
% % % % %     A5 = sparse(2*I1-1,I2,I3,size(K,1),size(K,1));
% % % % % 
% % % % %     FY_ = FY2(idx_,:)';    I3 = -(FY_(:));
% % % % %     A6 = sparse(2*I1-0,I2,I3,size(K,1),size(K,1));
% % % % % 
% % % % %     idx_ = 3*numnode+1:4*numnode ; 
% % % % %     FX_ = FX2(idx_,:)';    I3 = (FX_(:));
% % % % %     A7 = sparse(2*I1-1,I2,I3,size(K,1),size(K,1));
% % % % % 
% % % % %     FY_ = FY2(idx_,:)';    I3 = (FY_(:));
% % % % %     A8 = sparse(2*I1-0,I2,I3,size(K,1),size(K,1));
% % % % %     
% % % % %     K3 = K3 + A1 + A2  + A3  + A4  + A5  + A6  + A7 + A8;  
% % % % % % % % %     max(max(abs(K-K2)))
% % % %     max(max(abs(K2-K3)))
% % % %     pause
% %     toc 