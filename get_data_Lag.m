function [phi,B,en] = get_data_Lag ( x_l , con_cand , node , element , side ) 
        

        if strcmp(side,'L') || strcmp(side,'R')
            con_cand_ = [] ; cc =  1; 
            for ik = 1 : length(con_cand)
                sc_ik = element(con_cand(ik),:) ;
                x_diff = x_l-node(sc_ik,:) ; 
                x_diff(abs(x_diff)<1e-5)=0;
                if nnz(x_diff(:,2)) == 2 && length(unique((sign(x_diff(:,1)))))==2
                    con_cand_(cc) = con_cand(ik) ; 
                    cc = cc +  1 ; 
                end
            end
        end
        if strcmp(side,'B') || strcmp(side,'T')
            con_cand_ = [] ; cc =  1; 
            for ik = 1 : length(con_cand)
                sc_ik = element(con_cand(ik),:) ;
                x_diff = x_l-node(sc_ik,:) ;
                x_diff(abs(x_diff)<1e-5)=0;
                if nnz(x_diff(:,1)) == 2 && length(unique((sign(x_diff(:,2)))))==2
                    con_cand_(cc) = con_cand(ik) ; 
                    cc = cc +  1 ; 
                end
            end
        end
        if length(con_cand_)>2
%             error('sss')
% con_cand_
        end
        sctr = element(con_cand_(1),:) ;
        local = element_coordinate ( x_l , sctr , node );

        [N,dNdxi] = lagrange_basis('Q4',local);  % element shape functions
        J0 = node(sctr,:)'*dNdxi(:,:);                 % element Jacobian matrix
        invJ0 = inv(J0);
        dNdx  = dNdxi*invJ0;                      % derivatives of N w.r.t XY

        B1 = [] ; nn = 4 ; 
        B1(1,1:2:2*nn)      = dNdx(:,1)';
        B1(2,2:2:2*nn) = dNdx(:,2)';
        B1(3,1:2:2*nn)      = dNdx(:,2)';
        B1(3,2:2:2*nn) = dNdx(:,1)';

        sctrB1 = zeros(1,length(sctr)*2) ;
        sctrB1(1:2:end) = 2*sctr-1 ; 
        sctrB1(2:2:end) = 2*sctr ;

%         sctr = element(con_cand_(2),:) ; 
%         local = element_coordinate ( x_l , sctr , node );
% 
%         [N,dNdxi] = lagrange_basis('Q4',local);  % element shape functions
%         J0 = node(sctr,:)'*dNdxi(:,:);                 % element Jacobian matrix
%         xx = N(1:4)'*node(sctr,:) ;
%         invJ0 = inv(J0);
%         dNdx  = dNdxi*invJ0;                      % derivatives of N w.r.t XY
% 
%         B2 = [] ; nn = 4 ; 
%         B2(1,1:2:2*nn)      = dNdx(:,1)';
%         B2(2,2:2:2*nn) = dNdx(:,2)';
%         B2(3,1:2:2*nn)      = dNdx(:,2)';
%         B2(3,2:2:2*nn) = dNdx(:,1)';
% 
%         sctrB2 = zeros(1,length(sctr)*2) ;
%         sctrB2(1:2:end) = 2*sctr-1 ; 
%         sctrB2(2:2:end) = 2*sctr ;

        en = [sctrB1];% sctrB2] ; 
        B = [B1/1];% B2/2]; 
        phi = [N' ] ; 