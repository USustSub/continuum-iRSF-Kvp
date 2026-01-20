%% get fault properties
disp('evaluating faults')
d1 = log(L_/V0*exp(-1)) ;
d2 = log(L_/V0*exp(40)) ; 

m_type = zeros(4*numnode,2) ; 
T_type = zeros(4*numnode,1) ; 
m_type(:,2) = 1 ; 
for ij = 1 : numnode % loop over nodes 
    if rem(ij,numnode/10)==0
        disp(['ij -- > ' num2str(ij/numnode*100) ' %'])
    end
    node_c = node(ij,:) ;
    
    x_l = node_c - [dcX_L(ij) 0 ] ;
    x_r = node_c + [dcX_R(ij) 0 ] ;

    x_b = node_c - [0 dcY_B(ij) ] ;
    x_t = node_c + [0 dcY_T(ij) ] ;

    if (node_c(1) > 0 && node_c(1) < L && node_c(2) < D && node_c(2) > 0 )
        for ii = 1 : 4
            if ii == 1
                xx = x_l ; 
            elseif ii == 2
                xx = x_r ;
            elseif ii == 3
                xx = x_b ;
            elseif ii == 4
                xx = x_t ;
            end
            
            
            if strcmp(Ex,'Slip') == 1 || strcmp(Ex,'Herr') == 1 
                if norm(abs(xx(2)-D/2))<=fault_width/1+0.01*dY % && xx(1) > 0.2 && xx(1)<0.8
            %                 if norm(abs(xx(2)-D/2))<=L/100+0.01*dY % && xx(1) > 0.2 && xx(1)<0.8
                    DD = [ DD ; xx ij+(ii-1)*numnode] ; 
                    C0 = 0.00 ; 
                    HH(ij,ii*6-0) = C0 ; 
                end
            end
        
            if strcmp(Ex,'Preuss')
                rx = xx(1)-L/2;
                ry = xx(2)-D/2;
                CL = 2000 ; CL2 = CL/5; 
                ecc = rx^2/CL/CL + ry^2/CL2/CL2;
                if ecc <= 1
                    DD = [ DD ; xx ij+(ii-1)*numnode] ;
                    C0 = 0.00 ;
                    HH(ij,ii*6-0) = C0 ;
                end
%                 rx = xx(1)-L/2;
%                 ry = xx(2)-D/2-D/4;
%                 CL = 1000 ; CL2 = CL/5; 
%                 ecc = rx^2/CL/CL + ry^2/CL2/CL2;
%                 if ecc <= 1
%                     DD = [ DD ; xx ij+(ii-1)*numnode] ;
%                     C0 = 0.00 ;
%                     HH(ij,ii*6-0) = C0 ;
%                 end
            end
            
            
        if strcmp(Ex,'Junction') == 1  || strcmp ( Ex,'Preuss2') == 1 
            for ii___ = 1:n_faults% 4 10 18]
                counter3__
                for ij_ = 1 : 1 : size(lines_,1)/1-1
                    val_ = fault_width/1 ; 
        %                         val_ = dY/2;
                    v1 = lines_(ij_,:) ; 
                    v2 = lines_(ij_+1,:) ;
                    center_line = mean(lines_);
                    length_line = norm(lines_(2,:)-lines_(1,:));
                    Length_ = norm(v2-v1);
                    d2D = point_to_line_distance(xx, v1, v2);
                    
                    d2D_2 = point_to_line(xx, v1, v2);
                    
%                     if norm(d2D-d2D_2)>1e-8
%                         pause
%                     end
                    d_t1 = norm(v1-xx);     d_t2 = norm(v2-xx); 
                    T_type(ij+(ii-1)*numnode,1) = exp(d2) ;  
                    if d2D < val_ 
                        if  d_t1^2+d_t2^2-2*d2D^2<Length_^2
                            DD = [ DD ; xx ij+(ii-1)*numnode] ;
                            C0 = 0.00 ;
                            HH(ij,ii*6-0) = C0 ;
                            if n_faults~=-1
                                c_dist_ = norm(xx-center_line) ; 
                                if c_dist_<.3*length_line
                                    m_type(ij+(ii-1)*numnode,1) = 1 ;     
                                else
                                    m_type(ij+(ii-1)*numnode,1) = 2 ;     
                                end
                            end
                            m_type(ij+(ii-1)*numnode,2) = (1-d2D/val_) ;
                            T_type(ij+(ii-1)*numnode,1) = exp((val_*d1+d2D*(d2-d1))/val_) ; 
                        end
                    end
                end
            end
        end
        
        end
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
