% load  here.mat


if size(find(dg_c(inx(:))<0),1) > 0
    warning (['Bisection needed for ' num2str(size(find(dg_c(inx(:))<0),1)) ' nodes'])
if inc < 3

dg_c_n = dg_c(inx);
[oo_ , ~] = find(dg_c_n<0);
% dg_c_n = dg_c_n(oo_);
        
dg_a = 0 ;
dg_b = 0.5 ; 

dg_a_n = dg_c_n*0+dg_a; 
dg_b_n = dg_c_n*0+dg_b; 
get_Fa_Fb_n  ;      
dg_p_n = (dg_a_n + dg_b_n )/2; 
get_Fp_n ;
err = abs(Fp_n(oo_));
    while any(err > 1e5)
        f_t = Fa_n.*Fb_n;
        dg_b_n(f_t<0) = dg_p_n(f_t<0) ; 
        dg_a_n(f_t>0) = dg_p_n(f_t>0) ;
        dg_p_n = (dg_a_n + dg_b_n)/2; 
        get_Fp_n ;
        err = abs(Fp_n(oo_)); 
%         dg_b_n(1)
%         dg_a_n(1)
%             (err(1))
%         pause
    end
    dg_c(inx(oo_)) = dg_p_n(oo_);
    sin_phi_(inx(oo_)) = sin_phi_p_n(oo_);
    Vp(inx(oo_)) = Vp_p_n(oo_) ; 
    theta_n(inx(oo_)) = theta_new_p_n(oo_);

else
    
% %{
% load here.mat
for iii_  = 1 : length(inx)
% dg_c(inx(iii_))
% dg_c(inx(iii_))
    if dg_c(inx(iii_))<0

%     iii_/length(inx)
%         disp('Bisection correction due to negative root') 
        dg_a = 0 ;
        dg_b = 0.5 ; 
        get_Fa_Fb
        dg_p = (dg_a + dg_b)/2;
        get_Fp ;

        err = abs(Fp);
       while err > 1e-9
       if Fa*Fp<0 
           dg_b = dg_p;
       else
           dg_a = dg_p;          
       end
%         dg_a
%         dg_b
        
        dg_p = (dg_a + dg_b)/2; 
        get_Fp ;
       err = abs(Fp);
%        err
%        pause
       end
        dg_c(inx(iii_)) = dg_p;
        sin_phi_(inx(iii_)) = sin_phi_p;
        Vp(inx(iii_)) = Vp_p ; 
        theta_n(inx(iii_)) = theta_new_p;
    end
end


%}
end

end
%%
