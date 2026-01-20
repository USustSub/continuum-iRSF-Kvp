% if Time/60/60/24/365 > 35 && h__ == 1
%        dt = 0.001 ;
maxVp = max(Vp);
dt0 = dt ; 

if inc > 1
    k_ = 2/pi*Gve_(inx(1))/deltaY;
    P__e = (1-lambda__)*P__;
    zeta = 1/4 * ( k_*L_/a_./P__e - (b_-a_)/a_ ).^2 - k_*L_/a_./P__e ; 

    dtheta_max = 0*zeta ; 
    dtheta_max1 = a_*P__e./(k_*L_-(b_-a_).*P__e);
    dtheta_max2 = 1 - (b_-a_).*P__e/(k_*L_);
    dtheta_max(zeta>0) = dtheta_max1(zeta>0);
    dtheta_max(zeta<0) = dtheta_max1(zeta<0);
    dt_max = min(dtheta_max) ; 
    dt_max = 0.2*5/2 ; 
    dt_max = 0.25 ; 
%     if maxVp < 1e-4
%         dt_max = 0.05 ; 
%     end
    if maxVp < 1e-3 %&& inc > 1500
        dt_max = 0.5 ; 
    end
%     if maxVp < 1e-6
%         dt_max = 0.000005 ; 
%     end
%     if maxVp < 1e-7
%         dt_max = 0.00000005 ; 
%     end
else
        dt_max  = 0.2 ; 
end


% adjust dt for after seismic 
if inc > 5000 && inc < 5050
    dt_max = 1 ; 
%     if LD(end,2) > LD(end-1,2) && after_seismic == 1
%         dt = 1000 ; 
%         after_seismic = 0 ;
%         warning  ('entered after seismic')
%     end
end

dt = min([dt_max*L_/maxVp , 5*1e6*4]) ;

% dt = min(dt,1);

nDT = dt/dt0 ; 



disp(['Time step is ' num2str(dt) ' sec ' num2str(dt/60/60/24/365) ' years ' ' max Vp is ' num2str(maxVp) ])
disp(['nDT is ---> ' num2str(nDT)  '  dtheta_max is ---> ' num2str(dt_max) ])
%% Newmark constants:
beta = 2;
gama = 1.5 ; 
theta = 2;
if (beta>=(0.25*(0.5+gama)^2))&&(theta>=0.5)&&(gama>=0.5)
%     disp('Newmark Constants are OK!')
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

% h__ = 0 ; 
% end
