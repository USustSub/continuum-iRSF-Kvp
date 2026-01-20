
% filename = ['in/' Ex '_msh' num2str(nnx) '.mat' ] ; 
filename = ['../Run_new/in/'  'Herr_msh' num2str(nnx) '.mat' ] ; 

if isfile(filename) == 0 || rewrite_mesh == 1 
% mesh has not been generated before
[HH,mm_,dY,DD,DD_BC,cc_ij,B_x_b,B_y_b,B_xy_b,B_xy_r,B_y_r,B_x_r,...
B_x_t,B_y_t,B_xy_t,B_x_l,B_y_l,B_xy_l,...
Dmat_t1_l,Dmat_t1_r,Dmat_t1_b,Dmat_t1_t,...
D11,D12,D13,D14,D21,D22,D23,D24,...
D31,D32,D33,D34,D41,D42,D43,D44,...
maxS,SS,SS0,id1,Bx,By,En,Shape,Gve_,K__,dcX , dcY,dcX_L,dcX_R,dcY_B,dcY_T] = Pre_process_new_efficient_v3 () ; % get shape functions, ... 

% eval ( ['save ' filename ' HH mm_ dY DD DD_BC cc_ij B_x_b B_y_b B_xy_b B_xy_r B_y_r B_x_r B_x_t B_y_t B_xy_t B_x_l B_y_l B_xy_l Dmat_t1_l Dmat_t1_r Dmat_t1_b Dmat_t1_t D11 D12 D13 D14 D21 D22 D23 D24 D31 D32 D33 D34 D41 D42 D43 D44 maxS SS SS0 id1 Bx By En Shape Gve_ K__ dcX dcY dcX_L dcX_R dcY_B dcY_T'] ) ; 

else
% mesh already exists, load it
eval ( ['load ' filename ' HH mm_ dY DD DD_BC cc_ij B_x_b B_y_b B_xy_b B_xy_r B_y_r B_x_r B_x_t B_y_t B_xy_t B_x_l B_y_l B_xy_l Dmat_t1_l Dmat_t1_r Dmat_t1_b Dmat_t1_t D11 D12 D13 D14 D21 D22 D23 D24 D31 D32 D33 D34 D41 D42 D43 D44 maxS SS SS0 id1 Bx By En Shape Gve_ K__ dcX dcY dcX_L dcX_R dcY_B dcY_T'] ) ; 

if strcmp(Ex,'Preuss')
    SS = SS/1.5 ; 
    SS0 = SS0/1.5 ; 
end

end
node_dist = [deltaX deltaX deltaX deltaX ]; 

get_faults ; 


