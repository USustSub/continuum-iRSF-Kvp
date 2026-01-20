


vtfile =['out/paraview/dUx' num2str(inc) '.vtu'] ;
VTU_Fiber ([node node(:,1)*0],element,dU(1:2:end),vtfile);
vtfile =['out/paraview/dUy' num2str(inc) '.vtu'] ;
VTU_Fiber ([node node(:,1)*0],element,dU(2:2:end),vtfile); 
vtfile =['out/paraview/vx' num2str(inc) '.vtu'] ;
VTU_Fiber ([node node(:,1)*0],element,vu(1:2:end),vtfile); 
vtfile =['out/paraview/vy' num2str(inc) '.vtu'] ;
VTU_Fiber ([node node(:,1)*0],element,vu(2:2:end),vtfile); 
vtfile =['out/paraview/ax' num2str(inc) '.vtu'] ;
VTU_Fiber ([node node(:,1)*0],element,au(1:2:end),vtfile); 
vtfile =['out/paraview/ay' num2str(inc) '.vtu'] ;
VTU_Fiber ([node node(:,1)*0],element,au(2:2:end),vtfile); 

% VTKPostProcess(SS(:,1:2),tri,1,'Tri3',['stress/stress' num2str(inc)],SS(:,3:5),[SS(:,6) SS(:,7)])

vtfile =['out/paraview/Eii' num2str(inc) '.vtu'] ;
VTU_Fiber ([SS(:,1:2) SS(:,1)*0],tri,SS(:,6),vtfile);

vtfile =['out/paraview/SQJ' num2str(inc) '.vtu'] ;
VTU_Fiber ([SS(:,1:2) SS(:,1)*0],tri,SS(:,7),vtfile); 
