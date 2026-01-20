function VTU_Fiber (xx,yy,f1,output)



if size(yy,2) == 2 % L2 elements
    offset = [2:2:2*size(yy,1)]';     
    type = yy(:,1)*0+3;
elseif size(yy,2) == 3 % T3 element
    offset = [3:3:3*size(yy,1)]';     
    type = yy(:,1)*0+5;
elseif size(yy,2) == 4 % Q4 element
    offset = [4:4:4*size(yy,1)]';     
%     type = yy(:,1)*0+9; % Q4
    type = yy(:,1)*0+10; % T4
elseif size(yy,2) == 8 % Q8 element
    offset = [8:8:8*size(yy,1)]';     
    type = yy(:,1)*0+12;
end

fid = fopen(output,'wt');

% preamble
fprintf(fid, '<VTKFile type="UnstructuredGrid"  version="0.1"   > \n ');
fprintf(fid, '<UnstructuredGrid> \n ');
fprintf(fid,'<Piece  NumberOfPoints=" %6.0f " NumberOfCells=" %6.0f "> \n ', size(xx,1), size(yy,1) );

% nodal data
fprintf(fid, '<Points>  \n ');
fprintf(fid, '<DataArray  type="Float64"  NumberOfComponents="3"  format="ascii" >   \n ');
fprintf(fid,'%18.15f %18.15f %18.15f \n', xx' );
fprintf(fid, '</DataArray>  \n ');
fprintf(fid, '</Points>  \n ');

% Field values at nodes
fprintf(fid, '<PointData  Vectors="fields"> \n ');
% f1
fprintf(fid, '<DataArray  type="Float64"  Name="f1" NumberOfComponents="1" format="ascii"> \n ');
fprintf(fid,'%28.25f \n', f1' );
fprintf(fid, '</DataArray>  \n ');
% % f2
% fprintf(fid, '<DataArray  type="Float64"  Name="f2" NumberOfComponents="1" format="ascii"> \n ');
% fprintf(fid,'%20.20f \n', f2' );
% fprintf(fid, '</DataArray>  \n ');
% % f3
% fprintf(fid, '<DataArray  type="Float64"  Name="f3" NumberOfComponents="1" format="ascii"> \n ');
% fprintf(fid,'%20.20f \n', f3' );
% fprintf(fid, '</DataArray>  \n ');
% % f4
% fprintf(fid, '<DataArray  type="Float64"  Name="f4" NumberOfComponents="1" format="ascii"> \n ');
% fprintf(fid,'%20.20f \n', f4' );
% fprintf(fid, '</DataArray>  \n ');
% % f5
% fprintf(fid, '<DataArray  type="Float64"  Name="f5" NumberOfComponents="1" format="ascii"> \n ');
% fprintf(fid,'%20.20f \n', f5' );
% fprintf(fid, '</DataArray>  \n ');
% % f6
% fprintf(fid, '<DataArray  type="Float64"  Name="f6" NumberOfComponents="1" format="ascii"> \n ');
% fprintf(fid,'%20.20f \n', f6' );
% fprintf(fid, '</DataArray>  \n ');
% % f7
% fprintf(fid, '<DataArray  type="Float64"  Name="f7" NumberOfComponents="1" format="ascii"> \n ');
% fprintf(fid,'%20.20f \n', f7' );
% fprintf(fid, '</DataArray>  \n ');
% % f8
% fprintf(fid, '<DataArray  type="Float64"  Name="f8" NumberOfComponents="1" format="ascii"> \n ');
% fprintf(fid,'%20.20f \n', f8' );
% fprintf(fid, '</DataArray>  \n ');

fprintf(fid, '</PointData>  \n ');

   


% connectivity data 
fprintf(fid, '<Cells>  \n ');
fprintf(fid, '<DataArray  type="Float64"  Name="connectivity"  format="ascii">  \n ');
fprintf(fid,'%8.5f %8.5f  \n', yy'-1 );
fprintf(fid, '</DataArray>  \n ');
% offsett
fprintf(fid, '<DataArray  type="Float64"  Name="offsets"  format="ascii">  \n ');
fprintf(fid,'%8.5f   \n', offset' );
fprintf(fid, '</DataArray>  \n ');
% type
fprintf(fid, '<DataArray  type="Float64"  Name="types"  format="ascii">  \n ');
fprintf(fid,'%8.5f   \n', type' );
fprintf(fid, '</DataArray>  \n ');
fprintf(fid, '</Cells>  \n ');
fprintf(fid, '</Piece>  \n ');
fprintf(fid, '</UnstructuredGrid>  \n ');
fprintf(fid, '</VTKFile>  \n ');

fclose (fid) ;