function write_vtk1(option,value,idx, varargin)
% option = 'v', tetra; 's', surface
% value = the vector to be written 1xN.
% varargin: idx to indicate epi vs. endo vs. mfree node. Used when trying
% to write surface results from mfree outcome; 
% for example, varargin{1} = find(idx==3); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HOW TO RUN THIS FUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ===>>>>   write_vtk('s',t')====>>> t HAS TO BE TRANSPOSED COZ IT WAS Nx1, so it has
%                          to be 1xN. HERE IN THIS EXAMPLE MY S IS
%                          (413x971) so the t = 413x1 == >> after
%                          transposed is t' = 1x413.
% Input the name for the attribute
% t "This is just a name"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin > 4
    idx_mfree = varargin{1};
end

%heart struct location
path = '/home/wth4280/Documents/Research/data/Hackathon/geometry';
path_out = '/home/wth4280/Documents/Research/VTK/';
% load(strcat(path, 'Heart.mat'));
% cor = heart.node;
% [a,b]=size(cor)
% tet = heart.face;
% 

fid=fopen('heart_tet.cor','rb');
cor=fread(fid,[3,inf], 'double');
fclose(fid);

fid=fopen('heart_tet.tet','rb');
tet=fread(fid,[4,inf], 'double');
fclose(fid);

tet = tet - 1;

num_pts = size(cor,2);
num_tet = size(tet,2);

file_out = strcat(path_out,'Electric field');
file_out = strcat(file_out,num2str(idx),'.vtk');
fidw = fopen(file_out,'wt');

%header
fprintf(fidw, '# vtk DataFile Version 4.0\n');
fprintf(fidw, 'vtk paraview fie\n');
fprintf(fidw, 'ASCII\n');

%geometry
if option == 's'
    fprintf(fidw, 'DATASET POLYDATA\n');
else
    fprintf(fidw, 'DATASET UNSTRUCTURED_GRID\n');
end

fprintf(fidw, 'POINTS %d float\n', num_pts);
fprintf(fidw, '%f %f %f\n', cor);
fprintf(fidw, '\n');

if option == 's'
    fprintf(fidw, 'TRIANGLE_STRIPS %d %d\n', num_tet, 4*num_tet);
    fprintf(fidw, '%d %d %d %d\n', vertcat(3*ones(1,num_tet),tet));
else
    fprintf(fidw, 'CELLS %d %d\n', num_tet, 5*num_tet);
    fprintf(fidw, '%d %d %d %d %d\n', vertcat(4*ones(1,num_tet), tet));
    fprintf(fidw, 'CELL_TYPES %d\n', num_tet);
    fprintf(fidw, '%d\n', 10*ones(num_tet,1));
end

if isempty(value) == 1
    key = input('Add attribute data? (y/n) \n');
    if key == 'y'
        fprintf(fidw, 'POINT_DATA %d\n', num_pts);
        num_field = input('Input the number of attribute data to display\n');
        fprintf(fidw, 'FIELD attributedata %d\n', num_field);
        for i=1:num_field
            name_att = input('Input the filename for the attribute (including extension)\n','s');
            file = strcat(path,name_att);
            fid = fopen(file,'rb');
            val = fread(fid,'double');
            fclose(fid);
            
            if nargin > 4
                val = val(idx_mfree);
            end
            
            name_att = input('Input the name for the attribute\n', 's');
            fprintf(fidw, '%s %d %d float\n', name_att, 1, num_pts);
            fprintf(fidw, '%f\n', val);
            fprintf(fidw, '\n');
        end
        
        key1 = input('Add vector data?(y/n)\n');
        if key1 == 'y'
            num_vector = input('Input the number of vector data to display\n');
            for i=1:num_vector
                name_att = input('Input the filename for the attribute (including extension)\n','s');
                file = strcat(path,name_att);
                fid = fopen(file,'rb');
                val = fread(fid,[3,inf],'double');
                if nargin>4
                    val = val(:,idx_mfree);
                end
                fclose(fid);
                name_att = input('Input the name\n','s');
                
                if size(val,2) == num_tet %element dta
                    fprintf(fidw,'CELL_DATA %d\n', num_tet);
%                 else
%                     fprintf(fidw,'POINT_DATA %d\n', num_pts);
                 end
                fprintf(fidw,'VECTORS %s float\n',name_att);
                fprintf(fidw, '%f %f %f\n', val);
                fprintf(fidw, '\n');
            end
        end
    end
else
    fprintf(fidw, '\nPOINT_DATA %d\n', num_pts);
    fprintf(fidw, 'FIELD attributedata %d\n', 1);
    name_att = input('Input the name for the attribute\n', 's');
    fprintf(fidw, '%s %d %d float\n', name_att, 1, num_pts);
    fprintf(fidw, '%f\n', value);
    fprintf(fidw, '\n');
end



fclose(fidw);