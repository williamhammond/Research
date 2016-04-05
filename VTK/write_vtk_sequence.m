function u = write_vtk_sequence(option,value,interval,outputfile,varargin)
% option = 'v', tetra; 's', surface
% value: data matrix to be visualized; if [], read it from files -- will prompt for the entire file path including filenames; varargin{1} specifies which file to read (use only when TMP is saved in multile files)
% interval: the time interal at which the file is saved
% varargin: fileid, used in cases where "value" needs to be read from a set
% of single files
path='/home/wth4280/Documents/Research/VTK/';
addpath('/home/wth4280/Documents/Research/data/Hackathon/geometry');
if nargin>4
    fileid = varargin{1};
    
    value = [];
    path_in = input('Input path to input files \n','s');
    for i=1:length(fileid)
        %read tmp
        file_in = strcat(path_in,num2str(fileid(i)));
        file_in = strcat(file_in,'.bin');
        fid = fopen(file_in,'rb');
        value = [value fread(fid,'double')];
        u = value;
        fclose(fid);
    end
end

%name = input('Input file name for geometry (without extension) \n','s');

% file = strcat(path,name);
% file = strcat(file,'.cor');
% fid = fopen(file,'rb');
% cor = fread(fid,[3,inf],'double');
% fclose(fid);
% 
% if option == 'v'
%     file = strcat(path,name);
%     file = strcat(file,'.tet');
%     fid = fopen(file,'rb');
%     tet = fread(fid,[4,inf],'int');
%     fclose(fid);
% else
%     file = strcat(path,name);
%     file = strcat(file,'.tri');
%     fid = fopen(file,'rb');
%     tet = fread(fid,[3,inf],'int');
%     fclose(fid);
% end
fid=fopen('heart_tet.cor','rb');
cor=fread(fid,[3,inf], 'double');
fclose(fid);

fid=fopen('heart_tet.tet','rb');
tet=fread(fid,[4,inf], 'double');
fclose(fid);

tet = tet - 1;

num_pts = size(cor,2);
num_tet = size(tet,2);

%interpolate if needed
if nargin>4
    file = strcat(path,'Inter.bin');
   % file = strcat(file,num2str(num_pts));
    fid = fopen(file,'rb');
    T = fread(fid,[num_pts,inf],'double');
    fclose(fid);
    value = T*value;
end

% generate one vtk file for each time step
num_time = size(value, 2);

name_att = input('Input the name for the attribute\n', 's');

for i=1:interval:num_time
    file_out = strcat(path,outputfile);
    file_out = strcat(file_out,'_');
    if nargin > 5
        file_out = strcat(file_out, num2str(fileid(i)));
    else
        file_out = strcat(file_out, num2str(i));
    end
    file_out = strcat(file_out, '.vtk');
    
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
    
    fprintf(fidw, '\n');
    
    %scalar
    fprintf(fidw, 'POINT_DATA %d\n', num_pts);
    %fprintf(fidw, 'FIELD attributedata %d\n', 1);
    %fprintf(fidw, '%s %d %d float\n', name_att, 1, num_pts);
    %fprintf(fidw, '%f\n', value);
    %fprintf(fidw, '\n'); 
    fprintf(fidw, 'SCALARS time_sequence float %d\n', 1);
    fprintf(fidw, 'LOOKUP_TABLE default\n');
     fprintf(fidw, '%f\n', value(:,i));
%     fprintf(fidw, 'LOOKUP_TABLE mytable %d\n', 10);
%     fprintf(fidw, '%f %f %f %f\n', 0, 0, 1, 1);

%     fprintf(fidw, '%f %f %f %f\n', 0, 0.33, 1, 1);
%     fprintf(fidw, '%f %f %f %f\n', 0, 0.67, 1, 1);
%     fprintf(fidw, '%f %f %f %f\n', 1, 1, 1, 1);
%     fprintf(fidw, '%f %f %f %f\n', 0.33, 1, 0.67, 1);
%     fprintf(fidw, '%f %f %f %f\n', 0.67, 1, 0.33, 1);
%     fprintf(fidw, '%f %f %f %f\n', 1, 1, 0, 1);
%     fprintf(fidw, '%f %f %f %f\n', 1, 0.67, 0, 1);
%     fprintf(fidw, '%f %f %f %f\n', 1, 0.33, 0, 1);
%     fprintf(fidw, '%f %f %f %f\n', 1, 0, 0, 1);
%     
    fclose(fidw);
end
