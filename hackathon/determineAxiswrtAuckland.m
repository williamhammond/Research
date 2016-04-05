function determineAxiswrtAuckland(path_coarsemesh,path_fib)

%% Read heart surface data
filename = 'heart_sur';

file = [path_coarsemesh filename '.cor'];
fid = fopen(file,'rb');
cortet = fread(fid,[3, inf],'double');
fclose(fid);

file = [path_coarsemesh filename '.idx'];
fid = fopen(file,'rb');
idx = fread(fid,'int');
fclose(fid);

file = [path_coarsemesh filename '.mat'];
surf_dat = load(file);
cor = surf_dat.scirunfield.node;
face = surf_dat.scirunfield.face;

%% Redefine z-axis if needed 
prompt= 'Wills Method? (y/n)';
str = lower(input(prompt, 's'));
z_1 = zeros(1,3);
% Top and bottom of the left ventricle. Lets us align the z-axis through 
% the center.
if (str == 'y')
    fig = figure('DeleteFcn', 'doc datacursormode');
    trimesh(face',cor(1,:), cor(2,:), cor(3,:))

    xlabel('x');
    ylabel('y');
    zlabel('z');
    
    hold on;    hold off;


    cursorobj = datacursormode(fig);
    cursorobj.SnapToDataVertex = 'on'; % Snap to our plotted data, on by default

    disp('Select apex of the heart and hit enter')
    pause
    
    mypoints = getCursorInfo(cursorobj);
    mypoints = table({mypoints.Position}.', 'VariableNames', {'Position'});
    mypoints = mypoints.Position;
    
    z_1(1) = mypoints{1}(1);
    z_1(2) = mypoints{1}(2);
    z_1(3) = mypoints{1}(3);


    disp('Select points along the outer edge of the left ventricle')
    pause
    
    % Convert result from struct to cell
    mypoints = getCursorInfo(cursorobj);
    mypoints = table({mypoints.Position}.', 'VariableNames', {'Position'});
    mypoints = mypoints.Position;

    % Store the coordinate as a matrix rather than determineAxiswrtAuckland.ma cell
    outerPoints = zeros(length(mypoints),3);

    for i=1:length(mypoints)
        outerPoints(i,1) = mypoints{i}(1);
        outerPoints(i,2) = mypoints{i}(2);
        outerPoints(i,3) = mypoints{i}(3);
    end
    
    lv_c = sum(outerPoints) / length(outerPoints);
    
    disp('Select the end of the HLA')
    pause
    
    mypoints = getCursorInfo(cursorobj);
    mypoints = table({mypoints.Position}.', 'VariableNames', {'Position'});
    mypoints = mypoints.Position;
    
    x_1(1) = mypoints{1}(1);
    x_1(2) = mypoints{1}(2);
    x_1(3) = mypoints{1}(3);
    

    %% Find Y_axis
    x_1 = x_1 - lv_c;
    z_1 = z_1 - lv_c;
    
    x_1 = x_1 ./ norm(x_1);
    z_1 = z_1 ./ norm(z_1);
    
    y_1 = cross(x_1, z_1);  
    x_1 = cross(z_1, y_1);  
       
    z_1 = [-1 0 0 ; 0 -1 0; 0 0 -1] * z_1';
    z_1 = z_1';
    A = [x_1;y_1;z_1];
    cor_n = A * (cor - (lv_c' * (ones(1,size(cor,2)))));
    
    hold off;

    trimesh(face',cor_n(1,:), cor_n(2,:), cor_n(3,:))
    pause
    
end

%% decide rotation & scale depending on surface, find base SA slice
% find the slice that has a visible ventricles
z_min = min(cortet(3,:));
z_max = max(cortet(3,:));
z_base= z_max;
numSlice=10;
step=(z_max-z_min)/numSlice;

for i = 1:numSlice
    base = find(cortet(3,:) > z_max-i*step);
    cor_base = cortet(1:3,base);
    
    % display base
    figure(1)
    plot(cor_base(1,:),cor_base(2,:),'o');
    xlabel('x');
    ylabel('y');         
    
    prompt = 'Select Base in this slice? (y/n) ';
    str = lower(input(prompt, 's'));
    if (str=='y')        
        % % find the origin of the base
        c_bo = zeros(3,1); % base origin              
        input('Press enter and select the center of the base in figure');
        [c_bo(1), c_bo(2)] = ginput(1);
        c_bo(3) = z_base;
        break;
    end
end

%interaction to decide horizontal long axis, vertical long axis in local x-y plane
cv = zeros(3,1);
input('select points that extends from base to horizontal long axis of heart'); 
[cv(1), cv(2)] = ginput(1);
cv(3) = z_base;
hLA = cv - c_bo; % Why do we subtract here
hLA = hLA./norm(hLA);

%compute vertical Long Axis                   
key = input('rotate (1) 90 or (2) -90 around z axis(get vertical long axis from horizontal long axis)?');
if key == 1
    A = [0 -1 0; 1 0 0; 0 0 1];
else
    A = [0 1 0; -1 0 0; 0 0 1]; % select
end

vLA = A*hLA; 



% display
hold on;
quiver([c_bo(1),c_bo(1)],[c_bo(2),c_bo(2)],50.*[vLA(1),hLA(1)],50.*[vLA(2),hLA(2)]);

filename = 'heart_mfnode';
key = input('Save rotation matrix? (y/n)','s');
if key == 'y'
    file = strcat(path_fib,filename);
    file = strcat(file,'_4fibermapping.rot');
    fid = fopen(file,'wb');
    fwrite(fid,horzcat(hLA,vLA),'double');
    fclose(fid);
end

A(1,:) = hLA';
A(2,:) = vLA';
A(3,:) = cross(A(1,:),A(2,:));

%% Rotate heart to auckland's heart coordinate system
cor_n = A*(cor - c_bo*ones(1,size(cor,2)));
cor_ntet = A*(cortet - c_bo*ones(1,size(cortet,2)));
cor_n(3,:) = cor_n(3,:)-15; % why subtract 15

hold off;
figure;
hold on;
plot3(cor_n(1,:),cor_n(2,:),cor_n(3,:),'b+');
xlabel('x');
ylabel('y');
title('after conversion to coordinates of Auckland heart')
path_auk='/home/wth4280/Documents/Research/data/Auckland/'
file = [path_auk 'Auckland_f02_2_5.cor'];
fid = fopen(file,'rb');
aukcor = fread(fid,[3,inf],'double');
fclose(fid);
plot3(aukcor(1,:),aukcor(2,:),aukcor(3,:),'r*');
% cor_n(3,:) = cor_n(3,:)-15;
% cor(3,:) = cor(3,:)-15;
%cor_ntet(3,:) = cor_ntet(3,:);

%% Write data
file = strcat(path_fib,'heart_mfnode');
file = strcat(file,'.cor');
fid = fopen(file,'wb');
fwrite(fid,cor,'double');
fclose(fid);


file = strcat(path_fib,filename);
file = strcat(file,'_4fibermapping.cor');
fid = fopen(file,'wb');
fwrite(fid,cor_n,'double');
fclose(fid);

file = strcat(path_fib,filename);
file = strcat(file,'_4fibermapping.idx');
fid = fopen(file,'wb');
fwrite(fid,idx,'int');
fclose(fid);

file = strcat(path_fib,filename);
file = strcat(file,'_4fibermapping_tet.cor');
fid = fopen(file,'wb');
fwrite(fid,cor_ntet,'double');
fclose(fid);

fibA = A;
savefile1 = strcat(path_fib,'fibmatrix.mat');
save(savefile1,'fibA');

% file = strcat(path_sur,'heart_mfnode_azar');
% file = strcat(file,'_4fibermapping_lv_tet.cor');
% fid = fopen(file,'wb');
% fwrite(fid,lv_n,'double')
% fclose(fid);



    