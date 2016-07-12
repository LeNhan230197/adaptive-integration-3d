% Initialization
clear all;

xmin = 0;
xmax = 1;
ymin = 0;
ymax = 1;
zmin = 0;
zmax = 1;
mesh_size = [5 5 5]; % Initial meshgrid size
n_recursion_max = 2;

% Start profiling
profile on;
tic;

% Target function @f is defined is defined in a separate script, please 
% refer to f.m for more information. The function returns a list of vectors
% and a list of values corresponding to the vectors. These should be enough
% to get an integral value
[keys, vals] = adaptive_search_3d(@f, xmin, xmax, ymin, ymax, zmin, zmax, mesh_size, n_recursion_max);
combined_data = [keys, vals];

toc;
profile viewer;
profile off;

% Plotting
figure;
% plot(keys(:, 1), keys(:, 2), 'r.');
scatter3(keys(:, 1), keys(:, 2), keys(:, 3), 'filled');

% A less precise integration method

unit_size = max(diff(keys(:, 1))) * max(diff(keys(:, 2))) * max(diff(keys(:, 3)));
integral = 0;
for i = 1 : length(vals)
    key = keys(i, :);
    val = vals(i);
    
    % Only worry about the surface as the edge and corner are too small a
    % subset to worry at this stage
    if key(1) == xmin || key(2) == ymin || key(1) == xmax || key(2) == ymax || key(3) == zmin || key(3) == zmax
    	integral = integral + val * unit_size / 2;
    else
        integral = integral + val * unit_size;
    end
end

display(integral);


% % A demo of integration (2D): It can be made more precise.
% unit_size = max(diff(keys(:, 1))) * max(diff(keys(:, 2))) * max(diff(keys(:, 3)));
% 
% integral = 0;
% 
% for i = 1 : length(vals)
%     key = keys(i, :);
%     val = vals(i);
%     
%     % Corners are shared by 4 boxes
%     if (key(1) == xmin && key(2) == ymin) || (key(1) == xmax && key(2) == ymin) || (key(1) == xmin && key(2) == ymax) || (key(1) == xmax && key(2) == ymax)
%         integral = integral + val * unit_size / 4;
%         
%     % Edges are shared by 2 boxes
%     elseif key(1) == xmin || key(2) == ymin || key(1) == xmax || key(2) == ymax
%         integral = integral + val * unit_size / 2;
%     
%     % Rest are shared by 1 box
%     else
%         integral = integral + val * unit_size;
%     end
% end
% 
% display(integral);
    
    
    
    
    
    
    
