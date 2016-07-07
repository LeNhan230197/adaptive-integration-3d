% Update in V0.4: 
% Enable support in Matlab

% Define a wrapper function which Matlab requires
function [] = benchmark()

% Target function
function ret = f(x, y)
	% Detailed implementation of the target function goes here
    %pause(5);

    % ret=x.^2 + y.^2;
    % ret(ret<1)=0;
    % ret(ret>1)=1;s
    if x^2 + y^2 > 1
        ret = 1;
    else
        ret = 0;
    end 
end

function [X Y] = coord_meshgrid(xmin, xmax, ymin, ymax, mesh_size)

    % Update in V0.4: Matlab doesn't support underscore _, it is replaced
    % by t to denote temporary
	tx = linspace(xmin, xmax, mesh_size(1));
	ty = linspace(ymin, ymax, mesh_size(2));
	[X Y] = meshgrid(tx, ty);
end

% Function to convert meshgrid of coordinates to cellgrid of coordinates
function coord_cellgrid = mesh2cell(X, Y)
	mapper = @(x, y) [x y]; % Mapper function to convert coords to vector 
	coord_cellgrid = arrayfun(mapper, X, Y, 'un', 0); % Tested on GNU Octave
end

% Function to convert cellgrid of coordinates to meshgrid of coordinates
function [X Y] = cell2mesh(coord_cellgrid)
	mesh_size = size(coord_cellgrid); % Get the size of mesh
	
	% Generate index mesh to retrieve data from cellgrid
    % Update in V0.4: Matlab doesn't support underscore _, it is replaced
    % by t to denote temporary    
	ta = 1 : mesh_size(1); 
	tb = 1 : mesh_size(2); 
	[tA, tB] = meshgrid(ta, tb);
	
	% Mapper function to retrieve x and y respectively
	inverse_mapper_x = @(a, b) coord_cellgrid{a, b}(1); % Inverse map for X
	inverse_mapper_y = @(a, b) coord_cellgrid{a, b}(2); % Inverse map for Y
	
	% Apply mapper function to collectively retrieve X and Y from cellgrid
	X = arrayfun(inverse_mapper_x, tB, tA); % _B and _A has to be in this order
	Y = arrayfun(inverse_mapper_y, tB, tA); % to correctly retrieve X and Y	
end

% Function to convert cell to array
function arr = cell2array(coord_cellgrid)
	mesh_size = size(coord_cellgrid);
	
	% Reshape the cell to one colomn 
    % Update in V0.4: Matlab doesn't support underscore _, it is replaced
    % by t to denote temporary   
	tlength = mesh_size(1) * mesh_size(2); % Calculate needed array length
	cell_reshaped = reshape(coord_cellgrid, tlength, 1); % Built-in function
	
	% Convert cell to array
	arr = cell2mat(cell_reshaped);	% Built-in function to convert cell to matrix
end

% Initialization
xmin = 0;
xmax = 1;
ymin = 0;
ymax = 1;
mesh_size = [15 15]; % Initial meshgrid size
unit_size = (xmax - xmin) * (ymax - ymin) / (mesh_size(1) * mesh_size(2));

[X Y] = coord_meshgrid(xmin, xmax, ymin, ymax, mesh_size);

tic;
sum(sum(arrayfun(@f, X, Y) * unit_size, 1), 2)
toc;

% DO NOT USE weights as variable name
% all_weights = get_all_weights();
% all_weights.keys()
% all_weights.values()
%plot(data(:, 1), data(:, 2), 'r.');

end %end of benchmark function
