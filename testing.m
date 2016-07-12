% testing.m
% Script for unit testing

function [] = testing()

[X, Y, Z] = coord_meshgrid(0, 1, 0, 1, 0, 1, [2, 3, 4]);
% display(X);
% display(Y);
% display(Z);

cells = mesh2cell(X, Y, Z);
% display(cells);

[nX, nY, nZ] = cell2mesh(cells);

assert(isequal(X, nX));
assert(isequal(Y, nY));
assert(isequal(Z, nZ));

d2cells = subcell(cells,3);

display(d2cells);

function [X, Y, Z] = coord_meshgrid(xmin, xmax, ymin, ymax, zmin, zmax, mesh_size)
	tx = linspace(xmin, xmax, mesh_size(1));
	ty = linspace(ymin, ymax, mesh_size(2));
    tz = linspace(zmin, zmax, mesh_size(3));
	[X, Y, Z] = meshgrid(tx, ty, tz);
end % coord_meshgrid

% Function to convert meshgrid of coordinates to cellgrid of coordinates
function coord_cellgrid = mesh2cell(X, Y, Z)
	mapper = @(x, y, z) [x y z]; % Mapper function to convert coords to vector 
	coord_cellgrid = arrayfun(mapper, X, Y, Z, 'un', 0); 
end

function [X, Y, Z] = cell2mesh(coord_cellgrid)
	t_mesh_size = size(coord_cellgrid); % Get the size of mesh
	
	% Generate index mesh to retrieve data from cellgrid
    % Update in V0.4: Matlab doesn't support underscore _, it is replaced
    % by t to denote temporary    
	ta = 1 : t_mesh_size(1); 
	tb = 1 : t_mesh_size(2); 
    tc = 1 : t_mesh_size(3);
	[tA, tB, tC] = meshgrid(ta, tb, tc);
	
	% Mapper function to retrieve x and y respectively
	inverse_mapper_x = @(a, b, c) coord_cellgrid{a, b, c}(1); % Inverse map for X
	inverse_mapper_y = @(a, b, c) coord_cellgrid{a, b, c}(2); % Inverse map for Y
	inverse_mapper_z = @(a, b, c) coord_cellgrid{a, b, c}(3); % inverse map for Z
	
    % Apply mapper function to collectively retrieve X and Y from cellgrid
	X = permute(arrayfun(inverse_mapper_x, tA, tB, tC),[2,1,3]); 
	Y = permute(arrayfun(inverse_mapper_y, tA, tB, tC),[2,1,3]); 	
    Z = permute(arrayfun(inverse_mapper_z, tA, tB, tC),[2,1,3]); 
end % cell2mesh


function ret = subcell(cell_3d, mode)
    cell_size = size(cell_3d);
    lx = cell_size(1);
    ly = cell_size(2);
    lz = cell_size(3);
    if mode == 1 % first index fixed 
        
        % It's important to note the ordering
        subcells = cell(ly, 1);
        for i = 1 : ly
            % generate a new cell array
            new_cell = cell(lz, lx);
            
            for j = 1 : lz
                for k = 1 : lx
                    % copy elements to the new cell array
                    new_cell{j, k} = cell_3d{k, i, j};
                end
            end
            
            % append the new cell array to the list of cell arrays
            subcells{i} = new_cell;
        end
        
        ret = subcells;        

    elseif mode == 2 % second index fixed
        
        subcells = cell(lx, 1);
        for i = 1 : lx
            % generate a new cell array
            new_cell = cell(ly, lz);
            
            for j = 1 : ly
                for k = 1 : lz
                    % copy elements to the new cell array
                    new_cell{j, k} = cell_3d{i, j, k};
                end
            end
            
            % append the new cell array to the list of cell arrays
            subcells{i} = new_cell;
        end
        ret = subcells;      
        
    elseif mode == 3 % third index fixed
        subcells = cell(lz, 1);
        for i = 1 : lz
            % generate a new cell array
            new_cell = cell(lx, ly);
            
            for j = 1 : lx
                for k = 1 : ly
                    % copy elements to the new cell array
                    new_cell{j, k} = cell_3d{j, k, i};
                end
            end
            
            % append the new cell array to the list of cell arrays
            subcells{i} = new_cell;
        end
        
        ret = subcells;
    
    else
        fprintf('Invalid Mode');
        ret = 0;
    end
end % subcell

end % testing