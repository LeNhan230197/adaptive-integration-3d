% testing.m
% Script for unit testing

function [] = testing()

[X, Y, Z] = coord_meshgrid(0, 1, 0, 1, 0, 1, [2, 2, 2]);
display(X);
display(Y);
display(Z);

function [X, Y, Z] = coord_meshgrid(xmin, xmax, ymin, ymax, zmin, zmax, mesh_size)

    % Update in V0.1: Matlab doesn't support underscore _, it is replaced
    % by t to denote temporary
	tx = linspace(xmin, xmax, mesh_size(1));
	ty = linspace(ymin, ymax, mesh_size(2));
    tz = linspace(zmin, zmax, mesh_size(3));
	[X, Y, Z] = meshgrid(tx, ty, tz);
end % coord_meshgrid

end % testing