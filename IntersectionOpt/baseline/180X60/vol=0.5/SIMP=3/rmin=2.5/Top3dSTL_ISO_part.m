% A SIMPLE STL WRITER FOR TOP3D BY LIU (JAN 2015)
% -------------- USER INPUT ------------------
clear all
close all
% Input .mat file name
fin  = 'comp.mat';
% Output .stl file name
fout = 'comp.stl';
% Density cut off
cutoff = 0.5;
% Unit Length of an element
nx = 1; % mm
ny = 1; % mm
nz = 1; % mm
% For visualization
fcolor = 'r';
Alpha  = 1;
% -------------- END USER INPUT ------------------


% Load data points
load(fin, 'xPhys_thick')
buffer=xPhys_thick;
[nely, nelx, nelz] = size(buffer);

aux = zeros(nely+2, nelx+2, nelz+2);
aux(2:end-1,2:end-1,2:end-1) =buffer;

[X,Y,Z] = meshgrid(0:nx:nx*(nelx+1), 0:ny:ny*(nely+1), 0:nz:nz*(nelz+1));

[faces,verts] = isosurface(X-0.5,Z-0.5,Y-0.5,aux,cutoff);

% Facets
facets = single(verts);
facets = reshape(facets(faces',:)',3,3,[]);
% facets: (:,:,1) --> Vertices of face 1, 
% facets(:,1,1)   --> First vertice of face 1
V1 = squeeze(facets(:,2,:) - facets(:,1,:));
V2 = squeeze(facets(:,3,:) - facets(:,1,:));

% Normal vectors
normals = cross(V1, V2);
clear V1 V2
% Normal vectors normalization
normals = bsxfun(@times, normals, 1 ./ sqrt(sum(normals .* normals, 1)));

facets = cat(2, reshape(normals, 3, 1, []), facets);
clear normals

% Write STL file
%{
FORMAT:

UINT8[80] ? Header
UINT32 ? Number of triangles

foreach triangle
    REAL32[3] ? Normal vector
    REAL32[3] ? Vertex 1
    REAL32[3] ? Vertex 2
    REAL32[3] ? Vertex 3
    UINT16 ? Attribute byte count
end

%}

fid = fopen(fout, 'wb+');
fprintf(fid, '%-80s', 'MyTop3d');              % Title
fwrite(fid, size(facets, 3), 'uint32');        % Number of facets
facets = typecast(facets(:), 'uint16');        % Convert to unit16
facets = reshape(facets, 12*2, []);
facets(end+1, :) = 0;                          % Add color(0) to the end of each facet
fwrite(fid, facets, 'uint16');
fclose(fid);
fprintf('Wrote %d facets\n',size(facets, 3));

% Visualization
cla, hold on, view(30,30), rotate3d on, axis equal, axis([0 nx*nelx 0 nz*nelz 0 ny*nely]), box
set(gca, 'YDir', 'reverse', 'ZDir', 'reverse', 'ZtickLabel', flipud(get(gca,'Ztick')'));

p = patch('Faces', faces, 'Vertices', verts,...
    'FaceColor', fcolor, 'EdgeColor', 'none', 'FaceAlpha', Alpha);
camlight, lighting gouraud;
drawnow

% =========================================================================
% === This code was written by K Liu, Dept. of Mechanical Engineering   ===
% === Indiana University-Purdue University Indianapolis,                ===
% === Indiana, United States of America                                 ===
% === ----------------------------------------------------------------- ===
% === Please send your suggestions and comments to: kailiu@iupui.edu    ===
% === ----------------------------------------------------------------- ===
% === This code is an extension for Top3D program that helps users to   ===
% === generate STL file for 3D printing or post processing              ===
% === ----------------------------------------------------------------- ===
% === This code as well as Top3d program can be downloaded freely from  ===
% === the website: http://www.top3dapp.com/                             ===
% === ----------------------------------------------------------------- ===
% === Disclaimer:                                                       ===
% === The authors reserves all rights for the program.                  ===
% === The code may be distributed and used for educational purposes.    ===
% === The authors do not guarantee that the code is free from errors, and =
% === they shall not be liable in any event caused by the use of the code.=
% =========================================================================