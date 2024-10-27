%................................................................
% MATLAB codes for Finite Element Analysis
% problem5.m
% A.J.M. Ferreira, N. Fantuzzi 2019

%%
% clear memory
clear

% E; modulus of elasticity
% A: area of cross section
E = 210*10^9; A=1.5; EA=E*A;

a=readmatrix("conn.csv")
b=readmatrix("nodalcoords.csv")
% generation of coordinates and connectivities
elementNodes =a;
nodeCoordinates = b;
numberElements = size(elementNodes, 1);
numberNodes = size(nodeCoordinates, 1);
xx = nodeCoordinates(:, 1);
yy = nodeCoordinates(:, 2);

% for structure:
%   displacements: displacement vector
%   force : force vector
%   stiffness: stiffness matrix
GDof = 2 * numberNodes;
U = zeros(GDof, 1);
force = zeros(GDof, 1);

% applied load at node 2 
force =readmatrix("Force.csv");

% computation of the system stiffness matrix
[stiffness] = formStiffness2Dtruss(GDof, numberElements, ...
    elementNodes, numberNodes, nodeCoordinates, xx, yy, EA);

% boundary conditions and solution
prescribedDof = [1 2 42]';

% solution
displacements = solution(GDof, prescribedDof, stiffness, force);
us = 1:2:2 * numberNodes - 1;
vs = 2:2:2 * numberNodes;

% Drawing displacements
figure
hold on

% Original coordinates for plotting
originalCoords = nodeCoordinates;

% Deformation scaling
dispNorm = max(sqrt(displacements(us).^2 + displacements(vs).^2));
scaleFact = 10e8 * dispNorm; % Scale factor for visualization

% New coordinates after applying displacements
deformedCoords = originalCoords + scaleFact * [displacements(us), displacements(vs)];

% Plot original mesh
drawingMesh(originalCoords, elementNodes, 'L2', 'k.--'); % Original mesh in dashed line

% Plot deformed mesh
drawingMesh(deformedCoords, elementNodes, 'L2', 'k.-'); % Deformed mesh in solid line

axis equal
set(gca, 'fontsize', 18)
xlabel('X Coordinate')
ylabel('Y Coordinate')
title('Deformed and Original Mesh')
legend('Original Mesh', 'Deformed Mesh')
grid on

% output displacements/reactions
outputDisplacementsReactions(displacements, stiffness, ...
    GDof, prescribedDof)

% stresses at elements
stresses2Dtruss(numberElements, elementNodes, ...
    xx, yy, displacements, E)

% The drawingMesh function remains unchanged
