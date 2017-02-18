%% Part 1a
% First, the solution to Laplace's case for the simple 2D case, which is
% essentailly a 1D case, was found. Here, the charge density is uniform.

clear all;
close all;

W = 2;
L = 3;
V0 = 1;

dx = 0.25; % Mesh spacing along x
dy = 0.25; % Mesh spacing along y
nx = L/dx; % Number of points along x
ny = W/dy; % Number of points along y

%%
% The finite differences method can be implemented using a matrix
% calculation, where $GV=F$. $V$ is the voltages at the discrete points,
% $F$ is the "forcing" matrix that is used to set the boundary conditions,
% and G defines how the voltages are each point are related to the other
% voltages. The compute G, the the equation derived in class was used:
% 
% $$\frac{V_{x-1,y}-2V_{x,y}+V_{x+1,y}}{(\Delta x)^2} + \frac{V_{x,y-1}-2V_{x,y}+V_{x,y+1}}{(\Delta y)^2}=0$$
% 
% This equation is Laplace's equation where the second derivatives have
% been converted to discrete approximations. When generting $G$, the
% coefficients of $V_{x-1,y}$, $V_{x,y}$, etc., are needed. To avoid
% calculating these values at every point, they are calculated here just
% once:

a1 = -2*(1/dx^2 + 1/dy^2);
a2 = 1/(dx^2);
a3 = 1/(dy^2);

%%
% Here, the G matrix is generated. The mapCoordinate function takes a
% discrete (x,y) coordinate and converts it to an index in the V array.
% This mapping is performed such that the V array consists of successive
% rows (points with the same y-value).

G = zeros(nx*ny,nx*ny);

for x=2:(nx-1)
    for y=2:(ny-1)
        index = mapCoordinate(x,y,nx);
        G(index,index) = a1;
        G(index,mapCoordinate(x-1,y,nx)) = a2;
        G(index,mapCoordinate(x+1,y,nx)) = a2;
        G(index,mapCoordinate(x,y-1,nx)) = a3;
        G(index,mapCoordinate(x,y+1,nx)) = a3;
    end
end

%%
% Next, the F matrix is generated. The boundary is set to V0 along
% x = 0 and 0 along x = L.

F = zeros(nx*ny,1);

for y=1:ny
    index = mapCoordinate(1,y,nx);
    G(index,index) = 1;
    
    F(index) = V0;
    
    index = mapCoordinate(nx,y,nx);
    G(index,index) = 1;
end

%%
% For y = 0 and y = L, the derivative is set to 0. This is done by setting
% the boundary to the same voltage as the adjecent node.

for x=2:(nx-1)
    index = mapCoordinate(x,1,nx);
    G(index,index) = 1;
    G(index,mapCoordinate(x,2,nx)) = -1;
    
    index = mapCoordinate(x,ny,nx);
    G(index,index) = 1;
    G(index,mapCoordinate(x,ny-1,nx)) = -1;
end

%%
% After setting up the matrices, getting the solution is trivial. To
% convert back to a matrix where the voltages can be accessed with (x,y)
% coordinates, the reshape function is used.

soln = G\F;
soln = reshape(soln,[],ny)';

figure(1);
surf(linspace(0,L,nx),linspace(0,W,ny),soln);
xlabel('x');
ylabel('y');
title(sprintf('Finite Differences Solution with Grid Spacing of %.2f', dx));
set(gca, 'View', [45 45]);

%%
% The FD solution looks correct. It is linear along x, with V = V0 = 1 at
% x = 0 and V = 0 at x = L.