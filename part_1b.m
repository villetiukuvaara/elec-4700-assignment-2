%% Part 1b
% Next, the solution to Laplace's case for a 2D case with constant charge
% density was investigated. The the top and bottom boundaries (x = 0 and
% x = L) are at V0 and the other sides are at 0 V. First, the FD method was
% used with various grid sizes. This was compared to the analytic solution,
% which involves an infinite series.

clear all;
close all;

W = 2;
L = 3;
V0 = 1;

dx = 0.1; % Mesh spacing along x
dy = 0.1; % Mesh spacing along y
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
        G(index, mapCoordinate(x-1,y,nx)) = a2;
        G(index, mapCoordinate(x+1,y,nx)) = a2;
        G(index, mapCoordinate(x,y-1,nx)) = a3;
        G(index, mapCoordinate(x,y+1,nx)) = a3;
    end
end

%%
% Next, the F matrix is generated. The boundary is set to V0 along the
% sides where x = 0 and x = L.

F = zeros(nx*ny,1);

for x=1:nx
    index = mapCoordinate(x,1,nx);
    G(index,index) = 1;
    
    index = mapCoordinate(x,ny,nx);
    G(index,index) = 1;
end

for y=1:ny
    index = mapCoordinate(1,y,nx);
    G(index,index) = 1;
    
    F(index) = V0;
    
    index = mapCoordinate(nx,y,nx);
    G(index,index) = 1;
    
    F(index) = V0;
end

%%
% The corners have the boundary condition multiply defined. The
% analytical solution has V = 0 at the corners, so use this value.

F(mapCoordinate(1,1,nx)) = 0;
F(mapCoordinate(1,ny,nx)) = 0;
F(mapCoordinate(nx,1,nx)) = 0;
F(mapCoordinate(nx,ny,nx)) = 0;

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

%%
% The FD solution, above, looks good at first glance. For comparison, the
% analytic solution was calculated. It involves an infinite series, which
% was summed for various numbers of iterations. Below it is shown for 100
% iterations.
% Also, the average error after each iteration was calculated by finding
% the voltage error at each point and taking the average over all points.

analyticalSoln = zeros(ny, nx);
xx = repmat(linspace(-L/2,L/2,nx),ny,1);
yy = repmat(linspace(0,W,ny),nx,1)';
iterations = 100;
avgError = zeros(iterations,1);


for i=1:iterations
    n = 2*i - 1;
    analyticalSoln = analyticalSoln + 1./n.*cosh(n.*pi.*xx./W) ...
        ./cosh(n.*pi.*(L./2)./W).*sin(n.*pi.*yy./W);

    avgError(i) = mean(mean(abs(analyticalSoln.*4.*V0./pi - soln)));
end

analyticalSoln = analyticalSoln.*4.*V0./pi;

figure(2);
surf(linspace(0,L,nx),linspace(0,W,ny),analyticalSoln);
xlabel('x');
ylabel('y');
title(sprintf('Analytical Solution with %d iterations', iterations));

figure(3);
plot(1:i,avgError);
xlabel('Iteration');
ylabel('Average Error (V)');
title('Convergence of Analytical Solution');
grid on;

%%
% The analytical solution is seen to quickly converge within about 10
% iterations. After this point, the solution does not improve very much,
% and does not go below 0.005 V. This was for a mesh of 0.1. The simulation
% was repeated for different mesh sizes. With 0.5, the convergence was
% slower, and the error was about 0.15 V after 100 iterations. This appeals
% to intuition, since increasing the mesh size makes the second derivative
% approximation for FD poorer.
%
% The analytical solution is clearly easier to implement, as evidenced by
% the much shorter amount of code written for it. It took much more effort
% to generate to generate the G and F matrices for FD. However, finding the
% expression for the anaytical solution is a tricky task. Furthermore, this
% was for very
% simple geometry and boundary conditions, yet the analytical expression is
% not trivial. For more complex problems, the task of finding an analytic
% solution is very difficult or impossible. The FD method is more flexible
% and allows more complex geomtry or boundary conditions to be easily
% implemented.
