%% Chapter 17: Linear Partial Differential Equations in Chebfun2
% Alex Townsend, October 2014. 

%% Introduction 
% Chebfun has new capabilities for solving partial differential equations
% defined on bounded rectangular domains. Now one can solve many 2D linear 
% PDEs to a high accuracy by a single backslash command. The syntax for 
% PDEs closely follows the syntax for solving ODEs in Chebfun 
% (see guide chapter 7). 

%%
% Disclaimer: Though we have made good progress on solving PDEs in 
% Chebfun, we still have an enormous amount to learn.  Chebop2 prints a 
% warning on its first execution indicating the experimental nature of this 
% side of Chebfun. 

%% About chebop2 objects 
% Chebop2 objects represent linear partial differential operators acting on 
% bivariate functions (chebfun2 objects). A chebop2 can be assigned 
% boundary conditions and the associated PDE solved together with a 
% right-hand side. The solution is returned as a chebfun2 
% (see guide chapter 11).
%
% The underlying mathematics is based on a global spectral method. 
% In particular, it is an extension of the ultraspherical spectral method 
% [1] to two dimensions and separable representations of differential operators, 
% see [2].  

%% Chebop2 syntax
% Like a chebop, a chebop2 has a domain, an operator, and sometimes
% boundary conditions. For example, here is the Laplace operator 
% on the unit square $[-1,1]x[-1,1]$: 

L = chebop2( @(u) laplacian( u ) ); 

%%
% This chebop2 can be applied to a chebfun2 defined on $[-1,1]x[-1,1]$. 
% For example, we can verify that a harmonic function is in the kernel of 
% the Laplacian: 

f = chebfun2( @(x,y) exp(x).*sin(y) );   % harmonic function of 2 variables
norm( L * f )                            % equivalent to norm(laplacian(f))

%%
% A chebop2 can also be assigned boundary conditions. Homogeneous Dirichlet
% conditions are particularly convenient and can be assigned as follows: 

L.bc = 0; 

%% 
% More generally, boundary conditions can be supplied to the four sides of the 
% square. For example,

L.lbc = 1;                              % left boundary conditions
L.rbc = @(y,u) u - cos(pi*y) - 2;       % right boundary conditions
L.ubc = @(x,u) u - exp(-20*x.^2) + 1;   % top boundary conditions
L.dbc = 1;                              % bottom boundary conditions

%% 
% Neumann and Robin boundary conditions are also possible, see example 3 
% below.

%% 
% While it is possible to supply Dirichlet data that is 
% discontinuous at the corners of the square, it is rarely appropriate in 
% this setting since global spectral methods do require the solution 
% to be continuous --- the smoother the better.  

%% 
% When L is typed into the command line without a semicolon, a summary of 
% the chebop2 is printed: 

L 

%%
% The Chebop2 backslash command solves the corresponding PDE: 

u = L \ 0;           
contour( u ), hold on 
set(gca, FS, 16)

%% 
% The solution is a chebfun2 and its properties such as its maximum, 
% integral, and roots can be investigated via existing Chebfun2 code. Here, 
% we can verify the maximum principle by showing the maximum occurs on the
% boundary: 

[mx, loc] = max2( u ); 
plot(loc(1), loc(2), '.r', MS, 30)
title( sprintf('Maximum = %1.3f', mx), FS, 16), hold off

%% Examples  
% To show the potential and typically use of Chebop2, we give three
% examples. 

%% Example 1
% Consider Laplace's equation on a strip $[-10,10]\times[-1,1]$ with 
% homogeneous Dirichlet conditions and a forcing term $\cos(x) + \sin(xy)$. 

strip = [-10 10 -1 1]; 
L = chebop2( @(u) laplacian( u ), strip ); 
L.bc = 0; 
rhs = chebfun2(@(x,y) cos(x) + sin(x.*y), strip); 
u = L \ rhs; 
plot( u ), view(0,90), axis equal, set(gca, FS, 16)
axis(strip)

%% Example 2
% Consider the Helmholtz-like equation with a linearly varying coefficient 
% given by $\nabla^2u - 1000yu$ with constant Dirichlet boundary conditions of $1$.

N = chebop2( @(x,y,u) laplacian(u) - 1000*y.*u, [-1 1 -3 0]);
N.bc = 1; 
u = N \ 0; 
plot( u ), view(0,90), axis equal, set(gca, FS, 16)
text(1.15,-1.5,'gravity',FS,16), annotation('arrow',[.8 .8],[.48 .3])
axis([-1 1 -3 0])

%% Example 3  
% Consider the Klein-Gordon equation $u_{tt} = u_{xx} - 5u$ modelling 
% a string with dispersion, where the left endpoint is held fixed 
% and the right endpoint is attached to a freely moving loop around a 
% vertical wire. 

N = chebop2(@(u) diff(u,2,1) - diff(u,2,2) + 5*u, [-1 1 0 5]); 
N.lbc = 0; N.rbc = @(t,u) diff(u); 
N.dbc = @(x,u) [u-exp(-50*(x-.2).^2) ; diff(u)];
u = N \ 0; 
plot( u ),  view(0,90), axis equal, set(gca,FS,16)
xlabel('x',FS,16), ylabel('t',FS,16), axis([-1 1 0 5])

%% 
% The solution was displayed and solved in an unfamiliar way. First, 
% the solution was calculated for all "time" and visualized 
% as a surface (for $x$ in $[-1,1]$ and $t$ in $[0 5]$). Second, the solution was 
% computed by a spectral method, rather than a time-stepping scheme, and 
% hence time treated as a spatial variable.  
% In some settings, it can be more efficient to use a time-stepping scheme.
% See pde15s in Chebfun for a PDE solver based on time-stepping. 

%% Elliptic, Parabolic, and Hyperbolic partial differential equations 
% The Chebop2 backslash command uses a dense solver, as opposed to an 
% iterative one. For this reason, the classification of the type of PDE is 
% immaterial. For instance, the underlying methodologies for the wave, 
% heat, and Laplace equation are practically identical. What matters 
% is the smoothness of the resulting solution and the splitting rank of the 
% associated differential operator (see below). 

%% The splitting rank of a partial differential operator 
% Chebfun2 approximates bivariate functions by low rank approximants, i.e., 
% outer products of functions of one variable. In an analogous way, 
% Chebop2 approximates partial differential operators by sums of 
% tensor products of ordinary differential operators. The number of terms 
% needed is called the splitting rank, for example, 

rank( L )           % splitting rank

%% 
% It is this "low-rank" (separable) structure of operators associated
% with many partial differential equations and the special structure of 
% ultraspherical spectral discretizations that makes chebop2 more efficient
% than standard 2D spectral methods.  For example, chebop2 can resolve 
% solutions that require over a million degrees of freedom.   

%% 
% PDEs associated to differential operators with a splitting rank of 1 or 
% 2 are particularly efficient to solve in Chebop2. Fortunately, and quite 
% surprisingly, many well-known PDEs
% involve operators of splitting rank 2. These include the Laplacian, Helmholtz 
% operator, and the operators associated to the heat, wave, transport, 
% and Black-Scholes equation. PDEs with operators of higher splitting rank
% can still be solved effectively in Chebop2 with the same syntax. 

%% Future work 
% This is a rapidly developing part of Chebfun and we expect future 
% software and algorithmic developments in the near future. 

%% References 
% 
% [1] S. Olver and A. Townsend, A fast and well-conditioned spectral
% method, SIAM Review, 55 (2013), pp. 462-489.
% 
% [2] A. Townsend and S. Olver, The automatic solution of partial differential 
% equations using a global spectral method, submitted, 2014.  