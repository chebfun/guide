%% CHEBFUN GUIDE 11: NONLINEAR PDES VIA PDE15S (EXPERIMENTAL)
% Nick Hale, January 2011

%%
% Chebfun has an experimental feature PDE15S for solving a class of
% time-dependent Partial Differential Equations (PDEs) of the form
%        u_t(x,t) = F(u(x,t),x,t),                                      (1)
% where F may be a nonlinear function of u involving derivatives of u
% with respect to x.

%% Example
% Here is an example solving the nonlinear advection-diffusion equation
%        u_t = .02*u_xx + u_x                                          (2)
% with periodic boundary conditions on [-1 1] for t in [0, 3] and an
% initial condition
%        u(x,0) = exp(3*sin(pi*x)).

F = @(u,t,x,diff) .02*diff(u,2) - diff(u);
tt = 0:.05:3;
u0 = chebfun(@(x) exp(3*sin(pi*x)), [-1 1]);
uu = pde15s( F , tt , u0 , 'periodic' );

%%
% The variable tt specifies the time domain [tt(1) tt(end)] on which we
% want to solve the PDE, as well as the times at which we would like the
% value of the solution returned. The variable uu is a chebfun quasimatrix
% where uu(:,k) is the solution to the PDE at tt(k). We can then visualise
% the solution on a surface plot.

surf(uu,tt,'edgealpha',0)
xlabel('x'), ylabel('t')

%%
% Perhaps a cleaner way of visualising the solution is with a 'waterfall'
% plot, where now the lines from left to right show the columns of the
% quasimatrix uu at the times tt.

waterfall(uu,tt,'simple')
xlabel('x'), ylabel('t')

%% How does PDE15S work? 
% PDE15S uses a 'method of lines' type approach, where the solution in the
% spacial dimension is realised as a Chebyshev Spectral discretisation and
% the resulting 'ODE' u_t = F(u) is advanced in time using Matlab's ODE15S
% algorithm (hence the name, 'PDE15S').
%
% Boundary conditions are obviously important, and these are dealt with by
% passing a singular mass mmatrix into ODE15s. That is, we solve the system
%         [I 0] u_t = [F(u)]
%         [0 0]       [G(u)]
% where G(u) expresses the (possibly nonlinear and/or time-dependent)
% spacial boundary conditions.

%% Syntax
% Let's return to our original example
%
%   F = @(u,t,x,diff) -(1+0.6*sin(pi*x)).*diff(u);
%   tt = 0:.05:3;
%   u0 = chebfun(@(x) exp(3*sin(pi*x)), [-1 1]);   
%   uu = pde15s( F , tt , u0 , 'periodic' );
%
% The anonymous function F clearly represents the righthand side of (2),
% but the arguments F takes are a little confusing. For technical reasons
% we must pass a variable 'diff' to F to show what function represents the
% differential operator in the equation. (The precise reason for this is
% that when the function F is evaluated within PDE15S the variable u is a
% vector, and so diff(u) would simply compute a standard finite difference,
% which is not what is wanted here). Other than this, F simply takes the
% inputs u, t, and x as expected. F should be vectorised.

%%
% As mentioned above, the next input, tt, to PDE15S defines the solution
% times at which the output uu should be given. This is also used
% internally within PDE15S to determine the size of the 'time chunks' over
% which to solve at various spatial disrectisation sizes, and selecting
% this appropriately can have a siginificant effect of the efficiency of
% the routine. FOr future version we are hoping to automate this process.

%%
% The next input, u0, provides the initial condition at t = 0 and also the
% domain on which the PDE is defined. 

%%
% The final input contains the boundary conditions. These can take the form
% of some special key words, such as 'dirichlet', 'neumann', or 'periodic'
% as described in Chapter 7.4, or as a function handle with the same format
% as F. For example, here is the same PDE as solved above, but here with an
% initial condition u(x,0) = 0, a Neuman condition u_x(-1,t) at x=-1 for
% all t, and a nonlinear time-dependent condition
%       u(1,t) = (1+u^2)*sin(pi*u)-(1-exp(-t)).*cos(10*t);

u0 = 0*x;
bc.left = 'neumann';
bc.right = @(u,t,x,diff) (1+u.^2).*sin(pi*u)-(1-exp(-t)).*cos(10*t);
uu = pde15s( F , tt , u0 , bc );

%% Integro-Differential Equations
% PDE15S can also handle partial integro-differential equations, i.e.,
% where the rhs function F also involves definite and indefinite integrals
% of F. Below is an example, ('exported' from CHEBGUI -- see Chapter 12)
% solving the equation:
%                           /x       /1
%        u_t = 0.02*u_xx + | u dx * | u dx ,
%                          /-1      /-1
%
% on the time interval t = [0 4] with u(-1,t) = u(1,t) = 0, and 
%        u(x,0) = (1-x^2)*exp(-30*(x+.5)^2).

%%
% Notice how the definite and indefinite integral operations, 'sum' and
% 'cumsum' are passed in the arguments of pdefun in the same way as 'diff'.

% Create a domain and the linear function on it.
d = [-1,1];
x = chebfun('x',d);

% Construct a discretisation of the time domain to solve on.
t = 0:.1:4;

% Make the rhs of the PDE.
pdefun = @(u,t,x,diff,sum,cumsum) .02*diff(u,2)+cumsum(u).*sum(u);

% Assign boundary conditions.
bc.left = 'dirichlet';
bc.right = 'dirichlet';

% Create a chebfun of the initial condition.
sol0 = (1-x.^2).*exp(-30.*(x+.5).^2);

% Solve the problem using pde15s.
[t sol] = pde15s(pdefun,t,sol0,bc);

% Create plot of the solution.
waterfall(sol,t,'simple','linewidth',2)

%% Systems
% PDE15S also supports systems of equations. here is an example of a
% chemical reaction between species u and v which diffuse along a tube and
% interact to produce chemical w. The Neumann condition is automatically
% applied to each of the variables at both ends of the domain.

x = chebfun('x',[-1 1]);  
u0 = [ 1-erf(10*(x+0.7)) , 1 + erf(10*(x-0.7)) , 0 ];
f = @(u,v,w,diff)  [ .1*diff(u,2) - 100*u.*v , ...
                    .2*diff(v,2) - 100*u.*v , ...
                    .001*diff(w,2) + 2*100*u.*v ];
bc = 'neumann';     
uu = pde15s(f,0:.1:3,u0,bc);

%%
% Here we plot the growth of the new chemical v in time and space.
mesh(uu{3})
xlabel('x'), ylabel('t');

%% PDESET
% Various options can be passsed with the PDESET functions, which performs
% a similar role to ODESET for Matlab's ODE tools. See "help pdeset".

%%
% MORE HERE.

%% Passing the PDE as a chebop
% PDE15S can also be called with a chebop as the input.
%       u = pde15s( N , tt )
% where N.op represents the rhs function F, N.guess holds the initial
% condition. See "help chebop/pde15s".

%%
% MORE HERE.
