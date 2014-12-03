%% 10. Nonlinear ODEs, IVPs, and Chebgui
% Lloyd N. Trefethen, November 2009, latest revision December 2014

%%
% Chapter 7 described Chebfun's ``chebop'' capabilities
% for solving linear ODEs (ordinary
% differential equations) by the backslash command.  We will now describe
% extensions of chebops to nonlinear problems, as well as special
% methods used for ODE initial-value problems (IVPs) as opposed to
% boundary-value problems (BVPs).  Most of the design and
% implementation of these features was done by Asgeir Birkisson in
% collaboration with Toby Driscoll.

%% 10.1  Boundary-value problems: \ and `solvebvp`
% Chebfun contains overloads `bvp4c` and `bvp5c` of MATLAB codes
% of the same names.  However, these are not our recommended methods
% for solving BVPs, and we will not discuss them here.
% Instead, we present methods based on \ and its equivalent
% command `solvebvp`.

%%
% Recall that in Chapter 7, we realized linear
% operators as chebops constructed by commands like these:
L = chebop(-1, 1);
L.op = @(x,u) 0.0001*diff(u,2) + x.*u;

%%
% We could then solve the BVP
%
% $$ 0.0001 u'' + xu = 0, \qquad u(-1) = 0, ~~ u(1) = 1 $$
%
% as follows:
L.lbc = 0; L.rbc = 1;
u = L\0;
LW = 'linewidth'; lw = 1.6;
plot(u, 'm', LW, lw), ylim([-2.5 2.5])

%%
% What's going on in such a calculation is that |L| is a prescription for
% constructing matrices of arbitrary dimensions which are Chebyshev spectral
% approximations to the operator in question.  When backslash is executed, the
% problem is solved on successively finer grids until convergence is achieved.

%%
% The object |L| we have created is a chebop:
L

%%
% (Regrettably, in output displays like this the boundary
% conditions are hard to decipher.)
% Notice that Chebfun has detected that the chebop is linear.
% Doing this automatically is not a triviality! --- see
% [Birkisson & Driscoll 2013].

%%
% The same approach also works for nonlinear problems.
% For example, 
% in Section 7.9 we hand-coded a Newton iteration to solve the nonlinear
% BVP
%
% $$ 0.001u''-u^3 = 0,\qquad  u(-1) = 1,~~ u(1) = -1. $$
%
% Chebfun can solve such a problem in an automated way
% by means of the same syntax as before.
% Switching L to N to suggest nonlinearity, let us write
N = chebop(@(x,u) 0.001*diff(u,2) - u.^3);
N.lbc = 1; N.rbc = -1;
%%
% This gives us a chebop which Chebfun recognizes as nonlinear,
N
%%
% To solve the ODE, we can write
u = N\0;
clf, plot(u)

%%
% Note that this is the same result as in Section 7.9.
% How does Chebfun solve such problems?  That is a long story, which we shall
% not tell properly here.  In brief,
% a Newton iteration (or sometimes a damped Newton iteration) 
% is carried out in ``continuous mode'', that is,
% in a space of functions rather than vectors.  Recall that to find a zero
% of a scalar function, Newton's method requires a
% derivative at each iterative step,
% and to find a zero vector of a system of equations, it
% requires a Jacobian matrix.  Here, we seek
% a zero function of a nonlinear differential operator
% equation.  For this, Newton's method requires at each step the
% continuous analogue of a Jacobian matrix, which is a 
% Frechet derivative linear operator.  This Frechet
% derivative is realized in Chebfun by a
% continuous analogue of Automatic Differentiation using
% methods described in [Birkisson & Driscoll 2011].

%%
% Here is an example with a variable coefficient, a BVP
% due to George Carrier described in Sec. 9.7 of the book
% [Bender & Orzsag 1978].  On $[-1,1]$, we seek a function
% $u$ satisfying
%
% $$ \varepsilon u'' + 2(1-x^2) u + u^2 = 1, \qquad u(-1)=u(1) = 0 . $$
%
% with $\varepsilon = 0.01$.  Here is a Chebfun formulation
% and solution.
ep = 0.01;
N = chebop(-1, 1);
N.op = @(x,u) ep*diff(u,2) + 2*(1 - x.^2).*u + u.^2;
N.bc = 'dirichlet';
u = N\1; plot(u, 'm', LW, lw)

%%
% This solution is one of several valid solutions to this problem.
% To find another, we can specify a initial guess for the Newton
% iteration that differs from Chebfun's default (a polynomial
% function constructed to satisfy the boundary conditions).
x = chebfun('x');
N.init = 2*(x.^2 - 1).*(1 - 2./(1 + 20*x.^2));
[u, info] = solvebvp(N, 1);
plot(u,'m',LW,lw)

%%
% This time, instead of using `\`, we called the method
% |solvebvp| with two output arguments, and we specified two
% output arguments.
% second output is a MATLAB struct containing data showing the norms of the
% updates during the Newton iteration, revealing a slow
% initial phase followed by eventual rapid convergence.
nrmdu = info.normDelta;
semilogy(nrmdu,'.-k',LW,lw), ylim([1e-14,1e2])

%%
% Another way to get information about the Newton iteration with nonlinear
% backlash is by setting
cheboppref.setDefaults('plotting','on')

%%
% or
cheboppref.setDefaults('display','iter')

%%
% Type |help cheboppref| for details.  Here we shall not pursue this option
% and thus return the system to its factory state:
cheboppref.setDefaults('factory')

%%
% The heading of this section refers to the command |solvebvp|. When you apply
% backslash to a nonlinear chebop, it invokes the overloaded MATLAB command
% |mldivide|; this in turn calls a command |solvebvp| to do the actual work. By
% calling |solvebvp| directly, you can control the computation in ways not
% accessible through backslash. This situation is just like the relationship in
% MATLAB between |\| and |linsolve| for solving a linear system
% of equations. See the help documentation for details.

%% 10.2  Initial-value problems: \ and `solveivp`
% For IVPs, Chebfun contains overloads `ode113`, `ode45`, and `ode15s` of
% familiar MATLAB codes.  Again, however, these are
% not our recommended methods.
% Instead, we recommend \ and its equivalent `solveivp`.

%%
% For example, suppose we want to solve the nonlinear IVP
%
% $$ u' = u^2, \qquad t\in [0,1], \quad u(0)=0.95.  $$
%
% We can set up the problem like this:
N = chebop(0, 1);
N.op = @(t,u) diff(u) - u.^2;  
N.lbc = 0.95

%%
% Chebfun has detected that this is a nonlinear
% initial value problem.  We solve it and plot the
% solution:
u = N\0;                 
plot(u,'m',LW,lw)

%%
% A major change was introduced in Version 5.1 in how
% initial-value (and final-value) problems are solved.
% Before, Chebfun used the same global spectral representations
% as for BVPs.  This usually works fine for linear problems, but for
% nonlinear ones, it is inferior to the
% method of marching by Runge-Kutta or Adams formulas.
% In Chebfun Version 5.1,
% we have accordingly switched to solving IVPs numerically by
% `ode113` (by default), converting the resulting output to
% a Chebfun representation.  (This work, a substantial job since
% higher-order equations must be reformulated as first-order
% systems, was
% carried out by Asgeir Birkisson.)  If you wish to invoke the global
% spectral method instead, you can write
u2 = solvebvp(N,0);
%%
% For this problem the method converges, giving a
% solution that is close but not the same:
norm(u-u2)

%% 
% For many nonlinear problems, however, the `solvebvp` approach
% would not converge.

%%
% Here is another example of an IVP, which happens to be linear.
% To impose two boundary
% conditions at the left, we make |N.lbc| a function
% returning an array.
N = chebop(0, 10*pi);
N.op = @(t,u) diff(u,2) + u;
N.lbc = @(u) [u-1; diff(u)];
u = N\0;
plot(u, 'm', diff(u), 'c', LW, lw)

%%
% As a third example let us solve the van der Pol equation
% for a nonlinear oscillator.  Following the example in
% the MATLAB ODE doclumentation, we consider the problem
%
% $$ u'' = 1000(1-u^2)u' - u , \qquad u(0) = 2, ~~u'(0) = 0. $$
%
% Here is a solution
N = chebop(0,20);
N.op = @(t,u) diff(u,2) - 3*(1-u.^2).*diff(u) + u;
N.lbc = @(u) [u-2; diff(u)];
%pref = chebfunpref('ivpSolver',@ode15s);
tic, u = solveIVP(N,0); toc
plot(u), shg

%%
% This van der Pol equation is not stiff, because
% coefficient 3 is not very large.  If we increase 3,
% however, the solution becomes slow:
cheboppref.setDefaults('ivpSolver',@ode113)
N.op = @(t,u) .2*diff(u,2) - (1-u.^2).*diff(u) + u;
N.lbc = @(u) [u-2; diff(u)];
u = solveIVP(N,0);
plot(u), shg

%%
% Further examples from the old Guide.

%%
% Lorenz problem

%% 10.3 Graphical user interface: Chebgui
% Chebfun includes a GUI (Graphical User Interface) called
% `chebgui` for interactive
% solution of ODE, time-dependent PDE, and eigenvalue problems.  For
% many users, this ithe single most important part of Chebfun.
% We will not describe `chebgui` here, but we encourage readers
% to give it a try.
% Be sure to note the |Demo| menu, which contains dozens of preloaded examples,
% both scalars and systems. Perhaps most important of all is the "Export to
% m-file" button, which produces a Chebfun m-file corresponding to whatever
% problem is loaded into the GUI.  This feature enables one to get going quickly
% and interactively, then switch to a Chebfun program to adjust the fine points.
% To start exploring, just type `chebgui`.

%% 10.4 References
%
% [Bender & Orszag 1978] C. M. Bender and S. A. Orszag, _Advanced
% Mathematical Methods for Scientists and Engineers_, McGraw-Hill, 1978.
%
% [Birkisson 2014] A. Birkisson, _Numerical
% Solution of Nonlinear Boundary Value Problems for
% Ordinary Differential Equations in the Continuous
% Framework_, D. Phil. thesis, University of Oxford, 2014.
%
% [Birkisson & Driscoll 2011] A. Birkisson and T. A. Driscoll,
% "Automatic Frechet differentiation for
% the numerical solution of boundary-value problems",
% _ACM Transactions on Mathematical Software_, 38 (2012), 1-26.
%
% [Birkisson & Driscoll 2013] A. Birkisson and T. A. Driscoll,
% "Automatic linearity dection", preprint, `eprints.maths.ox.ac.uk`, 2013.
