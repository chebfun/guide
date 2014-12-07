%% 14. Chebfun2: Rootfinding and Optimisation
% Alex Townsend, March 2013, latest revision December 2014

%% 14.1 Zero contours 
% Chebfun2 comes with the capability to compute the zero contours of a
% function of two variables. For example, here we compute a representation
% of Trott's curve, an example from algebraic geometry [Trott 1997].
x = chebfun2(@(x,y) x); y = chebfun2(@(x,y) y); 
trott = 144*(x.^4+y.^4)-225*(x.^2+y.^2) + 350*x.^2.*y.^2+81;
r = roots(trott);
LW = 'linewidth'; MS = 'markersize';
plot(r,LW,1.6), axis([-1 1 -1 1]), axis square

%%
% The zero curves are represented as complex-valued chebfuns
% (see Section 13.4).   For example, 
r(:,1)

%%
% The zero contours of a function are computed by Chebfun2 to plotting
% accuracy: they are typically not accurate to machine precision. 

%% 14.2 |roots|
% Chebfun2 also comes with the capability of finding zeros of bivariate
% systems, i.e., the solutions to $f(x,y) = g(x,y) = 0$. 
% If the |roots| command is supplied with one chebfun2, it computes the zero
% contours of that function, as in the last
% section.  However, if it is supplied
% with two chebfun2 objects, as in |roots(f,g)|, then it 
% computes the roots of the bivariate system. Generically, these are
% isolated points.

%%
% What points on Trott's curve intersect the circle of radius $0.9$?
g = chebfun2(@(x,y) x.^2 + y.^2 - .9^2);
r = roots(trott,g); 
plot(roots(trott),'b',LW,1.6), hold on
plot(roots(g),'r',LW,1.6)
plot(r(:,1),r(:,2),'.k',LW,1.6,MS,20)
axis([-1 1 -1 1]), axis square, hold off

%%
% The solutions to bivariate polynomial systems and intersections of curves are
% typically computed to full machine precision.

%% 14.3 Intersections of curves
% The problem of determining the intersections
% of real parameterised complex curves can be expressed as a
% bivariate rootfinding problem.  For instance, here are the intersections
% between the 'splat' curve [Guettel 2010] and a 'figure-of-eight'
% curve. 
t = chebfun('t',[0,2*pi]);
sp = exp(1i*t) + (1+1i)*sin(6*t).^2;     % splat curve
figof8 = cos(t) + 1i*sin(2*t);           % figure of eight curve
plot(sp,LW,1.6), hold on
plot(figof8,'r',LW,1.6), axis equal

d = [0 2*pi 0 2*pi];
f = chebfun2(@(s,t) sp(t)-figof8(s),d);  % rootfinding
r = roots(real(f),imag(f));              % calculate intersections
spr = sp(r(:,2)); 
plot(real(spr),imag(spr),'.k',MS,20), ylim([-1.1 2.1])
hold off

%%
% Chebfun2 rootfinding is based on an algorithm described in
% [Nakatsukasa, Noferini & Townsend 2014].

%% 14.4 Global optimisation: |max2|, |min2|, and |minandmax2|
% Chebfun2 also provides functionality for global optimisation. Here is
% an example, where we plot the minimum and maximum as black dots.
f = chebfun2(@(x,y) sin(30*x.*y) + sin(10*y.*x.^2) + exp(-x.^2-(y-.8).^2));
[mn mnloc] = min2(f); 
[mx mxloc] = max2(f); 
plot(f), hold on 
plot3(mnloc(1),mnloc(2),mn,'.r',MS,40)
plot3(mxloc(1),mxloc(2),mx,'.r',MS,30) 
zlim([-6 6]), colormap('bone'), hold off

%% 
% If both the global maximum and minimum are required, it is roughly twice
% as fast to compute them at the same time by using the |minandmax2| command.
% For instance, 
tic; [mn mnloc] = min2(f);  [mx mxloc] = max2(f); t=toc;
fprintf('min2 and max2 separately = %5.3fs\n',t)
tic; [Y X] = minandmax2(f); t=toc;
fprintf('minandmax2 command = %5.3fs\n',t)

%%
% Here is a complicated function from the 2002 SIAM 100-Dollar,
% 100-Digit Challenge.  Chebfun2 computes its global minimum
% in a few seconds:
tic
f = cheb.gallery2('challenge');
[minval,minpos] = min2(f);
minval
toc

%%
% The result closely matches the correct solution, computed
% to 10,000 digits Bornemann et al. [2004]:
exact = -3.306868647475237280076113

%%
% Here is a contour plot of this wiggly function, with the
% minimum circled in black: 
colormap('default'), contour(f), hold on
plot(minpos(1),minpos(2),'ok',MS,20,LW,3), hold off

%% 14.5 Critical points
% The critical points of a smooth function of two variables can be located 
% by finding the zeros of $\partial f/ \partial y = \partial f / \partial x = 0$.   
% This is a rootfinding problem.  For example,
f = chebfun2(@(x,y) (x.^2-y.^3+1/8).*sin(10*x.*y));
r = roots(gradient(f));                       % critical points
plot(roots(diff(f,1,2)),'b',LW,1.2), hold on  % zero contours of f_x
plot(roots(diff(f)),'r',LW,1.2)               % zero contours of f_y
plot(r(:,1),r(:,2),'k.',MS,24)                % extrema
axis([-1,1,-1,1]), axis square

%% 
% There is a new command here called |gradient| that computes the gradient
% vector and represents it as a chebfun2v object.
% The |roots| command then solves for the isolated roots of the bivariate 
% polynomial system represented in the chebfun2v representing the gradient.  
% For more information about |gradient|, see Chapter 15. 
 
%% 14.6 Infinity norm
% The $\infty$-norm of a function is the maximum absolute value in its 
% domain. It can be computed by passing the 
% argument |inf| to the |norm| command.  
f = chebfun2(@(x,y) sin(30*x.*y));
norm(f,inf)

%% 14.7 References
%
% [Bornemann et al. 2004] F. Bornemann, D. Laurie, S. Wagon and
% J. Waldvogel, _The SIAM 100-Digit Challenge: A Study in
% High-Accuracy Nuemrical Computing_, SIAM 2004.
%
% [Guettel 2010] S. Guettel, 
% "Area and centroid of a 2D region", 
% `www.chebfun.org/examples/geom/Area.html`.
% 
% [Nakatsukasa, Noferini & Townsend 2014] Y. Nakatsukasa, V. Noferini
% and A. Townsend, "Computing the common zeros of two bivariate functions
% via Bezout resultants", _Numerische Mathematik_, to appear.
%%
% [Trott 2007] M. Trott, "Applying GroebnerBasis to three problems in
% geometry", _Mathematica in Education and Research_, 6 (1997), 15-28.
