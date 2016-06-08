%% Introduction to Diskfun
% Heather Wilber

%% Constructing a Diskfun
% Diskfun is a new part of Chebfun designed for computing with 
% bivariate functions defined on the unit disk. 
% When working with functions on the disk, it is often convenient to
% express them in terms of polar coordinates: Given a function $f(x,y)$ 
% expressed in Cartesian coordinates, we apply the following transformation
% of variables: 
% \begin{equation}
% x = \rho\cos\theta, \qquad y=\rho\sin\theta.
% This finds $f(\theta, \rho)$, where $\theta$ is the \textit{angular}
% variable and $\rho$ is the \textit{radial} variable. 

%%
% We can construct a function in Diskfun using either coordinate system.
% A function in Cartesian coordinates is constructed as follows: 
%%
g = diskfun(@(x,y) exp(-10*( (x-.3).^2+y.^2)));
norm(f-g)


%%
% To construct the same function using polar coordinates, we include the 
% flag 'polar' in the construction command. The result using either 
% coordinate system is the same up to approximately machine precision: 

%%

f = diskfun(@(t, r) exp(-10*( (r.*cos(t)-.3).^2+(r.*sin(t)).^2)), 'polar');
plot(f)


%%
% The object we have constructed is called a diskfun. 
%%
f


%% 
% The output describes the \textit{numerical rank} of $f$, as well as the
% approximate maximum value (the vertical scale).

%%
% Conveniently, we can evaluate $f$ using points written in either 
% polar or Cartesian coordinates. 

[  f(sqrt(2)/4, sqrt(2)/4) f( pi/4,1/2, 'polar') ]


%%
% We can also evaluate a univariate 'slice' of $f$. The result is returned 
% as a chebfun. Here, we plot three univariate
% functions in $theta$ at the fixed values $\rho = 1/4, 1/3, 1/2$. Instead
% of including the flag 'polar' each time we evaluate the function, we
% change a global parameter. (to do: add this ability into diskfun) 
%%
f1 = f(:,1/4);
f2 = f(:,1/3);
f3 = f(:,1/2);
plot(f1,'b', f2,'r', f3,'k')

%%
% As with the rest of Chebfun, Diskfun is designed to perform operations at
% essentially machine precision, and using Diskfun requires no special 
% knowledge concerning the underlying discretization procedures.  
% Currently, there are around 100 operations available. (check)
%

%% Basic operations
%  A suite of commands are available for computing with functions on the 
% disk. For example, the integral of the function  
% $g(x,y) = -x^2 -3xy-(y-1)^2$ over the unit disk is 
% found as follows: 

%%
f = diskfun(@(x,y) -x.^2 - 3*x.*y-(y-1).^2)
sum2(f)

%% 
% The exact answer is $-3\pi/2$: 
%%
-3*pi/2

% Since the area of the unit disk is $\pi$, the mean of $f$ should be
% $-3/2$: 
%%
mean2(f)
%% 
% We can also integrate along a contour. 
%%

f = diskfun(@(x,y) sin(3*x.*y)); 
plot(f)
%%
% The contour integral along the unit circle should be zero. The contour integral on 
% the unit circle is particularly important on the disk and is 
% computed using the 'unitcircle' flag . (to do: add 'unitcircle')
%%
integral(f, 'unitcircle')

%%
% We can also find global maximums and minimums. Here, we plot a function
% along with its maximum value: (TO DO:change output to 'cart'):

%%
 %f = @(t,r) sin(10*(r.^2)+r.*sin(t))+r.*cos(t);
 f = @(x,y) peaks(2*(x),2*y);
 f = diskfun(f);
 [j, k] = max2(f); 
 plot(f)
 view(3)
 hold on
 plot3( k(2).*cos(k(1)),k(2).*sin(k(1)),j, 'k.', 'Markersize', 20);
 hold off
 
 %%
 % We can visualize a diskfun in many ways.  Here,
 % plot a range of contours, with the zero contours in black:
%% 
contour(f, 30)
hold on
contour(f, [0 0], '-k', 'Linewidth', 2)
hold off
 
%%
% The roots of a function can also be found explicitly. They are stored as 
% an array of chebfuns. Each cell in the array consists of two chebfuns 
% that parametrize the zero contour. 
%%
r=roots(f);
plot(f)
hold on
for j=1:length(r)
    rj=r{j}; 
    plot(rj(:,1), rj(:,2), '-k', 'Linewidth', 2)
end
hold off
 
 %%
 % Differentiation on the disk with respect to the polar variable $\rho$
 % can lead to singularities, even for smooth functions. 
 % For example, the function $f(\rho, \theta) = \rho^2$ is a smooth 
 % function on the disk, but $\partial f/ \partial \rho = 2 \rho$ is not 
 % smooth. For this reason, differentiation in
 % diskfun is done with respect to the Cartesian coordinates, $x$ and $y$.
 % TO DO: (pick a better example!)
 %%
 % Here we examine a set of harmonic conjugate pairs.
    u = diskfun(@(x,y) x.*y-x+y);
    v = diskfun(@(x,y) -1/2*x.^2+1/2*y.^2-x-y);
 %%
 % We can check that these functions satisfy the Cauchy-Riemann equations. 
 % First, we compute $\partial u / \partial y$: 
 
 %%
 dyu = diff(u, 2, 1); 
 %%
 % Next, we find $\partial v / \partial x$: 
 
 dxv = diff(v, 1, 1); 
 
 %%
 % It should be true that $u_y +v_x =0$: 
 
 norm(dyu+dxv)
 
 %% 
 % Likewise, we check that $\partial u /\partial x$ =\partial v / \partial
 % y$: 
 
 norm(diff(u, 1, 1) -diff(v, 2, 1))

 
 %%
 % TO DO:(maybe replace with this example instead): create a check for this using recurrence relationships
 
 f = diskfun(@(x,y) besselj(0, 5*y).*besselj(0, 20*x).*exp(-x.^2-y.^2))
 plot(f)
 %%
 %
 plot(diff(f, 1, 5))
 
 
 %% Vector Calculus
 % To compute with vector fields defined on the disk, diskfunv objects are
 % used. Each diskfunv object consists of two components, and each
 % component is a diskfun object. We can construct a diskfunv in several
 % ways. 
 
 %%
 % Here, we create a diskfunv using a function declaration for each
 % component: 
 
 
 %%
 % Alternatively, we can create a diskfunv using a diskfun for each
 % component. 
 
 %% 
 % Diskfun 
 
 %%
 

%% Algebra on Diskfuns
% We can add, subtract, or multiply two diskfuns together. 

f = diskfun(@(th, r) -10*cos(((sin(pi*r).*cos(th) + 10*sin(2*pi*r).*sin(th)))/4), 'polar')
g = diskfun(@(x,y) exp(-10*( (x-.3).^2+y.^2)));
subplot(2,2,1)
plot(g)
axis square
view(2)

subplot(2,2,2)
h = g+ f;
plot(h)
axis square
view(2)

subplot(2, 2, 3)
h = g - f;
plot(h)
axis square
view(2)
 
subplot(2, 2, 4)
h = g.*f; 
plot(h)
axis square 
view(2)


%% What is a diskfun?
% The underlying algorithm for constructing diskfuns adaptively selects and stores a
% collection of 1D circular and radial ``slices" of a function $f$ on the unit disk to create a representation of f$ that is
% compressed, low rank approximation to $f$. The idea used to construct
% this low rank approximation is similar to what is used in
% These slice are formed through the selection of pivot values sampled from
% the function, and rely on symmetry features that enforce smoothness over the pole of the
% disk [1]. The numerical rank of a diskfun corresponds to the number of slices it is composed of. We can view the slices and pivots using the plot command. 

%%
g =diskfun(@(th, r) -cos(((sin(pi*r).*cos(th) + sin(2*pi*r).*sin(th)))/4), 'polar');
plot(g)
hold on
plot(g, 'k*-', 'Linewidth', 1.5)

%%
% There are 9 circular slices, and 9 radial slices. The asteriks
% represent the pivot values. There are twice as many pivots because they
% are sampled symmetrically to ensure smoothness. 

%% References
% [1] A. Townsend, H. Wilber, and G. Wright, Numerical computations with 
% functions defined on the spherefun (and disk), in preparation, 2015. 

