%% 20.  Diskfun
% Heather Wilber, June 2016

%% 
% Diskfun is a new part of Chebfun designed for computing with 
% bivariate functions defined on the unit disk. 
% When working with functions on the disk, it is often convenient to
% express them in terms of polar coordinates: Given a function $f(x,y)$ 
% expressed in Cartesian coordinates, we apply the following transformation
% of variables: 
% \begin{equation}
% x = \rho\cos\theta, \qquad y=\rho\sin\theta.
% \end{equation}
% This finds $f(\theta, \rho)$, where $\theta$ is the \textit{angular}
% variable and $\rho$ is the \textit{radial} variable. 

%%
% We can construct a function in Diskfun using either coordinate system;
% the default is Cartesian.
% A function in Cartesian coordinates can be constructed as follows: 
%%
g = diskfun(@(x,y) exp(-10*( (x-.3).^2+y.^2)));



%%
% To construct the same function using polar coordinates, we include the 
% flag 'polar' in the construction command. The result using either 
% coordinate system is the same up to approximately machine precision: 

%%

f = diskfun(@(t, r) exp(-10*( (r.*cos(t)-.3).^2+(r.*sin(t)).^2)), 'polar');
plot(f)
view(3)
norm(f-g)


%%
% The object we have constructed is called a diskfun. 
%%
f


%% 
% The output describes the {\em numerical rank} of $f$, as well as the
% approximate maximum absolute value (the vertical scale).

%%
% Conveniently, we can evaluate $f$ using points written in either 
% polar or Cartesian coordinates. This setting is stored as the 'coords' 
% parameter. Here, it has been set to polar coordinates: 
%%
f.coords


%%
% To evaluate in Cartesian coordinates, we must include the 'cart' flag: 

[   f(sqrt(2)/4, sqrt(2)/4, 'cart')    f(pi/4,1/2)  ]

%%
% Alternatively, we can change the default `coords' setting:
%%
f.coords = 'cart'

f(sqrt(2)/4, sqrt(2)/4)
%%
% We can also evaluate a univariate 'slice' of $f$, either radial or 
% angular. The result is returned as a chebfun, 
% either nonperiodic or periodic (respectively).
% Here, we plot three angular slices at
% the fixed radii $\rho = 1/4, 1/3, 1/2$. 
%%
f.coords='polar';
f1 = f(:,1/4);
f2 = f(:,1/3);
f3 = f(:,1/2);
plot(f1,'b', f2,'r', f3,'k')


%%
% As with the rest of Chebfun, Diskfun is designed to perform operations at
% essentially machine precision, and using Diskfun requires no special 
% knowledge concerning the underlying discretization procedures.  
% Currently, there are more than 100 operations available.

%% Basic operations
% A suite of commands are available for computing with functions on the 
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
% We can also find global maxima and minima. Here, we plot a function
% along with its maximum value.

%%
 f = @(x,y) cos(15*((x-.2).^2+(y-.2).^2)).*exp(-((x-.2)).^2-((y-.2)).^2);
 f = diskfun(f);
 [j, k] = max2(f) 
 plot(f)
 colorbar 
 hold on
 plot3( k(1),k(2),j, 'k.', 'Markersize', 30);
 hold off
 
 %%
 % If it is preferred, the location of the maximum can be returned in polar
 % coordinates: 
 
 [jp, kp] = max2(f, 'polar');
 
 [kp(2)*cos(kp(1)) kp(2)*sin(kp(1))]
 
%%
% We can visualize a diskfun in many ways.  Here is a contour plot,
% with the zero contours displayed in black:
%% 
contour(f, 'Linewidth', 1.2)
colorbar
hold on
contour(f, [0 0], '-k', 'Linewidth', 2)
hold off
 
%%
% The roots of a function (1D contours)
% can also be found explicitly. They are stored as 
% a cell array of chebfuns. Each cell consists of two chebfuns 
% that parametrize the contour. 
%%
r = roots(f);   %note: these look pretty terrible right now. Need to figure out
                % how to fix the gap
                % it may be that we just need to pick a different function
plot(f)
hold on
for j = 1:length(r)
    rj = r{j}; 
    plot(r{j}(:,1), r{j}(:,2), '-k', 'Linewidth', 2)
end
hold off
%%
% As with chebfun, we can add, subtract, or multiply diskfuns 
% together. 

g = diskfun(@(th, r) -40*(cos(((sin(pi*r).*cos(th)...
    + sin(2*pi*r).*sin(th)))/4))+39.5, 'polar');
plot(g)
snapnow

h = g+ f;
plot(h)
snapnow

h = g - f;
plot(h)
snapnow
 
h = g.*f; 
plot(h)

%%
% Differentiation on the disk with respect to the polar variable $\rho$
% can lead to singularities, even for smooth functions. 
% For example, the function $f(\rho, \theta) = \rho^2$ is smooth 
% on the disk, but $\partial f/ \partial \rho = 2 \rho$ is not 
% smooth. For this reason, differentiation in
% diskfun is done with respect to the Cartesian coordinates, $x$ and $y$.
% TO DO: (pick a better example.)

%%
% Here we examine a  harmonic conjugate pair of functions.
   u = diskfun(@(x,y) x.*y-x+y);
   v = diskfun(@(x,y) -1/2*x.^2+1/2*y.^2-x-y);
%%
% We can check that these functions satisfy the Cauchy-Riemann equations. 
% First, we compute $\partial u / \partial y$: 
 
%%
dyu = diffy(u);  
%%
% Next, we find $\partial v / \partial x$: 
 
dxv = diffx(v); 
 
%%
% It should be true that $u_y +v_x =0$: 
 
norm(dyu+dxv)
 
%% 
% Likewise, we check that $\partial u /\partial x$ =\partial v / \partial
% y$: 
 
norm(diffx(u) -diffy(v))

 
%%
% TO DO:(maybe replace with this example instead): 
% Here, we examine some derivatives of a function involving the Bessel
% function.
f = diskfun(@(x,y) besselj(0, 5*y).*besselj(0, 5*(x-.1)).*exp(-x.^2-y.^2));
plot(f)
snapnow

plot(diffx(f))
title('derivative of f with respect to x')
snapnow

plot(diffy(f))
title('derivative of f with respect to y')
snapnow

plot(laplacian(f))
title('Scalar laplacian of f')
snapnow

%% Poisson's equation and the cylindrical harmonic functions
% We can use Diskfun to find solutions to Poisson's equation on the disk. 
% In this example, we find $u(\theta, \rho)$ in 
% \[ \Delta^2 u = f, \qquad f(\theta, 1) = 1, \]
% where $(\theta, \rho) \in [-\pi, \pi, 0, 1]$ and 
% $f = sin\left( 21 \pi \left(1 + \cos(\pi \rho)
% \right) \rho^2-2\rho^5\cos \left( 5(t-.11)\right) \right. 
%%
f = @(t,r) sin(21*pi*(1+cos(pi*r)).*(r.^2-2*r.^5.*cos(5*(t-0.11))));
rhs = diskfun(f, 'polar'); % rhs 
bc = @(t) 0*t+1;  %boundary condition
u = diskfun.poisson(f, bc, 512) % solve for u with 512 x 512 degrees of freedom

plot(rhs)
title('f'), snapnow

plot(u)
title('u'), snapnow

%%
% The accuracy of the solver can be examined by considering the
% {\em cylindrical harmonic functions}\footnote{strictly speaking, we
% consider the cylindrical harmonic functions with a fixed height},
% which are the eigenfunctions of
% Laplace's equation on the disk. These functions form an orthogonal 
% basis that is analogous to the trigonometric basis for 1D functions on
% the unit circle. Cylindrical harmonic functions with a fixed height are 
% defined by two parameters; they are of the form 
% \[V_n^\ell(\theta, \rho) = $A_n^\ell J_n(\omega_\ell \rho)e^{in\theta}\], 
% where $A_n^\ell$ is a normalization constant, $J_n$ is an $n$th order 
% Bessel function, and $\omega_\ell$ is the $\ell$th zero of $J_n$.
% They are easily constructed in Diskfun: 

%%
f = diskfun.diskharm(5, 6); % n = 5, ell = 6

plot(f)


%%
% The exact solution to $\Delta^2u = f$ with homogeneous Dirichlet boundary 
% conditions is (get correct eigenvalue C) $Cf$. 

%%
% TO DO
%C= eigenvalue;
%norm(u -C*f)



%% Vector Calculus
% Vector-valued functions and operations on the disk are performed using 
% {\em Diskfunv}. Here, we define a diskfun and find its gradient, which is
% returned as a diskfunv object.

%%
psi = @(x,y) 10*exp(-10*(x+.3).^2-10*(y+.5).^2)+10*...
    exp(-10*(x+.3).^2-10*(y-.5).^2)+ 10*(1-x.^2-y.^2)-20;

phi = @(x,y) 10*exp(-10*(x-.6).^2-10*(y).^2);

f = diskfun(psi)+diskfun(phi);
u = grad(f)


%%
% The vector-valued function $\vec{u}(x,y)$ consists of two components, 
% ordered as $x$ and $y$, respectively. Each of these is stored as 
% a diskfun. The Cartesian coordinate system is used because it is 
% common to do so in application, and this ensures that
% each component is a smooth function on the disk. The unit vectors 
% in the polar and radial directions are discontinuous at the origin of 
% the disk, and working with them can lead to singularities.
%%
plot(f)
hold on
quiver(u, 'k')
hold off

 %%
 % Once a diskfunv is created, we can apply several operations.  
 % For example, the divergence is given by
 %%
 D = div(u); 
 contour(D, 'Linewidth', 1.5)
 hold on
 quiver(u, 'k')
 hold off
 
%% 
 % Several vector calculus operations are available.
 % For example, since $\vec{u}$ is a gradient field, we expect 
 %$\Delta \cross \vec{u} =  0$
 %%
 v = curl(u);
 norm(v) % wow not actually so good. 
 

 %%
 % We can perform a variety of algebraic and
 % calculus-based operations using diskfunvs. For a complete listing of the
 % available operations, type { \tt methods diskfunv}.
 

%% Constructing a diskfun
% The underlying algorithm for constructing a diskfun adaptively selects and 
% stores a collection of 1D circular and radial ``slices" of a function $f$
% on the unit disk to create a compressed representation of $f$. 
% We compute this representation using a method of low rank approximation.

%%
% These slice are formed through the selection of pivot values sampled from
% the function, and rely on symmetry features that enforce smoothness over the pole of the
% disk [1]. The numerical rank of a diskfun corresponds to the number of slices it is composed of. We can view the slices and pivots using the plot command. 

%%
g
plot(g)
snapnow
clf
plot(g, 'k.-', 'Linewidth', 1.5, 'Markersize', 15)

%%
% There are 9 circular slices and 9 radial slices. The astersiks
% represent the pivot values. There are twice as many pivots because they
% are sampled symmetrically to ensure smoothness. 
% TO DO: explain how slices relate to numerical rank, plot column and row
% slices, give Fourier-Chebyshev series and show coeffs2...reference BMC
% and spherefun. 
%% References
% [1] A. Townsend, H. Wilber, and G. Wright, Numerical computation with 
% functions defined on the sphere (and disk), submitted, 2016.

