%% 20. Diskfun
% Heather Wilber, August 2016

%%
LW = 'Linewidth';
MS = 'Markersize';
axis off

%% Introduction
% Diskfun is a new part of Chebfun designed for computing with
% 2D scalar and vector-valued functions defined on the unit disk. Conceptually,
% it is an extension of Chebfun2 to the polar setting, designed to accurately
% and efficiently perform over 100 operations with functions on the disk,
% including differentiation, integration, vector calculus, and rootfinding, among
% many other things. Diskfun was developed in tandem with Spherefun, and
% the two are algorithmically closely related [2,3].

%%
% To get started, we simply call the Diskfun constructor. In this example,
% we consider a Gaussian function.

g = diskfun(@(x,y) exp(-10*((x-.3).^2+y.^2)));
plot(g)
view(3)

%%
% When working with functions on the disk, it is sometimes convenient to
% express them in terms of polar coordinates: Given a function $f(x,y)$
% expressed in Cartesian coordinates, we apply the following transformation
% of variables:
% \begin{equation}
% x = \rho\cos\theta, \qquad y=\rho\sin\theta.
% \end{equation}
% This finds $f(\theta, \rho)$, where $\theta$ is the {\em angular}
% variable and $\rho$ is the {\em radial} variable.

%%
% To construct $g$ using polar coordinates, we include the
% flag 'polar'. The result using either
% coordinate system is the same up to essentially machine precision:

f = diskfun(@(t, r) exp(-10*((r.*cos(t)-.3).^2+(r.*sin(t)).^2)), 'polar');
norm(f-g)

%%
% The object we have constructed is called a diskfun, with a lower case ``D".
% We can find out more about a diskfun by printing it to the command line.
f

%%
% The output describes the {\em numerical rank} of $f$ (discussed below),
% as well an approximate maximum absolute value of $f$ (the vertical scale).

%%
% To evaluate a diskfun, we can evaluate using either polar or Cartesian
% coordinates.  (To evaluate in polar coordinates, we needs to include
% the 'polar' flag.)

[   f(sqrt(2)/4, sqrt(2)/4)    f(pi/4,1/2, 'polar')  ]

%%
% We can also evaluate a univariate ``slice" of $f$, either radial or
% angular. The result is returned as a chebfun, either a nonperiodic or
% periodic function, respectively. Here, we plot three angular slices at
% the fixed radii $\rho = 1/4$, $1/3$, and $1/2$.

c = f( : , [1/4 1/3 1/2] , 'polar');
plot(c(:,1), 'r', c(:,2), 'k', c(:,3), 'b')
title( 'Three angular slices of a diskfun' )

%%
% Where ever possible, we interpret commands with respect to the function
% in Cartesian coordinates. So, for example, the command {\tt diag} returns
% the radial slice $f(x,x)$ as a nonperiodic chebfun, and {\tt trace} is
% the integral of $f(x,x)$ over its domain.

d = diag( f );
plot( d )
title( 'The diagonal slice of f' ), snapnow
trace_f = trace( f ) % The trace is the integral of f(x,x) over its domain  
sum( d )

%%
% Just like the rest of Chebfun, Diskfun is designed to perform operations at
% essentially machine precision, and using Diskfun requires no special
% knowledge about the underlying algorithms or discretization procedures.
% Those interested in such details can find an in-depth description of
% how Diskfun works in [3].

%% Basic operations
% A suite of commands are available in Diskfun, and here we describe only 
% a few. A complete listing can be found by typing {\tt methods diskfun} 
% in the MATLAB command line.

%%
% We start by adding, subtracting, and multiplying diskfuns together:

g = diskfun(@(th, r) -40*(cos(((sin(pi*r).*cos(th)...
    + sin(2*pi*r).*sin(th)))/4))+39.5, 'polar');
f = diskfun( @(x,y) cos(15*((x-.2).^2+(y-.2).^2))...
    .*exp(-((x-.2)).^2-((y-.2)).^2));

plot( g )
title( 'g' )
snapnow

plot( f )
title( 'f' )
snapnow

h = g + f;
plot(h)
title( 'g + f' )
snapnow

h = g - f;
plot( h )
title( 'g - f' )
snapnow

h = g.*f;
plot( h )
title( 'g x f' )


%%
% In addition to algebraic operations, we can also solve unconstrained 
% global optimzation problems. In this example, we use the command 
% {\tt max2} to plot $f$ along with its maximum and minimum values.

[jM, kM] = max2( f )
[jm, km] = min2( f )
plot( f ), hold on
colorbar
plot3(kM(1), kM(2), jM, 'k.', MS, 30);
plot3(km(1), km(2), jm, 'r.', MS, 30);
hold off

%%
% There are many ways to visualize a function on the disk and Diskfun offers
% several different ways. For example, here is a contour plot of $g$, with 
% the zero contours displayed as black lines:

contour(g, 'Linewidth', 1.2), hold on
colorbar
contour(g, [0 0], '-k', 'Linewidth', 2)
hold off

%%
% The roots of a function (1D contours) can also be found explicitly and 
% computed with as functions. The contours are stored as a cell array of 
% chebfuns. Each cell consists of two chebfuns that parametrize the contour.

r = roots(g);
plot(g), hold on
for j = 1:length(r)
    rj = r{j};
    plot(r{j}(:,1), r{j}(:,2), '-k', LW, 3)
end
hold off

%%
% One can also perform calculus on diskfuns. For instance, the integral of 
% the function $g(x,y) = -x^2 - 3xy-(y-1)^2$ over the unit disk can be 
% computed using the {\tt sum2} command. We know that the exact answer is
% $-3\pi/2$, which we use to compare the accuracy against.

f = diskfun(@(x,y) -x.^2 - 3*x.*y-(y-1).^2)
sum2(f)
-3*pi/2

%%
% Differentiation on the disk with respect to the polar variable $\rho$
% can lead to singularities, even for smooth functions.
% For example, the function $f( \theta, \rho) = \rho^2$ is smooth
% on the disk, but $\partial f/ \partial \rho = 2 \rho$ is not
% smooth. For this reason, differentiation in Diskfun is only done with 
% respect to the Cartesian coordinates, $x$ and $y$.

%%
% Here, we examine a  harmonic conjugate pair of functions, $u$ and $v$.
% We can use Diskfun to check that they satisfy the Cauchy-Riemann equations.
% Geometrically, this implies that the contour lines of $u$ and $v$
% intersect orthogonally. We can inspect this.

u = diskfun(@(x,y) exp(x).*sin(y));
v = diskfun(@(x,y) -exp(x).*cos(y));
dyu = diffy(u);  % u_y
dxv = diffx(v);  % u_x

norm(dyu+dxv)             % Check u_y =- v_x
norm(diffx(u) - diffy(v)) % Check u_x = v_y

contour(u, 30, 'b'), hold on
contour(v, 30, 'm')
hold off
title( 'Contour lines for u and v' )
snapnow

%%
% In the next example, we compute the derivatives of a function
% involving the Bessel function.

f = diskfun(@(x,y) besselj(0, 5*y).*besselj(0, 5*(x-.1)).*exp(-x.^2-y.^2));
plot( f )
axis off
snapnow
title( 'f' )

plot( diffx( f ) )
axis off
title( 'derivative of f with respect to x' )
snapnow

plot( diffy( f ) )
axis off
title( 'derivative of f with respect to y' )
snapnow

plot( laplacian( f ) )
axis off
title( 'Scalar laplacian of f' )
snapnow

%% Poisson's equation and the cylindrical harmonics
% We can use Diskfun to compute solutions to Poisson's equation on the disk.
% In this example, we compute the solution $v(\theta, \rho)$ for Poisson's
% equation with a Dirichlet boundary condition, so that
% $$ \Nabla^2 v = f, \qquad f(\theta, 1) = 1. $$
% Here, $(\theta, \rho) \in [-\pi, \pi] \times [0, 1]$ and
% $f = \sin \left( 21 \pi \left( 1 + \cos( \pi \rho)
% \right) \rho^2-2\rho^5\cos \left( 5(t-.11)\right) \right)$. The solution is
% returned as a diskfun, so we can immediately plot it, evaluate it, 
% find its zero contours, or perform other operations.

f = @(t,r) sin(21*pi*(1+cos(pi*r)).*(r.^2-2*r.^5.*cos(5*(t-0.11))));
rhs = diskfun(f, 'polar'); % Right-hand side
bc = @(t) 0*t+1;           % Boundary condition
v = diskfun.poisson(f, bc, 512) % Solve for u for 512 x 512 dofs

plot( rhs )
axis off
title( 'f' ), snapnow

plot( v )
axis off
title( 'v' ), snapnow

%%
% The accuracy of the solver can be examined by considering the
% {\em cylindrical harmonics}\footnote{Strictly speaking, we
% consider the cylindrical harmonic functions with a fixed height.},
% which are the eigenfunctions of Laplace's equation on the disk. These 
% functions form an orthogonal basis that is analogous to the trigonometric 
% basis for 1D functions on the unit circle. Cylindrical harmonic functions 
% with a fixed height are defined by two parameters; they are of the form
% $$[V_n^\ell(\theta, \rho) = A_n^\ell J_n(\omega_\ell \rho)e^{in\theta}$$,
% where $A_n^\ell$ is a normalization constant, $J_n$ is an $n$th order
% Bessel function, and $\omega_\ell$ is the $\ell$th zero of $J_n$.
% They are easily constructed in Diskfun:

v = diskfun.harmonic(4, 5); % n = 5, ell = 5
plot(v)
title('Cylindrical harmonic with parameters (4,5)')
view([28,53])
axis off


%%
% Because it is an eigenfunction, $\Nabla^2v = -\lambda v$. The eigenvalues
% are related to the zeros of the Bessel functions. In this case,
% $\sqrt{\lambda} \approx 20.8269329569623$.

%%
bc = @(t) 0*t; %need to look into this
lam = 20.8269329569623;
u = diskfun.poisson(-lam^2*v,bc, 100);

%%
% We expect that $u = v$:
%%
norm(u-v)

%% Vector calculus
% Since the introduction of Chebfun2, Chebfun has supported computation with
% vector-valued functions, including functions in $2D$ (Chebfun2v), $3D$
% (Chebfun3v), and spherical geometries (Spherefunv). Similarly, Diskfunv
% allows one to compute with vector-valued functions on the disk.
% Currently, there are dozens commands available in Diskfunv, including
% vector-based algebraic commands such as {\tt cross} and {\tt dot},
% as well as commands that map vector-valued functions to
% scalar-valued functions (e.g., {\tt curl}, {\tt div} and {\tt jacobian})
% and vice-versa (e.g., {\tt grad}), and commands for performing calculus
% with vector fields (e.g., {\tt laplacian}).

%%
% To get started, we create a diskfun representing a function
% $f = \psi + \phi$, and then compute its gradient.
% The result is returned as a vector-valued object called a diskfunv, with 
% a lower case ``D".

psi = @(x,y) 10*exp(-10*(x+.3).^2-10*(y+.5).^2)+10*...
    exp(-10*(x+.3).^2-10*(y-.5).^2)+ 10*(1-x.^2-y.^2)-20;
phi = @(x,y) 10*exp(-10*(x-.6).^2-10*(y).^2);

f = diskfun( psi ) + diskfun( phi );
u = grad( f )

%%
% The vector-valued function $\mathbf{u}$ consists of two components,
% ordered with respect to unit vectors in the directiosn of $x$ and $y$,
% respectively. Each of these is stored as a diskfun. The Cartesian
% coordinate system is used because this ensures that each component is a
% smooth function on the disk. The unit vectors in the polar and radial
% directions are discontinuous at the origin of the disk, and working
% with them can lead to singularities.  We can view the vector field using 
% a quiver plot:

plot(f)
hold on
quiver(u, 'k')
hold off

%%
% Once a diskfunv object is created, we have overloaded dozens of commands 
% to compute with them. For example, here is a contour plot of the 
% divergence for $\mathbf{u}$.

D = div( u );
contour(D, 'Linewidth', 1.5), hold on
quiver(u, 'k')
hold off

%%
% Since $\mathbf{u}$ is the gradient of $f$, we can verify that
% $\Nabla \cdot \mathbf{u} = \Nabla^2 f$:

norm( div(u) - lap(f) )

%%
% Additionally, since $\mathbf{u}$ is a gradient field,
% \Nabla \cross \mathbf{u} =  0$. We can verify this with the {\tt curl}
% command.

v = curl(u);
norm( v )

%%
% Diskfunv objects can be created through calling the constructor directly with
% function handles or diskfuns input for each component, or by vertically
% concatenating two diskfuns. Here, we demonstrate this by forming a 
% diskfunv $\mathbf{V}$ that represents the numerical surface curl for a 
% scalar-valued function $g$, which can be interpreted as 
% $\Nabla \times [0, 0, g]$.

g = diskfun( @(x,y) cosh(.25.*(cos(5*x)+sin(4*y.^2)))-2 );
dgx = diffx( g );
dgy = diffy( g );

V = diskfunv(dgy, -dgx);  % call constructor
norm(V - [dgy; -dgx])     % equivalent to vertical concatenation

plot( g ), hold on
quiver(V, 'w')
title( 'The numerical surface curl of g' )
hold off

%%
% This construction is equivalent to using the command ${\tt curl}$
% on the scalar function $g$:

norm( V - curl(g) )

%%
% To see all the available commands in the vector part of Diskfun, 
% type {\tt methods diskfunv} in the MATLAB command line.

%% Constructing a diskfun
% The above sections describe how to use Diskfun, and this section provides 
% a brief overview of how the algorithms in Diskfun work. This can be useful 
% for understanding various aspects of approximation involving functions 
% on the disk. More details can be found in [3], and also in the closely 
% related Spherefun chapter (Chapter 17) of the guide.

%%
% Like Chebfun2 and Spherefun, Diskfun uses a variant of Gaussian
% elimination to form low rank approximations to functions. This often results
% in a compressed representation of the function, and it also facilitates
% the use of highly efficient algorithms that work primarily on 1D
% functions.

%%
% To construct a diskfun from a function $f$, we consider an extended version
% of $f$, denoted by $\tilde{f}$, which is formed by taking $f(\theta, \rho)$ 
% and letting $\rho$ range over [-1, 1]$, as opposed to $[0, 1]$.
% This is the disk analogue to the so-called double Fourier sphere method 
% discussed in Chapter 17. Also, see [1,4]. The function $\tilde{f}$ has a 
% special structure, referred to as a block-mirror-centrosymmetric (BMC) 
% structure.  By forming approximants that preserve the BMC structure of 
% $\tilde{f}$, smoothness near the origin is guaranteed. To see the BMC 
% structure, we construct a diskfun $f$ and use the {\tt pol2cart} command:

f= @(th, r) -cos(((sin(pi*r).*cos(th) + sin(2*pi*r).*sin(th)))/4);
f = diskfun(f, 'polar');
%tf = pol2cart(f, 'dfs') %to do: thinking about whether a doubled version
%should be easily accesible (it just does CDR)
%plot(tf)
%title('The BMC function associated with f')

%%
% A structure-preserving method of Gaussian elimination (GE) (see [3])
% adaptively selects a collection of 1D circular and radial
% ``slices" that are used to approximate $\tilde{f}$. Each circular slice
% is a periodic function in $\theta$, and is represented by a trigonometric
% interpolant (or trigfun, see Chapter 11). Each radial slice, a function in
% $\rho$, is represented as a chebfun.  These slices form a low rank
% representation of $f$,
% $$ f(\theta, \rho) \approx \sum_{j=1}^{n} d_jc_j(\rho)r_j(theta)$$,
% where $\{d_j\}_{j=1}^{n}$ are pivot values associated with the GE
% procedure.

%%
% The {\tt plot} command can be used to display the ``skeleton" of $f$:
% the locations of the slices that were adaptively selected and sampled
% during the GE procedure.
clf
plot(f, '.-', MS, 20)

%% References
%%
% [1] B. Fornberg, A Practical Guide to Pseudospectral Methods, 
% Cambridge University Press, 1998.
%%
% [2] A. Townsend, H. Wilber, and G.B. Wright, Computing with functions in
% spherical and polar geometries I. The sphere, SISC, to appear, (2016).
%%
% [3] A. Townsend, H. Wilber, and G.B. Wright, Computing with functions in
% spherical and polar geometries II. The disk, submitted, (2016).
%%
% [4] L. N. Trefethen, Spectral methods in MATLAB, SIAM, 2000. 