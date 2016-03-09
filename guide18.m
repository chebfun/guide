%% 18. Chebfun3
% Behnam Hashemi and Nick Trefethen, November 2015

%%
% NOTE: Chebfun3 is not yet publicly available,
% but we expect to release it during 2016.  Some of
% the details described below may change.

%% 18.1  Introduction
% The Chebfun project began in 2003, and Chebfun2, for 2D functions, was first
% released in 2013. Chebfun3, for 3D functions, has been created by Behnam
% Hashemi and will be released sometime in 2016.

%%
% Chebfun3 aims to compute with functions in a 3D box
% $[a,b]\times[c,d]\times [e,g]$, which by default is the cube $[-1,1]^3$.
% Our aim in creating it has been to make it analogous to Chebfun2 wherever
% possible, but of course, there are some differences in capabilities and
% in underlying algorithms and representations.  Like Chebfun2, Chebfun3
% can carry out a wide range of computations on functions and in particular
% is good at integration, differentiation, and computation of minima and
% maxima.  

%%
% For a simple starting example, here is a 3D Runge function:
f = chebfun3(@(x,y,z) 1./(1+x.^2+y.^2+z.^2));

%%
% At (0,.5,.5), it takes the value $2/3$:
format long, f(0,.5,.5)

%%
% The triple integral over the whole box is computed by `sum3`,
sum3(f)

%% 
% The exact result is $4.28685406230184188268\dots.$

%%
% Since the volume of the box is $8$, the mean value is
% $1/8$ of this result:
mean3(f)

%%
% The maximum value is $1$:
%
max3(f)

%%
% The most basic Chebfun3 plotting command is `slice`, which by default
% shows the function on various slices in the three directions.
slice(f)

%%
% Another plotting capability is `isosurface` which, by default, plots 
% three isosurfaces marking quartiles:
clf, isosurface(f)

%%
% So far, there are about 70 methods that can be applied to chebfun3 objects. For a
% complete list type `methods chebfun3`.
% (By the way, notice that in print we use the plural form "chebfun3 objects",
% because the expression "chebfun3s" is confusing, though informally
% in conversation we speak of "chebfun3's".)

%% 18.2 Anatomy of a chebfun3
% First, a quick reminder. In 1D, Chebfun represents functions by
% polynomials and piecewise polynomials, all derived via interpolation
% in Chebyshev points.  In 2D, Chebfun2 does not have piecewise 
% representations, only global ones. These are constructed in a low-rank
% fashion: a function $f(x,y)$ is approximated by functions of 
% rank $1,2,3,\dots$ until approximately 16-digit precision is achieved 
% [Townsend & Trefethen 2013b]. Each rank 1 piece in this representation 
% is an outer product $d c(y) r(x)$, where $d$ is a scalar, and the 
% univariate functions $c$ and $r$ are represented as chebfuns.  At each 
% step of the increasing-rank process, the "pivot location" in the domain 
% rectangle that determines the choice of $c$ and $r$ is determined by an
% approximation to the largest remaining function value in the rectangle.  
% This is analogous to an approximate form of complete pivoting in Gaussian 
% elimination for the representation of a matrix via ranks $1,2,3,\dots$ [Townsend & Trefethen 2013a].

%%
% There is no unique obvious generalization of this 2D idea to 3D or
% higher dimensions.  On the contrary, there are at least half a dozen
% ideas in the literature. (This "literature" is mostly about discrete
% tensors rather than multivariate functions [Bebendorf 2011], but one paper about functions
% worthy of note is [Gorodetsky, Karaman and Marzouk 2015].)  Thus we faced
% some fundamental decisions in the design of Chebfun3, and along the way
% we explored a number of different possibilities, ranging from a
% straightforward 3D tensor product to various compressed-rank
% representations.

%%
% Chebfun3 represents functions in the _Tucker format_, which means that
% $f(x,y,z)$ is written as a 3D product involving an
% $\infty\times r_1$ quasimatrix $c$ ("columns") in the $x$ direction, an
% $\infty\times r_2$ quasimatrix $r$ ("rows") in the $y$ direction, and an
% $\infty\times r_3$ quasimatrix $t$ ("tubes")
% in the $z$ direction.  Here $r_1, r_2, r_3$
% are positive integers whose size might typically be 10 or 20 for
% the functions that Chebfun3 can efficiently represent.  These
% three quasimatrices are combined in a product with coefficients
% specified by an $r_1\times r_2\times r_3$ _core tensor_ called `core`.  
% In short, we write:
% $$ f(x,y,z) \approx core \times_1 c(x) \times_2 r(y) \times_3
% t(z), $$ or more fully:
%
% $$ f(x,y,z) \approx \sum_{i=1}^{r_1} \sum_{j=1}^{r_2} 
% \sum_{k=1}^{r_3} core(i,j,k) c_i(x) r_j(y) t_k(z). $$

%%
% We can get some information about the Tucker representation of
% our 3D Runge function $f$ by typing `f` without a semicolon:
f

%%
% Let's look at a few even simpler examples.  Suppose that $f$ happens
% to depend only on $x$, like this:
f = chebfun3(@(x,y,z) exp(x))

%%
% From this output we see that $f$ has rank 1 in some sense
% (the word "rank" has no unique definition for 3D representations).
% To be precise, the core tensor in this case only needs to
% be a $1\times 1\times 1$ scalar, and the same would apply
% for any function that depends on just one of $x, y$, or $z$.
% To find out what degree polynomial is used to capture the
% dependences in the three directions we can use `length`
% [whether this is the right name is under discussion -- it
% differs from Chebfun2 usage]:
[m, n, p] = length(f)

%%
% Similarly here is what it looks like if the function is
% $\exp(y)$ instead of $\exp(x)$:
[m, n, p] = length(chebfun3(@(x,y,z) exp(y)))

%%
% These results compare in the expected way with
% the 1D standard chebfun of $\exp(x)$.
length(chebfun(@(x) exp(x)))

%%
% To see a little more of the structure of a chebfun3, let
% us cook up an example with 2 columns, 2 rows, and 2 tubes.
% Here is such a function:
f = chebfun3(@(x,y,z) exp(x).*(log(2+y).*exp(z))+sin(y)/1e6)

%%
% The representation of $f$ involves an $2\times 2\times 2$
% core tensor, an $\infty\times 2$
% quasimatrix for the $x$ direction, and also $\infty\times 2$
% quasimatrices for the $y$ and $z$ directions. Here is the
% core tensor:
format short, f.core

%%
% Here are plots of the coefficients of the three
% quasimatrices:
plotcoeffs(f,'.-')

%%
% These explorations give an idea of what a chebfun3 looks like.
% However, they don't explain how the system constructs such an
% object.  We will not give details here; a paper is in preparation 
% [Hashemi & Trefethen 2016]. Very briefly, we use what we call a "slice" 
% decomposition as an intermediate step: chebfun3 constructs a sum in which
% each term is a product of a chebfun in one direction with a chebfun2 in 
% the other two directions.  Thus the construction process does not
% treat $x$, $y$, and $z$ exactly symmetrically.  At the end, however,
% the constructed representation is converted to the Tucker
% form described above.

%% 18.3  Computing with chebfun3 objects
% Of course, Chebfun is all about computing with functions, not
% just representing them.  For example, here are two 3D functions:
f = chebfun3(@(x,y,z) sin(x+y.*z));
g = chebfun3(@(x,y,z) cos(15*exp(z))./(5+x.^3+2*y.^2+z));

%%
% Here are the maxima of $f$, $g$, and $fg$:
format long
max3(f)
max3(g)
max3(f.*g)

%%
% Of course, one can find the location as well as the
% value of an extremum:
[maxval,maxpos] = max3(f+g)

%%
% Here is the integral over the cube of $f\exp(g)$:
sum3(f.*exp(g))

%%
% If we execute just `sum`, it integrates over just one
% dimension, by default $x$, so the output is a 2D function, i.e.,
% a chebfun2:
close all, contourf(sum(f),20), colorbar

%%
% There is also a `sum2` command for integration over two
% dimensions, by default $x$ and $y$, giving as output a
% 1D function, i.e., a chebfun:
plot(sum2(exp(g+2*f)))  

%%
% Here is a line integral over a 3D spiral
curve = chebfun(@(t) [cos(t) sin(t) t/(8*pi)], [0, 8*pi]);
close all, plot3(curve(:,1), curve(:,2), curve(:,3) ), title('Helix')
f = chebfun3(@(x,y,z) x+y.*z );
I = integral( f, curve )
exact = -sqrt(1+(8*pi)^2)/(8*pi)

%% 18.4 Getting inside a chebfun3
% Suppose $f$ is a chebfun3.
% We can examine its columns, rows, and tubes by executing
% `f.cols`, `f.rows`, and `f.tubes`.  For example, let us look at the
% columns associated with the chebfun3 $g$ just considered.
% This is a quasimatrix with 9 columns:
size(g.cols)

%%
% Here is a plot of the columns:
plot(g.cols)

%%
% The tubes are more interesting (also the rows):
plot(g.tubes)

%%
% A plot of _all_ the coefficients of $g$ looks like this.
plotcoeffs(g,'.-')

%%
% If we look just at the columns we get the corresponding piece:
clf, plotcoeffs(g.cols,'.-')

%% 18.5  Periodic chebfun3 objects 
% Chebfun3 can use trigonometric functions instead of polynomials for 
% representing smooth functions which are triply periodic. 
% To create a trig-based chebfun3 object, we can use the 'trig' (or 'periodic') 
% flag in the Chebfun3 constructor. For example, the function 
% $f(x,y,z) = \tanh(3\sin x) - \sin(y+1/2) + \cos(6z)$ on $[-\pi, \pi)^3$ 
% can be constructed as follows:
f = chebfun3(@(x,y,z) tanh(3*sin(x))-(sin(y+1/2)).^2+cos(6*z), [-pi pi -pi pi -pi pi], 'trig');

%% 
% The text 'trig' in the following display of columns of $f$ indicates
% representation in terms of trigonometric functions:
f.cols

%% 
% Here is the length of $f$ and a plot of its coefficients:
[m, n, p] = length(f)
plotcoeffs(f,'.-')

%%
% As we see, $f$ is resolved to machine precision using trigonometric 
% interpolants through 143, 5, and 13 equally spaced samples over $[-\pi,\pi)^3$
% in the $x$, $y$, and $z$ directions, respectively. In other words,
% trigonometric polynomials of degrees 71, 2, and 6 are enough to resolve 
% this function. Let's compare with the length of $f$ if we ignore its periodicity:
fCheb = chebfun3(@(x,y,z) tanh(3*sin(x))-sin(2*y+1/2)+cos(6*z), [-pi pi -pi pi -pi pi]);
[m_fCheb, n_fCheb, p_fCheb] = length(fCheb)
%% 
% As expected, a smooth periodic function can be represented with 
% trigfun factor quasimatrices using fewer samples than standard chebfun 
% factor quasimatrices.

%% 18.6 Derivative and Laplacian
% Like Chebfun and Chebfun2, Chebfun3 is good at calculus,
% having commands `diffx`, `diffy`, `diffz`, and `lap` (identical
% to `laplacian`).  There is also a general command `diff` that
% can take appropriate arguments to specify dimensions.

%%
% For example, the following function is a harmonic function:
f = chebfun3(@(x,y,z) 1./sqrt(x.^2 + y.^2 + (2-z).^2));

%%
% So its Laplacian will be the constant 2:
Lf = lap(f);
Lf(rand(3,1),rand(3,1),rand(3,1))
%%
% Let's compare the Laplacian of $f$ with the divergence of its gradient,
% which will be discussed in the next subsection.
norm(Lf - div(grad(f)))

%% 18.7 3D vector fields
% Consider a vector-valued function of three variables like 
%
% $$F(x,y,z) = (f(x,y,z); g(x,y,z); h(x,y,z)).$$ 
%
% We can represent such functions using `chebfun3v`. Similar to `chebfun2v`, 
% we can construct a chebfun3v object, either by explicitly calling the 
% constructor or by vertical concatenation of chebfun3 objects. We already
% hinted at this by using `grad(f)' in the last subsection. Let's check its size:
size(grad(f))

%%
% Chebfun can plot 3D vector fields using the `quiver3` command, which 
% draws a bunch of arrows. It first sets up a grid of points for 
% which to plot arrows and then draws the vector field. Suppose we wish
% to plot the vector field $F(x,y,z) = -yi + xj + zk$. To do this we type
F = [chebfun3(@(x,y,z) -y); chebfun3(@(x,y,z) x); chebfun3(@(x,y,z) z)];
close all, quiver3(F,0)
view([2 2 40])

%%
% The `0' used in the quiver3 above just ensures that MATLAB does not 
% rescale the vectors.

%%
% According to the fundamental theorem of calculus for line integrals, also 
% known as the gradient theorem, the line integral of a gradient vector field 
% along a smooth curve depends only on the endpoints of the path.
f = chebfun3(@(x,y,z) sin(x+20*y+z.^2).*exp(-(3+y.^2)) , [-5*pi, 5*pi, -5*pi, 5*pi, -5*pi, 5*pi]);
F = grad(f);
curve = chebfun(@(t) [t.*cos(t) t.*sin(t) t], [0, 5*pi]); 
plot3(curve(:,1), curve(:,2), curve(:,3),'b'), 
title('Conical spiral'), shg
I_spiral = integral(F,curve)
ends = f(5*pi*cos(5*pi), 5*pi*sin(5*pi), 5*pi) - f(0, 0, 0)

%% 
% We can determine if a given vector field is conservative using `curl`:
norm( curl(F) )

%%
% Moreover, it is well-known that the line integral of a conservative vector 
% field is independent of path. Let's compare the numerical values of the line 
% integral of our vector field over two paths with the same endpoints:
curve2 = chebfun(@(t) [t*cos(5*pi) t*sin(5*pi) t], [0, 5*pi]);
hold on, plot3(curve2(:,1), curve2(:,2), curve2(:,3),'r-.'), hold off, shg
I_line = integral(F,curve2), I_line - I_spiral

%% 18.8 Rootfinding
% Rootfinding in Chebfun3 is still under development.  What follows
% describes the code `root`, which attempts to find just one single root
% of a chebfun3v or equivalently of three chebfun3 objects.

%%
% One usually expects three 3D functions to have a finite number of common
% roots, though it is possible that such a system of equations would have an 
% infinite number of roots. Given $f, g$, and $h$, rootfinding in Chebfun3
% is done by the following two steps: In step 1, an initial guess is computed 
% from the tensor of values of the chebfun3 object $objFun:= f^2 + g^2 + h^2$. 
% This gives us a root with roughly 3 accurate digits. Step2, then attempts 
% to improve the accuracy to machine precision using a few Newton iterations. 
% [Research needed here!]
% Here is an example: 

f = chebfun3(@(x,y,z) y-x.^2);
g = chebfun3(@(x,y,z) z-x.^3);
h = chebfun3(@(x,y,z) cos(exp(x.*sin(-2+y+z))));

%%
% The intersection of the zero level surfaces of the first two functions $f$
% and $g$ creates a _twisted cubic_ which is an interesting curve introduced in most 
% textbooks on algebraic geometry. See e.g., [Cox, Little & O'Shea 2015].
close all, isosurface(f,0,'g')
hold on, isosurface(g,0,'b')
view([-2,5,5])

%%
% Let us compute the only common root of $f$, $g$ and $h$ in the default cube:
r = root(f,g,h)

%%
% The root is accurate to machine epsilon:
format short
res1 = f(r(1),r(2),r(3))
res2 = g(r(1),r(2),r(3))
res3 = h(r(1),r(2),r(3))

%%
% Here is a plot of all the three zero level surfaces together with the
% common root we just found:
isosurface(h,0,'r')
plot3(r(1),r(2),r(3),'yh','markersize',30)
view([-8,8,5]), alpha(0.9)

%% 18.9 Higher-order SVD
% The higher-order SVD (HOSVD) of a discrete tensor was introduced in 
% [De Lathauwer, De Moor & Vandewalle 2000]. For an order-3 tensor, 
% it uses SVDs of the three modal unfolding matrices to compute a
% factorization involving
% a core tensor with the three matrices of left 
% singular vectors of the unfolded tensor. Chebfun3 contains a continuous
% analogue of the HOSVD. Here is an example:
f = chebfun3(@(x,y,z) sin(x+y+z));
[Score,Scols,Srows,Stubes] = hosvd(f);
%%
% First of all, the frontal slices 
% of the core tensor have decreasing Frobenius norms which we may consider as
% a 3D analogue of the decay of singular values of a matrix:
norm(Score(:,:,1),'fro'), norm(Score(:,:,2),'fro')

%%
% Also, the columns of the three factor quasimatrices `Scols`, `Srows`, and 
% `Stubes` are orthonormal. For example, the departure from orthogonality 
% in the columns of the quasimatrix `Scols` is:
norm( eye(size(Scols,2)) - Scols'*Scols )

%%
% Moreover, the core tensor S is _all orthogonal_, i.e., its frontal slices 
% are orthogonal:
norm(Score(:,:,1).*Score(:,:,2)) 

%%
% Its lateral slices are orthogonal:
norm(squeeze(Score(:,1,:).*Score(:,2,:)))
%%
% And also its horizontal slices are orthogonal:
norm(squeeze(Score(1,:,:).*Score(2,:,:)))

%% 18.10 Changing the accuracy with chebfun3eps
% Chebfun has always had a parameter that describes its
% target accuracy, traditionally named `eps` though in the
% process of being renamed to `chebfuneps`.  For 1D computations,
% we do not recommend that users normally change this parameter
% from its factory value of machine precision (unless dealing with
% noisy functions), because the speedups
% to be obtained are usually not very large.  See section
% 8.8 of this Guide.

%%
% In two dimensions, and even more in three dimensions, the
% potential gains from loosening the tolerance become much
% greater.  Many users of Chebfun3 may find, for example, that
% they want to work with 10 digits of accuracy rather than 16 --
% the speedup in many cases is on the order of a factor of 10.
% For this reason, as of very recently, Chebfun allows users
% to set different tolerances `chebfuneps`, `chebfun2eps`, and
% `chebfun3eps` for computations in 1D, 2D and 3D.  At this moment
% of writing -- though it may change! -- the factory values of
% these parameters are all machine epsilon.

%%
% To illustrate the speedups, here we compute chebfun3 objects
% for a Runge function at four different accuracies.  We caution
% users that the actual accuracy achieved may be a digit or two
% less than requested.

disp('   eps     time      rank m-length n-length p-length')
ff = @(x,y,z)1./(0.01+x.^2+y.^2+z.^2);
for ep = 10.^(-16:4:-4)
  f = chebfun3(ff,'eps', ep);
  tic, f = chebfun3(ff,'eps', ep); t = toc;
  r = rank(f);
  [m,n,p] = length(f);
  fprintf('%8.1e  %6.4f %7d %7d %7d %7d\n', ep,t,r,m,n,p)
end 

%%
% As the above example indicates, one way to control the
% accuracy of a chebfun3 construction is by explicitly specifying
% `eps` in a call to the constructor.  Alternatively -- and
% necessarily, if you are doing follow-on operations on previously
% constructed chebfun3 objects -- you can change the parameter
% globally as described in Chapter 8.  For example, let us check the 
% speed and accuracy of a certain computation using the factory value of 
% `chebfun3eps`, i.e., machine precision. It runs slowly, and gives high 
% accuracy:
ff = @(x,y,z) tanh(2*(x+y+z));

tic
f = chebfun3(ff);
g = cos(f+1).^2;
error = g(.5,.6,.7) - cos(ff(.5,.6,.7)+1)^2
toc

%% 
% Now let us set `chebfun3eps` to 1e-10 and run the same computation.  
% Faster and less accurate!
chebfun3eps 1e-10

tic
f = chebfun3(ff);
g = cos(f+1).^2;
error = g(.5,.6,.7) - cos(ff(.5,.6,.7)+1)^2
toc

%% 
% The following command reverts `chebfun3eps` to the factory value:
chebfun3eps factory

%% 18.11 Chebfun3t for pure tensor product comparisons
% Chebfun3, like Chebfun2 before it, exploits low-rank compression
% of functions where possible.  The question of how much one gains
% from this on average is controversial and certainly unresolved.
% One can construct examples where the gain is as large as you like,
% and other examples where there is no gain at all (and indeed where
% the low-rank algorithms take longer).

%%
% To enable interested users to explore these tradeoffs, Chebfun
% offers a certain fraction of the Chebfun3 functionality implemented
% in a completely different, more straightforward but not
% rank-compressed fashion.  The command `chebfun3t` will construct
% a chebfun3t object that is represented as a multivariate Chebyshev
% polynomial defined by a tensor of Chebyshev coefficients.
% Here is an example where Chebfun3 is much faster than Chebfun3t:
ff = @(x,y,z) sin(120*(x+y+z));
tic, chebfun3(ff), toc
tic, chebfun3t(ff), toc

%%
% On the other hand here is an example where we get no benefit:
ff = @(x,y,z) tanh(3*(x+y+z));
tic, chebfun3(ff), toc
tic, chebfun3t(ff), toc

%%
% Let us emphasize that the main tool we are offering for
% computation with functions in 3D is Chebfun3, not Chebfun3t.
% The latter, only partially implemented, as only provided to
% enable certain comparisons.

%% 18.12 References 
% 
% [Bebendorf 2011] M. Bebendorf, "Adaptive cross approximation of 
% multivariate functions", _Constructive Approximation_, 34 (2011) 149-179.
%
% [Cox, Little & O'Shea 2015] D. Cox, J. Little, and D. O'Shea, _Ideals, 
% Varieties, and Algorithms_, 4th Edition, Springer, 2015.
%
% [De Lathauwer, De Moor & Vandewalle 2000] L. De Lathauwer, B. De Moor and 
% J. Vandewalle, "A multilinear Singular Value Decomposition", _SIAM 
% Journal on Matrix Analysis and Applications_, 21 (2000) 1253-1278.
%
% [Golub & Van Loan 2013] G. H. Golub, and C. F. Van Loan, _Matrix Computations_,
% 4th Edition, Johns Hopkins University Press, 2013.
%
% [Gorodetsky, Karaman & Marzouk 2015] A. A. Gorodetsky, S. Karaman, 
% and Y. M. Marzouk, "Function-train: a continuous analogue of the tensor-train 
% decomposition", arXiv:1510.09088v1, 2015.
%
% [Hashemi & Trefethen 2016] B. Hashemi and L. N. Trefethen, "Chebfun in 
% three dimensions", manuscript in preparation.
%
% [Townsend & Trefethen 2013a] A. Townsend and L. N. Trefethen, "Gaussian 
% elimination as an iterative algorithm", _SIAM News_, 46, March 2013.
%
% [Townsend & Trefethen 2013b] A. Townsend and L. N. Trefethen, "An extension
% of Chebfun to two dimensions", _SIAM Journal on Scientific Computing_, 35
% (2013), C495-C518.
