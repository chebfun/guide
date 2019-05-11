%% 3. Rootfinding and Minima and Maxima
% Lloyd N. Trefethen, October 2009, latest revision May 2019

%% 3.1 |roots|
% Chebfun comes with a global rootfinding capability -- the ability to find
% all the zeros of a function in its region of definition.  For example,
% here is a polynomial with two roots in $[-1,1]$:
  x = chebfun('x');
  p = x^3 + x^2 - x;
  r = roots(p)

%%
% We can plot $p$ and its roots like this:
  plot(p), grid on
  hold on, plot(r,p(r),'.r')

%%
% Of course, one does not need Chebfun to find roots of a polynomial. The
% MATLAB |roots| command works from a polynomial's coefficients and computes
% estimates of all the roots, not just those in a particular interval.
  roots([1 1 -1 0])

%%
% A more substantial example of rootfinding involving a Bessel function was
% considered in Sections 1.2 and 2.4.  Here is a similar calculation for
% the Airy functions Ai and Bi, modeled after the page on Airy functions
% at WolframMathWorld.
  Ai = chebfun(@(x) airy(0,x),[-10,3]);
  Bi = chebfun(@(x) airy(2,x),[-10,3]);
  hold off, plot(Ai,'r')
  hold on, plot(Bi,'b')
  rA = roots(Ai); plot(rA,Ai(rA),'.r')
  rB = roots(Bi); plot(rB,Bi(rB),'.b')
  axis([-10 3 -.6 1.5]), grid on

%%
% Here for example are the three roots of Ai and Bi closest to 0:
  [rA(end-2:end) rB(end-2:end)]

%%
% Chebfun finds roots by a method due to Boyd and Battles [Boyd 2002,
% Boyd 2014, Battles 2006].  If the
% chebfun is of degree greater than about $50$, it is
% broken into smaller pieces recursively.  On each small piece zeros are
% then found as eigenvalues of a "colleague matrix", the analogue for
% Chebyshev polynomials of a companion matrix for monomials [Specht 1960,
% Good 1961]. This method can be startlingly fast and accurate.  For
% example, here is a sine function with $11$ zeros:
  f = chebfun('sin(pi*x)',[0 10]);
  lengthf = length(f)
  tic, r = roots(f); toc
  r

%%
% A similar computation with 101 zeros comes out equally well:
  f = chebfun('sin(pi*x)',[0 100]);
  lengthf = length(f)
  tic, r = roots(f); toc
  fprintf('%22.14f\n',r(end-4:end))

%%
% And here is the same on an interval with 1001 zeros.
  f = chebfun('sin(pi*x)',[0 1000]);
  lengthf = length(f)
  tic, r = roots(f); toc
  fprintf('%22.13f\n',r(end-4:end))

%%
% Here are the 130 roots from the "fish fillet" example in
% the Chebfun gallery:
f = cheb.gallery('fishfillet');
hold off, plot(f)
tic, r = roots(f); toc
hold on, plot(r,f(r),'.r'), hold off

%%
% With the ability to find zeros, we can solve a variety of nonlinear
% problems.  For example, where do the curves $x$ and $\cos(x)$ intersect?  Here
% is the answer.
  x = chebfun('x',[-2 2]);
  hold off, plot(x)
  f = cos(x);
  hold on, plot(f,'k')
  r = roots(f-x)
  plot(r,f(r),'or')

%%
% All of the examples above concern chebfuns consisting of a single fun.
% If there are several funs, then roots are included at jumps as necessary.
% The following example shows a chebfun with $26$ pieces.  It has nine zeros:
% one within a fun, eight at jumps between funs.
  x = chebfun('x',[-2 2]);
  f = x^3 - 3*x - 2 + sign(sin(20*x));
  hold off, plot(f), grid on
  r = roots(f);
  hold on, plot(r,0*r,'.r')
  
%%
% If one prefers only the "genuine" roots, omitting those at jumps, they
% can be computed by using the `nojump` flag:
  hold off, plot(f), grid on
  r = roots(f,'nojump');
  hold on, plot(r,0*r,'.r')

%% 3.2 |min|, |max|, |abs|, |sign|, |round|, |floor|, |ceil|
% Rootfinding is more central to Chebfun than one might at first imagine,
% because a number of commands, when applied to smooth chebfuns, must
% produce non-smooth results, and it is rootfinding that tells us where to
% put the discontinuities. For example, the |abs| command introduces
% breakpoints wherever the argument goes through zero.  Here we see that |x|
% consists of a single piece, whereas |abs(x)| consists of two pieces.
  x = chebfun('x')
  absx = abs(x)
  subplot(1,2,1), plot(x)
  subplot(1,2,2), plot(absx)

%%
% We saw this effect already in Section 1.4.
% Another similar effect shown in that section occurs with
% |min(f,g)| or |max(f,g)|.  Here, breakpoints are introduced
% at points where |f-g| is zero:
  f = min(x,-x/2), subplot(1,2,1), plot(f)
  g = max(.6,1-x^2), subplot(1,2,2), plot(g), ylim([.5,1])

%%
% The function |sign| also introduces breaks, as illustrated in the last
% section. The commands |round|, |floor|, and |ceil| behave like this too.
% For example, here is $\exp(x)$ rounded to nearest multiples of $0.1$:
  g = exp(x);
  clf, plot(g)
  gh = round(10*g)/10;
  hold on, plot(gh,'jumpline','-');
  grid on

%% 3.3 Local extrema
% Local extrema of smooth functions can be located by finding
% zeros of the derivative.  For example, here is a variant of
% the Airy function again, with all its extrema in its range of
% definition located and plotted.
  f = chebfun('exp(real(airy(x)))',[-15,0]);
  clf, plot(f)
  r = roots(diff(f));
  hold on, plot(r,f(r),'.r'), grid on

%%
% Chebfun users don't have to compute the derivative explicitly
% to find extrema, however.  An alternative is to type
  [ignored,r2] = minandmax(f,'local');

%%
% which returns both interior local extrema and also
% the endpoints of |f|.
% Similarly one can type |min(f,'local')| and |max(f,'local')|.

%%
% These methods will find non-smooth extrema as well as smooth ones,
% since these correspond to "zeros" of the derivative where
% the derivative jumps from one sign to the other.  Here is an
% example.
  x = chebfun('x');
  f = exp(x)*sin(30*x);
  g = 2-6*x^2;
  h = max(f,g);
  clf, plot(h)

%%
% Here are all the local extrema, smooth and nonsmooth:
  [ignored,extrema] = minandmax(h,'local');
  hold on, plot(extrema,h(extrema),'.r')

%%
% Suppose we want to pick out the local extrema that are actually local
% minima.  We can do that by hand by checking for the second derivative to be
% positive:
  h2 = diff(h,2);
  maxima = extrema(h2(extrema)>0);
  plot(maxima,h(maxima),'ok')

%%
% Or we could do it implicitly with |local|,
  [maxval,maxpos] = min(h,'local')
  plot(maxpos,maxval,'.k')

%% 3.4 Global extrema: max and min
% If |min| or |max| is applied to a single chebfun, it returns its global
% minimum or maximum.  For example:
  f = chebfun('1-x^2/2');
  [min(f) max(f)]

%%
% Chebfun computes such a result by checking the values of |f| at
% endpoints and at zeros of the derivative.

%%
% As with standard MATLAB, one can find the location of the extreme point
% by supplying two output arguments:

  [minval,minpos] = min(f)
%%
% Note that just one position is returned even though the minimum is
% attained at two points.  This is consistent with the behavior of standard
% MATLAB.

%%
% This ability to do global 1D optimization in Chebfun is rather
% remarkable.  Here is a nontrivial example.
  f = chebfun('sin(x)+sin(x^2)',[0,15]);
  hold off, plot(f,'k')
%%
% The length of this chebfun is not as great as one might imagine:
  length(f)
%%
% Here are its global minimum and maximum:
  [minval,minpos] = min(f)
  [maxval,maxpos] = max(f)
  hold on
  plot(minpos,minval,'.b')
  plot(maxpos,maxval,'.r')

%%
% For larger chebfuns, it is inefficient to compute the global minimum and
% maximum separately like this -- each one must compute the derivative and
% find all its zeros. The alternative |minandmax| code mentioned above
% provides a faster alternative:
  [extremevalues,extremepositions] = minandmax(f)

%% 3.5 |norm(f,1)| and |norm(f,inf)|
% The default, $2$-norm form of the |norm| command was considered in Section
% 2.2. In standard MATLAB one can also compute $1$-, $\infty$-, and Frobenius
% norms with |norm(f,1)|, |norm(f,inf)|, and |norm(f,'fro')|.  These have been
% overloaded in Chebfun, and in the first two cases, rootfinding is part of
% the implementation.  (The $2$- and Frobenius norms are equal for a single
% chebfun but differ for quasimatrices; see Chapter 6.) The $1$-norm
% |norm(f,1)| is the integral of the absolute value, and Chebfun computes
% this by adding up segments between zeros, at which $|f(x)|$ will typically
% have a discontinuous slope. The $\infty$-norm is computed from the
% formula $\|f\|_\infty = \max(\max(f),-\min(f))$.

%%
% For example:
  f = chebfun('sin(x)',[103 103+4*pi]);
  norm(f,inf)
  norm(f,1)

%% 3.6 Roots in the complex plane
% Chebfuns live on real intervals, and the funs from which they are made
% live on real subintervals.  But a polynomial representing a fun may have
% roots outside the interval of definition, which may be complex.
% Sometimes we may want to get our hands on these roots, and the |roots|
% command makes this possible in various ways through the flags |'all'|,
% |'complex'|, and |'norecursion'|.

%%
% The simplest example is a chebfun that is truly intended to
% correspond to a polynomial.  For example, the chebfun
f = 1+16*x^2;

%%
% has no roots in $[-1,1]$:
roots(f)

%%
% We can extract its complex roots with the command
roots(f,'all')

%%
% The situation for more general chebfuns is more complicated.  For example,
% the chebfun
g = exp(x)*f(x);

%% 
% also has no roots in $[-1,1]$,
roots(g)

%%
% but it has plenty of roots in the complex plane:
roots(g,'all')

%%
% Most of these are spurious. What has happened is that |g| is represented by
% a polynomial chosen for its approximation properties on $[-1,1]$.
% Inevitably that polynomial will have roots in the complex plane, even if
% they have little to do with |g|.
% (See the discussion of the Walsh and Blatt-Saff theorems
% in Chapter 18 of [Trefethen 2013].)

%%
% One cannot expect Chebfun to solve this problem perfectly -- after all,
% it is working on a real interval, not in the complex plane, and analytic
% continuation from the one to the other is well known to be an ill-posed
% problem.  Nevertheless, Chebfun may do a pretty good job of selecting
% genuine complex (and real) roots near the interval of definition if you
% use the |'complex'| flag:
roots(g,'complex')

%%
% We will not go into detail here of how this is done, but the idea is that
% associated with any fun is a family of "Chebfun ellipses" in the complex
% plane, with foci at the endpoints, inside which one may expect reasonably
% good accuracy of the fun.
% Assuming the interval is $[-1,1]$ and the fun
% has length $L$, the Chebfun ellipse associated with a parameter
% $\delta\ll 1$ is
% the region in the complex plane bounded by the image under the map
% $(z+1/z)/2$ of the circle $|z|=r$, where $r$ is defined by
% the condition $r^{-L}=\delta$.  (See
% Chapter 8 of [Trefethen 2013].)  The command |roots(g,'complex')| first
% effectively does |roots(g,'all')|, then returns only those roots lying
% inside a particular Chebfun ellipse -- we take the one corresponding to
% delta equal to the square root of the Chebfun tolerance
% parameter |eps|.

%%
% One must expect complex roots of chebfuns to lose accuracy as one moves
% away from the interval of definition.  Here's an example:
  f = exp(exp(x)).*(x-0.1i).*(x-.3i).*(x-.6i).*(x-1i);
  roots(f,'complex')

%%
% Notice that the first three imaginary roots are located with about
% $13$, $10$, and $6$ digits of accuracy,
% while the fourth does not appear in the list at all.

%%
% Here is a more complicated example:
F = @(x) 4+sin(x)+sin(sqrt(2)*x)+sin(pi*x);
f = chebfun(F,[-100,100]);

%%
% This function has a lot of complex roots lying in strips on either side
% of the real axis.
  r = roots(f,'complex');
  hold off, plot(r,'.')

%%
% If you are dealing with complex roots of complicated chebfuns like this,
% it may be safer to add a flag |'norecursion'| to the roots call, at the cost
% of slowing down the computation. This turns off the Boyd-Battles
% recursion mentioned above, lowering the chance of missing a few roots
% near interfaces of the recursion.  If we try that here we find that the
% results look almost the same as before in a plot:
r2 = roots(f,'complex','norecursion');
hold on, plot(r,'om')

%%
% However, the accuracy has improved:
norm(F(r))
norm(F(r2))

%%
% To find poles in the complex plane as opposed to zeros, see Section 4.8 
% and also [Austin, Kravanja & Trefethen 2015].
% More advanced methods of rootfinding and polefinding are based on
% rational approximations rather than polynomials, an area where Chebfun
% has significant capabilities; see the next chapter of this guide,
% Chapter 28 of [Trefethen 2013], and [Webb 2013].

%% 3.7 References
%
% [Austin, Kravanja & Trefethen 2015] A. P. Austin,
% P. Kravanja, and L. N. Trefethen, "Numerical algorithms
% based on analytic function values at roots of unity",
% _SIAM Journal on Numerical Analysis_, to appear.
%
% [Battles 2006] Z. Battles, _Numerical Linear Algebra for
% Continuous Functions_, DPhil thesis, Oxford University
% Computing Laboratory, 2006.
%
% [Boyd 2002] J. A. Boyd, "Computing zeros on a real interval
% through Chebyshev expansion and polynomial rootfinding",
% _SIAM Journal on Numerical Analysis_, 40 (2002), 1666-1682.
%
% [Boyd 2014] J. A. Boyd, _Solving Transcendental Equations:
% The Chebyshev Polynomial Proxy and Other Numerical
% Rootfinders, Perturbation Series, and Oracles_, SIAM, 2014.
%
% [Good 1961] I. J. Good, "The colleague matrix, a Chebyshev
% analogue of the companion matrix", _Quarterly Journal of 
% Mathematics_, 12 (1961), 61-68.
%
% [Specht 1960] W. Specht, "Die Lage der Nullstellen eines Polynoms.
% IV", _Mathematische Nachrichten_, 21 (1960), 201-222.
%
% [Trefethen 2013] L. N. Trefethen, _Approximation Theory and
% Approximation Practice_, SIAM, 2013.
%
% [Webb 2013] M. Webb, "Computing complex singularities of
% differential equations with Chebfun",
% _SIAM Undergraduate Research Online_, 6 (2013), 
% http://dx.doi.org/10.1137/12S011520.
