%% 1. Getting Started with Chebfun
% Lloyd N. Trefethen, October 2009, latest revision December 2015

%% 1.1  What is a chebfun?
% A chebfun is a function of one variable defined on an interval $[a,b]$. The
% syntax for chebfuns is almost exactly the same as the usual MATLAB syntax
% for vectors, with the familiar MATLAB commands for vectors overloaded in
% natural ways. Thus, for example, whereas |sum(f)| returns the sum of the
% entries when |f| is a vector, it returns a definite integral when |f| is a
% chebfun.

%%
% Chebfun with a capital C is the name of the software system.

%%
% The aim of Chebfun is to "feel symbolic but run at the speed of numerics".
% More precisely, our vision is to achieve for functions what floating-point
% arithmetic achieves for numbers: rapid computation in which
% each successive operation is carried out exactly apart from a rounding error
% that is very small in relative terms [Trefethen 2007].

%%
% The implementation of Chebfun is based on the mathematical fact that
% smooth functions can be represented very efficiently by polynomial
% interpolation in Chebyshev points, or equivalently, thanks to the Fast
% Fourier Transform, by expansions in Chebyshev polynomials.  For a simple
% function, 20 or 30 points often suffice, but the process is stable and
% effective even for functions complicated enough to require 1000 or
% 1,000,000 points. Chebfun makes use of adaptive procedures that aim to
% find the right number of points automatically so as to represent each
% function to roughly machine precision, that is, about
% 15 digits of relative
% accuracy.  (Originally Chebfun stored function values at Chebyshev
% points; in Version 5 it switched to storing Chebyshev expansion 
% coefficients.)

%%
% The mathematical foundations of Chebfun are for the most part well
% established by results scattered throughout the 20th century.  A key
% early figure, for example, was Bernstein in the 1910s. Much of the
% relevant material can be found collected in the Chebfun-based book
% _Approximation Theory and Approximation Practice_ [Trefethen 2013].

%%
% Chebfun was originally created by Zachary Battles and Nick Trefethen at
% Oxford during 2002-2005 [Battles & Trefethen 2004].  Battles left the
% project in 2005, and soon four new members were added to the team:
% Ricardo Pachon (from 2006), Rodrigo Platte (from 2007), and Toby Driscoll
% and Nick Hale (from 2008). In 2009, Asgeir Birkisson and Mark Richardson
% also became involved, and other contributors included Pedro Gonnet, Joris
% Van Deun, and Georges Klein.  Nick Hale served as Director of the
% project during 2010-2014.  The Chebfun Version 5 rewrite was directed
% by Nick Hale during 2013-2014, and the team included Anthony
% Austin, Asgeir Birkisson, Toby Driscoll, Hrothgar, Mohsin
% Javed, Hadrien Montanelli, Alex Townsend,
% Nick Trefethen, Grady Wright, and Kuan Xu.
% October 2014 brough new arrivals Jared Aurentz, Behnam
% Hashemi, and Mikael Slevinsky.
% Further information about Chebfun history is available at the Chebfun
% web site, http://www.chebfun.org
% where one can also find a discussion of other software projects related
% to Chebfun.
% This Guide is based on Chebfun Version 5.3, released
% in December 2015.

%% 1.2  Constructing simple chebfuns
% The |chebfun| command constructs a chebfun from a specification such as a
% string or an anonymous function.  If you don't specify an interval, then
% the default interval $[-1,1]$ is used. For example, the following command
% makes a chebfun corresponding to $\cos(20x)$ on $[-1,1]$ and plots it.
  f = chebfun('cos(20*x)');
  plot(f), ylim([-1.2,1.2])

%%
% From this little experiment, you cannot see that |f| is represented by a
% polynomial.  One way to see this is to find the length of |f|:
  length(f)

%%
% Another is to remove the semicolon that suppresses output:
  f

%%
% These results tell us that |f| is represented by a polynomial interpolant
% through 51 Chebyshev points, i.e., a polynomial of degree 50.  These
% numbers have been determined by an adaptive process.  We can see the data
% points by plotting |f| with the |'.-'| option:
  plot(f,'.-'), ylim([-1.2 1.2])
  
%%
% The formula for $N+1$ Chebyshev points in $[-1,1]$ is
% $$ x(j) = -\cos(j \pi/N), \quad  j = 0:N, $$
% and in the figure we can see that the points are clustered accordingly
% near $1$ and $-1$. Note that in the middle of the grid, there are about 5
% points per wavelength, which is evidently what it takes to represent this
% cosine to 15 digits of accuracy.  For intervals other than $[-1,1]$,
% appropriate Chebyshev points are obtained by a linear scaling.

%%
% The curve between the data points is the polynomial interpolant, which 
% can be
% evaluated by the barycentric formula introduced by Salzer [Berrut &
% Trefethen 2004, Salzer 1972].  This method of evaluating polynomial
% interpolants is stable and efficient even if the degree is in the
% millions [Higham 2004].  In recent years Chebfun actually evaluates
% polynomials from their Chebyshev series rather than by barycentric
% interpolation; the difference in the two methods is little.

%%
% What is the integral of $f$ from $-1$ to $1$?  Here it is:
  sum(f)

%%
% This number was computed by integrating the polynomial (Clenshaw-Curtis
% quadrature -- see Section 2.1), and it is interesting to compare it
% to the exact answer from calculus:
  exact = sin(20)/10

%%
% Here is another example, now with the chebfun defined by an anonymous
% function instead of a string. In this case the interval is specified as
% $[0,100]$.
  g = chebfun(@(t) besselj(0,t),[0,100]);
  plot(g), ylim([-.5 1])

%%
% The function looks complicated, but it is actually a polynomial
% of surprisingly small degree:
  length(g)

%%
% Is it accurate?  Well, here are three random points
% in $[0,100]$:
  format long
  x = 100*rand(3,1)

%%
% Let's compare the chebfun to the true Bessel function at these points:
  exact = besselj(0,x);
  error = g(x) - exact;
  [g(x) exact error]

%%
% If you want to know the first 5 zeros of the Bessel function, here they
% are:
  r = roots(g); r = r(1:5)

%%
% Notice that we have just done something nontrivial and potentially
% useful.  How else would you find zeros of the Bessel function so readily?
% As always with numerical computation, we cannot expect the answers to be
% exactly correct, but they will usually be very close. In fact, these
% computed zeros are accurate to close to machine precision:
  besselj(0,r)

%%
% Most often we get a chebfun by operating on other chebfuns. For example,
% here is a sequence that uses plus, times, divide, and power operations on
% an initial chebfun |x| to produce a famous function of Runge:
  x = chebfun('x');
  f = 1./(1+25*x.^2);
  length(f)
  clf, plot(f)

%% 1.3  Operations on chebfuns
% There are more than 200 commands that can be applied to
% a chebfun.  For a list of many of them you can type |methods|:
  methods chebfun

%%
% To find out what a command does, you can use |help|.
  help chebfun/sum

%%
% Most of the commands in the list above exist in ordinary MATLAB; some
% exceptions are |domain|, |restrict|, |chebpoly|, and |remez|.
% We have already seen |length| and |sum| in action.  In fact we have
% already seen |subsref| too, since that is the MATLAB command for (among
% other things) evaluating arguments in parentheses.
%  Here is another example of its use:
  f(0.5)

%%
% Here for comparison is the true result:
  1/(1+25/4)

%% 
% In this Runge function example, we have also implicitly seen |times|,
% |plus|, |power|, and |rdivide|, all of which have been overloaded from
% their usual MATLAB uses to apply to chebfuns.

%%
% In the next part of this tour we shall explore many of these commands
% systematically.  First, however, we should see that chebfuns are not
% restricted to smooth functions.

%% 1.4  Piecewise smooth chebfuns
% Many functions of interest are not smooth but piecewise smooth.  In this
% case a chebfun may consist of a concatenation of smooth pieces, each with
% its own polynomial representation.  Each of the smooth pieces is called a
% "fun".  This enhancement of Chebfun was developed initially by Ricardo
% Pachon during 2006-2007, then also by Rodrigo Platte starting in 2007
% [Pachon, Platte and Trefethen 2010]. Essentially funs are the "classic
% chebfuns" for smooth functions on $[-1,1]$ originally implemented by
% Zachary Battles in Chebfun Version 1.

%%
% Later we shall describe the options in greater detail, but for the moment
% let us see some examples.  One way to get a piecewise smooth function is
% directly from the Chebfun constructor, taking advantage of its capability of
% automatic edge detection.  For example, in the default "splitting off"
% mode a function with a jump in its derivative produces a warning message,
  f = chebfun('abs(x-.3)');

%%
% The same function can be successfully captured with splitting on:
  f = chebfun('abs(x-.3)','splitting','on');

%%
% The |length| command reveals that |f| is defined by four data points,
% two for each linear interval:
  length(f)

%%
% We can see the structure of |f| in more detail by typing |f| without a
% semicolon:
  f

%%
% This output confirms that f consists of two funs, each defined by two
% points and two corresponding function values.
% The functions live on intervals defined by breakpoints at
% $-1$, $1$, and a number very close to $0.3$. 

%%
% Another way to make a piecewise smooth chebfun is to construct it
% explicitly from various pieces.  For example, the following command
% specifies three functions $x^2$, $1$, and $4-x$, together with a vector of
% endpoints indicating that the first function applies on $[-1,1]$, the
% second on $[1,2]$, and the third on $[2,4]$:
  f = chebfun({@(x) x.^2, 1, @(x) 4-x},[-1 1 2 4]);
  plot(f)

%%
% We expect |f| to consist of three pieces of lengths 3, 1, and 2, and this
% is indeed the case:
  f

%%
% Our eyes see pieces, but to Chebfun, |f| is just another function.  For
% example, here is its integral.
  sum(f)

%%
% Here is an algebraic transformation of |f|, which we plot in another color
% for variety.
  plot(1./(1+f),'r')

%% 
% Some Chebfun commands naturally introduce breakpoints in a chebfun. For
% example, the |abs| command first finds zeros of a function and introduces
% breakpoints there.  Here is a chebfun consisting of 6 funs:
  f = abs(exp(x).*sin(8*x));
  plot(f)

%%
% And here is an example where breakpoints are introduced by the |max|
% command, leading to a chebfun with 13 pieces:
  f = sin(20*x);
  g = exp(x-1);
  h = max(f,g);
  plot(h), ylim([0 1.2])

%%
% As always, |h| may look complicated to a human, but to Chebfun it is just a
% function.  Here are its mean, standard deviation, minimum, and maximum:
  mean(h)

%%
  std(h)

%%
  min(h)

%%  
  max(h)

%%
% A final note about piecewise smooth chebfuns is that the automatic edge
% detection or "splitting" feature, when it is turned on, may subdivide
% functions even though they do not have clean point singularities, and
% this may be desirable or undesirable depending on the application.  For
% example, considering $\sin(x)$ over $[0,1000\pi]$ with splitting on, we end up
% with a chebfun with many pieces:
  tic, f = chebfun('sin(x)',[0 1000*pi],'splitting','on'), toc

%%
% In this case it is more efficient -- and more interesting mathematically
% -- to omit the splitting and construct one global chebfun:

  tic, f2 = chebfun('sin(x)',[0 1000*pi]), toc

%%
% Splitting on and off are discussed further in Section 8.3.

%% 1.5  Infinite intervals and infinite function values
% A major change from Chebfun Version 2 to Version 3 was the generalization of
% chebfuns to allow certain functions on infinite intervals or which
% diverge to infinity; the initial credit for these innovations belongs to
% Nick Hale, Rodrigo Platte, and Mark Richardson, and many later
% developments are due to Kuan Xu.
% For example, here is a function on the whole real axis,
  f = chebfun('exp(-x.^2/16).*(1+.2*cos(10*x))',[-inf,inf]);
  plot(f)

%%
% and here is its integral:
  sum(f)

%%
% Here's the integral of a function on $[1,\infty)$:
  sum(chebfun('1./x.^4',[1 inf]))

%%
% Notice that several digits of accuracy have been lost here.  Be careful! --
% operations involving infinities in Chebfun are not always as accurate
% and robust as their finite counterparts.

%%
% Here is an example of a function that diverges to infinity,
% which we can capture with the |'exps'| flag; see Chapter 7
% for details:
  h = chebfun('(1/pi)./sqrt(1-x.^2)','exps',[-.5 -.5]);
  plot(h)

%%
% In this case the integral comes out just right:
  sum(h)

%%
% For more on the treatment of infinities in Chebfun, see Chapter 9.

%% 1.6  Periodic functions
% Until 2014, Chebfun used only nonperiodic representations, based on
% Chebyshev polynomials.  Beginning with Version 5, there is a new capability
% of representing sufficiently smooth periodic functions by trigonometric
% polynomials instead, that is, Fourier series.  Such an object is still
% called a chebfun, but it is a periodic one,
% and the signal to invoke such capabilities is the string |`trig`|.
% For abbreviation, we call a periodic chebfun a ``trigfun''.
% This section gives a quick introduction, and more details can be
% found in Chapter 11.

%%
% Trigfuns were initiated by Grady Wright in the first half of 2014
% [Wright et al. 2015].
% A very interesting project along the same lines was
% carried out independently by Kristyn McLeod and Rodrigo Platte at Arizona
% State University [McLeod 2014].

%%
% For example, here is a periodic function on $[-\pi,\pi]$ represented
% in the usual way by a Chebyshev series.
ff = @(t) sin(t) + cos(2*t) - cos(t)/3 + cos(100*t)/6;
f = chebfun(ff,[-pi,pi]);
max(f)
plot(f)

%%
% Its length, very roughly, is $100 \pi$, 
length(f)

%%
% Here is the same function represented by a Fourier series:
f2 = chebfun(ff,[-pi,pi],'trig')
max(f2)
plot(f2,'m')

%%
% Its length is now only about $200$ (exactly 201). This improvement
% by a factor of about $\pi/2$ is typical.
length(f2)

%%
% Sampling at a few arbitrary points confirms that the
% two functions agree closely:
xx = [1/3 sqrt(2) exp(1)];
f(xx) - f2(xx)

%%
% Readers may be interested to compare |plotcoeffs| applied
% to the first and second versions of $f$.  Rather than display
% that here we shall turn to a simpler example involving a 
% shorter Fourier series.  Consider the function
f = chebfun('7 + sin(t) + exp(1)*cos(2*t)',[-pi,pi],'trig')

%%
% Here are the coefficients of $f$ as an expansion in sines and
% cosines:
[a,b] = trigcoeffs(f)

%%
% Here they are as an expansion in complex
% exponentials:
c = trigcoeffs(f)

%%
% Bookkeeping of Fourier coefficients can often be a headache. If
% these examples don't make the patterns clear, details can be found with
% |help trigcoeffs|.

%%
% For a mathematically less trivial example,
% here is the cosine expansion of a function whose Fourier series coefficients
% are known to be values of a Bessel function:
f = chebfun('exp(cos(t))',[-pi pi],'trig');
[a,b] = trigcoeffs(f);
n = floor(length(f)/2);
exact = 2*besseli(0:n,1); exact(1) = exact(1)/2;
disp('        computed             exact')
disp([a exact'])

%% 1.7  Rows, columns, and quasimatrices
% MATLAB doesn't only deal with column vectors: there are also row vectors
% and matrices.  The same is true of Chebfun. The chebfuns shown so far
% have all been in column orientation, which is the default, but one can
% also take the transpose, compute inner products, and so on:

%%
  x = chebfun(@(x) x)
%%
  x'
%%
  x'*x
%%
% One can also make matrices whose columns are chebfuns
% or whose rows are chebfuns, like this:
  A = [1 x x.^2]

%%
  A'*A

%%
% These are called _quasimatrices_, and they are discussed in
% Chapter 6.

%% 1.8  Chebfun features not in this Guide
% Some of Chebfun's most remarkable features haven't made it
% into this edition of the Guide.  Here are some of our favorites:

%%
% o |leg2cheb| and |cheb2leg| for fast Legendre-Chebyshev conversions (also
% |legvals2chebcoeffs|, |chebcoeffs2legpts|, and ten more)

%%
% o |conv| and |circconv| for convolution,

%%
% o The |'equi'| flag to the Chebfun constructor for equispaced data,

%%
% o |polyfit| for least-squares fitting in the continuous context,

%%
% o |inv| for computing the inverse of a chebfun,

%%
% o |pde15s| for PDEs in one space and one time variable.

%%
% To learn about any of these options, try the appropriate |help|
% command.  Just as a taster, here's a hint of how fast Chebfun
% can convert a ten-thousand coefficient Chebyshev expansion
% to Legendre coefficients and back again using an
% algorithm from [Hale & Townsend 2013]:
tic
ccheb = randn(10000,1);
cleg = cheb2leg(ccheb);
ccheb2 = leg2cheb(cleg);
norm(ccheb-ccheb2,inf)
toc

%% 1.9  Chebfun example galleries
% MATLAB has long had a `gallery` command to generate interesting
% matrices, and with version 5.1, Chebfun has introduced an analogous
% `gallery` command to generate interesting functions.

%%
% Here is what is currently available:
help cheb.gallery

%%
% For example, here is a chebfun representing the Airy function,
plot(cheb.gallery('airy')), ylim([-.8 .8])
title('Airy function')

%%
% In this instance the underlying code fits in a line,
fa = @airy; p = chebfun(fa, [-40,40]); 

%%
% Some examples make use of more complicated code, like this 
% approximation to a Daubechies wavelet scaling function
% (accurate to about 3 digits of accuracy; the underlying
% function is a fractal):
plot(cheb.gallery('daubechies')), ylim([-0.5 1.5])
title('Daubechies scaling function')

%%
% To find out how a gallery example was generated, take a look
% at the code with `type +cheb/gallery` or `edit +cheb/gallery`.

%%
% Like the MATLAB `gallery` command, `cheb.gallery` produces a
% plot if you call it without specifying output variables.
% To illustrate, let us finish with an example the
% Chebfun team enjoys from the appendix to
% [Trefethen 2013], "Six myths of polynomial interpolation and
% quadrature":
cheb.gallery('zigzag')

%%
% This function looks piecewise linear, but in fact, it is a 
% polynomial of degree 10000.  This serves no purpose from
% an approximation point of view -- one would never represent
% this function in this manner -- but it illustrates the robustness of
% high-degree polynomial approximation.

%%
% If you call `cheb.gallery` without any input arguments, it selects
% a gallery function at random.

%%
% Other collections worth exploring are `cheb.gallerytrig` for
% periodic functions and `cheb.gallery2` for 2D functions.

%% 1.10  How this Guide is produced
% This guide is produced in MATLAB using the |publish| command with a 
% style sheet somewhat different from the usual; the output of |publish|
% is then processed by Markdown.
% To publish a chapter for yourself, make sure the
% chebfun guide directory is in your path and then type, for example,
% |open(publish('guide1'))|.  The formatting may not be exactly 
% right but it should certainly be intelligible.

%% 1.11  References
% [Battles & Trefethen 2004] Z. Battles and L. N. Trefethen, "An extension
% of MATLAB to continuous functions and operators", _SIAM Journal on
% Scientific Computing_, 25 (2004), 1743-1770.
% 
% [Berrut & Trefethen 2005] J.-P. Berrut and L. N. Trefethen, "Barycentric
% Lagrange interpolation", _SIAM Review 46_, (2004), 501-517.
%
% [Hale & Townsend 2013]  N. Hale and A. Townsend,
% A fast, simple, and stable Chebyshev--Legendre
% transform using an asymptotic formula, _SIAM Journal on Scientific
% Computing_, 36 (2014), A148-A167.
%
% [Higham 2004] N. J. Higham, "The numerical stability of barycentric
% Lagrange interpolation", _IMA Journal of Numerical Analysis_, 24 (2004),
% 547-556.
%
% [McLeod 2014] K. N. McLeod, "Fourfun: A new system for automatic
% computations using Fourier expansions," _SIAM Undergraduate Research
% Online_, 7 (2014), `http://dx.doi.org/10.1137/14S013238`.
%
% [Pachon, Platte & Trefethen 2010] R. Pachon, R. B. Platte and L. N.
% Trefethen, "Piecewise-smooth chebfuns", _IMA J. Numer. Anal._, 30 (2010),
% 898-916.
%
% [Salzer 1972] H. E. Salzer, "Lagrangian interpolation at the Chebyshev
% points cos(nu pi/n), nu = 0(1)n; some unnoted advantages", _Computer
% Journal_ 15 (1972), 156-159.
%
% [Trefethen 2007] L. N. Trefethen, "Computing numerically with functions
% instead of numbers", _Mathematics in Computer Science_ 1 (2007), 9-19.
% Revised and reprinted in _Communications of the ACM_ 58 (2014),
% 91-97.
%
% [Trefethen 2013] L. N. Trefethen, _Approximation Theory and
% Approximation Practice_, SIAM, 2013.
% 
% [Wright et al. 2015] G. B. Wright, M. Javed, H. Montanelli, and
% L. N. Trefethen, Extension of Chebfun to periodic functions,
% _SIAM J. Sci. Comp._, to appear.
