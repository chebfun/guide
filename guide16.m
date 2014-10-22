%% 16. Periodic Chebfuns
% Grady B. Wright, October 2014

%%
LW = 'linewidth'; lw = 1.6; MS = 'MarkerSize'; ms = 10;

%%
% One of the major new features of Chebfun version 5 is the ability to 
% use trigonometric interpolants instead of Chebyshev interpolants for 
% representing smooth periodic functions.  These trig-based chebfuns, or
% "trigfuns" as we like to refer to them, can be created with the use of
% the `'periodic'` (or `'trig'`) flag in the chebfun constructor.  For
% example, the periodic function $f(x) = \cos(8\sin(x+1/7))$ on
% $[-\pi,\pi]$ can be constructed and plotted as follows:
dom = [-pi,pi];
f = chebfun(@(x) cos(8*sin(x+1/7)),dom,'trig')
plot(f,LW,lw);

%%
% The text `'trig'` on the right in the displayed information for $f$
% above indicates that it is represented using a trigonometric interpolant,
% which we discuss in more detail in the next section.  For now we note that the
% the length of 61 for $f$ displayed above means that it is resolved to
% machine precision using a trigonometric interpolant through 61 equally
% spaced samples of $f$ over $[-\pi,\pi)$, or equivalently, using 30
% trigonometric (or Fourier) modes.

%%
% In this chapter we review some of the functionality available for
% trigfuns as well as some theory of trigonometric interpolation.  We start
% with a discusion of the latter first.  To simplify some of the
% terminology below, we will refer to trigonometric based chebfuns as
% trigfuns and Chebyshev polynomial based chebfuns as
% simply chebfuns.

%% 16.1  Trigonometric series and interpolants
% The classical trigonometric series of a periodic function
% $u$ defined on $[-\pi,\pi]$ is given formally as
%
% $$ \mathcal{F}[u] = \sum_{k=-\infty}^{\infty} a_k e^{ikx} $$
%
% where the coefficients are given as
%
% $$ a_k = \frac{1}{2\pi} \int_{-\pi}^{\pi} u(x)e^{-ikx} dx. $$
%
% Alternatively, we can express the series in terms of sines and cosines:
%
% $$ \mathcal{F}[u] = \sum_{k=0}^{\infty} a_k \cos(k x) +
% \sum_{k=1}^{\infty} b_k \sin(k x) $$
%
% where 
%
% $$ a_k = \frac{1}{\pi} \int_{-\pi}^{\pi} f(x)\cos(kx) dx, \quad
% b_k = \frac{1}{\pi} \int_{-\pi}^{\pi} f(x)\sin(kx) dx, $$
%
% for $k>0$ and
%
% $$ a_0 = \frac{1}{2\pi} \int_{-\pi}^{\pi} f(x) dx. $$
%
% However, in what follows we will use the complex exponential form the the series.
%
% Note that the series above is often refered to as the _Fourier series_ of
% $u$, but we will use the term _trigonometric series_ as advocated by
% Zygmund [Zygmund, 1959].  Also, note that similar expressions for the series and
% coefficients hold for more general intervals $[a,b]$ by shifting and
% scaling appropriately. In fact, the default interval for trigfuns is
% $[-1,1]$ to keep things consistent with the default for chebfuns.
%
% The approximation theory results for Trigonometric series are classic
% and can be found in, for example, [Zygmund, 1959].  Here we review some
% of these properties.
%
% The convergence of the trigonometric series of $u$ depends on the
% smoothness of $u$ on $[-\pi,\pi]$ and its periodic extension. Let $q_N$
% be the truncated series
%
% $$ q_N(x) = \sum_{k=-(N-1)/2}^{(N-1)/2} a_k e^{ikx}, $$
%
% when $N$ is odd and 
%
% $$ q_N(x) = \sum_{k=-N/2}^{N/2} a_k e^{ikx}, $$
%
% when $N$ is even. If $u$ is periodic, continuous, and of bounded variation on 
% $[-\pi,\pi]$ then as $N\rightarrow \infty$
% 
% $$\| u - q_N \| \rightarrow 0,$$
% 
% where $\|\cdot\|$ is the maximum norm. The following estimate shows that
% the difference between $u$ and $q_N$ is dominated by the tail of the 
% series:
%
% $$\| u - q_N \| \leq \sum_{|k| > \lfloor N/2 \rfloor} |a_k|. $$
%
% The decay rate of the coefficients $a_k$ also depend on the smoothness
% of $u$ and can readily be obtained using integration by parts.  A classical
% result says that if $u$ is $(\ell-1)$-times continuously differentiable
% on $[-\pi,\pi]$, with each of these derivatives being periodic, and 
% $\ell$th derivative is of bounded variation on $[-\pi,\pi]$ then 
% 
% $$ |a_k| = O(|k|^{-\ell}),\; k=\pm 1, \pm 2,\ldots $$
%
% This estimate can be used in combination with the bound above on 
% $\| u - q_N \|$ to obtain the error estimate:
%
% $$ \| u - q_N \| \leq C N^{-\ell},  $$
% 
% for some constant $C > 0$ independent of $N$. If $u$ and its periodic
% extension are $C^{\infty}$ then $q_N$ converges to $u$ faster than any
% polynomial power as $N\rightarrow\infty$.  If $u$ is _analytic_ on
% $[-\pi,\pi]$ then the convergence rate is _exponential_, i.e. 
% $\| u - q_N \| = O\left(b^N\right)$ for $0 < b < 1$.
%
% Computing the coefficients $a_k$ exactly is, of course, impossible
% in practice, so a quadrature formula is instead applied to obtain an
% approximate solution. Letting $x_j = -\pi + 2\pi j/N$, for $j=0,\ldots,N-1$, 
% and applying the Trapezoidal rule gives the following approximation:
%
% $$ a_k \approx c_k := \frac{1}{N} \sum_{j=0}^{N-1} u(x_j) e^{-i k x_j}. $$
% 
% The coefficients $c_k$ are typically referred to as the _Discrete Fourier Coefficients_ and when
% they are replaced by $a_k$ in the truncated trigonometric series $q_N$, the resulting
% series is called the _discrete Fourier series_ of $u$.  The form of this series 
% depends on the parity of $N$.  For $N$ odd we have
%
% $$ p_N(x) = \sum_{k=-(N-1)/2}^{(N-1)/2} c_k e^{ikx}. $$
%
% But, when $N$ is even we have 
%
% $$ p_N(x) = \sum_{k=-N/2+1}^{N/2-1} c_k e^{ikx}\,+ c_{N/2} \cos(N/2 x), $$
%
% which is often rewritten as
%
% $$ p_N(x) = {\sum_{k=-N/2}^{N/2}} \phantom{a}^{\prime} c_k e^{ikx}, $$
% 
% where $c_{-N/2} = c_{N/2}$ and the prime means the first and
% last terms in the sum are halved.  The reason the $\sin(N/2 x)$ term is
% missing is that this function vanishes when evaluated at $x_j$, so that
% there is no contribution of this mode in the coefficient $c_{N/2}$.

%%
% The discrete Fourier series $p_N$ of a function $u$ has the interesting property that it
% _interpolates_ $u$ at the equally spaced grid points $x_j = -\pi+2\pi
% j/N$, for $j=0,\ldots,N-1$. It is one example of an interpolant whose
% coefficients can be computed by a quadrature rule (the Trapezoidal rule
% in this case). 

%%
% A note on terminology: We prefer to refer to $p_N$ as the _trignometric 
% interpolant_ of $u$ to make the connection to interpolation crystal clear.  
% Trigonometric interpolation is the underlying technology used for periodic
% chebfuns (trigfuns) just as polynomial interpolants at Chebyshev nodes are the
% underlying technology used for non-periodic chebfuns (chebfuns).

%%
% Besides simplicity, there are two (maybe three) main reasons that trigonometric 
% interpolants at equally spaced points with coefficients determined by the
% trapezoidal rule are immensely useful:
%
% 1. The coefficients $c_k$ can be computed in $O(N\log N)$ operations
% using the fast Fourier transform [Van Loan 1992] and the computation is
% numerically stable [Henrici 1986].
%
% 2. [Perhaps we should say something about the series being able to be
% evaluate stably using Horner's scheme?].
%
% 3. The approximation properties of the interpolants are very similar to
% those of the truncated trigonometric series approximation of a function
% $u$; see, for example, [Canuto et al. 2006/7], [Hesthaven et. al. 2007], 
% or [Trefethen 2000] for a more indepth discussion.
%
% The latter result follows from the Poisson summation formula 
% (or aliasing formula), which relates the interpolation coefficients 
% (discrete Fourier coefficients) $c_k$ to the trigonometric series
% coefficients (continuous Fourier coefficients) $a_k$ as follows:
%
% $$ c_k = a_k + \sum_{\ell=1}^{\infty} \left(a_{k+\ell N} + a_{k-\ell N}\right).$$
%
% This formula shows that the decay rate of $c_k$ is the same as that 
% of $a_k$, which we know depend on the smoothness of $u$.  Furthermore, it
% follows that if $u$ is sufficiently smooth (e.g. continuous with one
% derivative of bounded variation), each $c_k$ converges to $a_k$ as
% $N\to\infty$ (up to rounding errors). 

%%
% The Poisson summation formula shows that the same estimates as quoted above
% for the truncated trigonometric series, also apply to the trigonometric
% interpolant, only the constants will be different.  In fact, one can show
% that the trigonometric interpolant can at most differ from the 
% truncated series approximation with the same number of terms by a factor
% of at most two (provided $u$ is continuous with one derivative of bounded
% variation).  We can illustrate this nicely in Chebfun for the function
% $u(x) = |x|^3$ over $[-\pi,\pi]$ as follows:

%%
% First we construct the periodic version of $u$ approximated to machine
% precision:
u_exact = @(x) abs(sin(x)).^3;
u = chebfun(u_exact,[-pi,pi],'trig');

%%
% Next we construct the truncated trigonometric series approximation with
% $N=11$ terms using the "`trunc'" option:
q11 = chebfun(u_exact,[-pi,pi],'trunc',11,'trig');

%%
% We follow this by constructing the $N=11$ point trigonometric interpolant
% by overriding the default adaptive construction:
p11 = chebfun(u_exact,[-pi,pi],11,'trig');

%%
% The differences between the true function and the two approximations are
% illustrated on the plot below.
plot(q11-u,'b-',p11-u,'r-',LW,lw), legend('u-q_{11}','u-p_{11}')

%%
% The ratio of the max-norms of these differences is given as 
norm(p11-u,inf)/norm(q11-u,inf)

%%
% and is clearly within the bounds predicted by the theory.

%%
% In general, if $N$ is sufficiently large and $u$ is sufficiently smooth,
% then the there will be no numerical differences between the coefficients
% $a_k$ and $c_k$, and thus also no differences between the trigonometric
% interpolant and truncated series. The Chebfun system is meant to operate
% in this regime.

%% 16.2 Basic operations
% Many of the basic opertions for computing with functions in Chebfun 
% can also be applied directly to trigfuns with only one small change: the
% initial construction of the function using the `'periodic'` flag.

%%
% The operations of addition, subtraction, multiplication, division, and
% function composition can all be directly applied to a trigfun.  However,
% one should be aware that operations performed should result in a smooth and periodic
% function. For example the following illustrates some of these operations:
g = chebfun(@(x) sin(x),dom,'trig');
f = tanh(cos(1+2*g).^2)-0.5
plot(f, LW, lw)

%%
% The max, min, and roots of a trigfun can be computed just as in Chebfun.
% For example, for $f$ defined above we have
[maxf,xmaxf] = max(f);
[minf,xminf] = min(f);
rootsf = roots(f);
maxf
minf
rootsf

%%
% These can be visualized as
plot(f, LW, lw), hold on
plot(xmaxf,maxf,'gs',xminf,minf,'md',rootsf,0*rootsf,'ro',MS,ms)
legend('f','max f','min f','zeros f','location','southwest')
hold off;

%%
% The derivative of trigfun is computed using `diff`:
df = diff(f);
plot(df, LW, lw)

%%
% Similarly, the definite integral is computed using `sum`:
intf = sum(f)

%%
% Comparing this latter approximation to the standard Chebfun result, we
% see the approximations are identical:
g = chebfun(@(x) sin(x),dom);
f = tanh(cos(1+2*g).^2)-0.5;
sum(f)

%% 16.3 Complex-valued functions
% Complex-valued trigfuns are also possible and can be quite useful for
% computing certain integrals over closed contours in the complex-plane.
% To see why, consider a smooth and closed contour $\Gamma$ in the complex plane and
% suppose we wish to compute
%
% $$ \int_{\Gamma} f(z) dz. $$
%
% If we parametrize $\Gamma$ using a real varaible, say $t$, then since the
% contour is closed, the integrand becomes periodic in $t$ and we get
%
% $$ \int_{\Gamma} f(z) dz = \int_{a}^{b} f(z(t)) z'(t) dt. $$


%%
% Assuming $f$ and the parameterization of $\Gamma$ are smooth, the
% integrand can then be represented quite efficiently using a trigfun.

%%
% Here is a simple example for the function:
ff = @(z) (1-2*z)./(z.*(z-1).*(z-3));

%%
% Suppose we want to integrate this function along heart-shaped contour
% given by the parameterization
%
% $$ z(t) = 2\sin^3(t) + i/8(13\cos(t)-5\cos(2t)-2\cos(3t)-\cos(4t)), $$
%
% where $0 \leq t \leq 2\pi$.

%%
% To do this we first construct a trigfun representing the contour:
zz = @(t) 2*sin(t).^3 + 1i/8*(13*cos(t)-5*cos(2*t)-2*cos(3*t)-cos(4*t));
z = chebfun(zz,[0,2*pi],'periodic');
plot(z,LW,lw), axis equal, title('Contour of integration')

%%
% The integrand is then constructed by a simple composition:
f = ff(z);

%%
% This is how the real and imaginary parts of the integrand
% look on the contour:
subplot(1, 2, 1)
plot(real(f), LW, lw)
title('real part')
subplot(1, 2, 2)
plot(imag(f), LW, lw)
title('imaginary part')

%%
% To compute the integral, we recall that
%
% $$ \int_{\Gamma} f(z) dz = \int_{0}^{2\pi} f(z(t)) z'(t) dt. $$
%
% We therefore first compute $z'(t)$:
dz = diff(z);
%%
% The contour integral is then computed simply as
s = sum(f.*dz)

%%
% Comparing this to the true answer of $-5 \pi i/3$ (computed using
% residues), we see that the Chebfun result is very good:
norm(-5/3*pi*1i-s)

%%
% Here is one more example involving a function with an essential singularity at the
% origin and a contour integral on the unit circle.
ff = @(z) exp(1./z).*sin(1./z);
z = chebfun(@(t) exp(1i*t),[0 2*pi],'trig');
f = ff(z);
dz = diff(z);
s = sum(f.*dz)

%%
% The result nicely matches the exact result of $2\pi i$:
exact = 2i*pi

%% 16.4 Circular convolution
% The circular (or periodic) convolution of two smooth periodic functions
% $f(x)$ and $g(x)$, with period $T$ is defined as
%
% $$ (f*g)(x) := \int_{x_0}^{x_0 + T} g(\xi)f(x-\xi)d\xi $$
%
% where $x_0$ is any aribtrary point. These types of convolutions can be
% computed for 2 trigfuns using the `circconv` function. 

%%
% The example below 
% demonstrates this function in combination with the additional feature 
% that allows fourfuns to be constructed from function values. The latter
% is demonstrated first:
rng('default'), rng(0);
n = 201;
x = trigpts(n);
func_vals = exp(sin(2*pi*x)) + 0.05*randn(n,1);
f = chebfun(func_vals,dom,'trig')

%%
% Here $f$ interpolates the noisy `func_vals` at 201 equally spaced points
% from $[-\pi,\pi)$ using trigonometric polynomials. The high frequencies in this
% function can be smoothed by convolving it with a mollifier, in this case
% we use a (normalized) Gaussian with variance 0.1.
sigma = 0.1;
g = chebfun(@(x) 1/(sigma*sqrt(2*pi))*exp(-0.5*(x/sigma).^2),dom,'trig');

%%
% Note that the resulting respresentation of $g$ is actually the periodic 
% extension of the Gaussian over $[-\pi,\pi]$.  Since $g$ decays quickly
% over the interval $[-\pi,\pi]$, it can be represented to machine precision
% over this interval using the trigonometric interpolant.

%% 
% The circular convolution of $f$ and $g$ is computed and visualized as follows:
h = circconv(f,g);
clf
plot(g,'b',f,'r',h,'k',LW,lw);
legend('Molifier g','Noisy function f','Smoothed function h');
hold off;

%%
% Note that `circconv` has different functionality than the chebfun 
% `conv` function which applies to non-periodic functions defined over
% potentially different intervals.

%% 16.5 trigfuns vs. chebfuns
% Trigonometric interpolants have a resolution power of 2 points per 
% wavelength whereas Chebyshev interpolants require approximately $\pi$
% points per wavelength.  This means that smooth periodic functions can 
% be represented as trigfuns using less information (i.e. fewer samples)
% than standard chebfuns.

%%
% To illustrate this consider the chebfun and trigfun representations of
% the function $f(x) = \cos(11\sin(x-1/\pi))$ over $[-\pi,\pi]$:
ff = @(x) cos(11*sin(3*(x-1/pi)));
f_cheb = chebfun(ff,dom);
f_trig = chebfun(ff,dom,'trig');

%%
% The ratio of length of the two representations should be approximately
% $\pi/2$ and this is indeed what we find:
ratio = length(f_cheb)/length(f_trig)
pi/2

%%
% While this compression in the representation of a smooth periodic
% function is nice, it's not as important as the benefits gained in
% computing (high) derivatives of smooth periodic functions represented
% by trigfuns.  For example, consider the function
f = chebfun(@(x)cos(10*sin(x)),dom,'trig');
plot(f,LW,lw);

%%
% All odd derivatives of $f$ vanish at $\pm \pi$.  We can check how good the
% trigfun approximates the 9th derivative at $\pi$ is using
df9 = diff(f,9);
df9(pi)/norm(df9,inf)

%%
% Here we have normalized by the max-norm of the 9th derivative since the
% derivatives are growing in magnitude.  We can compare this to the 
% standard chebfun representation as
f_cheb = chebfun(@(x)cos(10*sin(x)),dom);
df9_cheb = diff(f_cheb,9);
df9_cheb(pi)/norm(df9,inf)

%%
% The trigfun result clearly provides a better approximation to the high
% derivatives of a smooth periodic function.

%%
% Trying to construct a trigfun from a non-periodic or non-smooth function
% will typically result in a warning being issued and an "unhappy" trigfun,
% as illustrated for the unit step function below:
f = chebfun(@(x) 0.5*(1+sign(x)),dom,'trig')
plot(f,LW,lw);

%%
% The length of $f$ is 65536, which is the maximum number of samples used
% in the construction process to try to resolve $f$. The famous Gibbs'
% phenomenon can be seen near the discontinuity in the plot of $f$. Of course,
% Chebfun can be used to represent this function in non-periodic mode (i.e. using
% Chebyshev interpolants) with the option of `splitting on`:
f = chebfun(@(x) 0.5*(1+sign(x)),dom,'splitting','on')
plot(f,LW,lw);

%%
% Splitting is not an option for trigfuns. In the next section we return to
% the step function example and show how one can also build up
% trigonometric series approximations using a specified number of terms.

%%
% The standard operations of addition, subtraction, multiplication, and
% division are currently supported between a trigfun and a chebfun. The
% result of the operation, however, is always a chebfun since these can be
% used to represent non-periodic, and non-smooth functions.  To illustrate
% this consider the periodic function $f(x) = \tanh(\cos(3(x-1/7))$ and 
% the non-periodic unit step function $g(x)$:
f = chebfun(@(x) tanh(cos(3*(x-1/7))),dom,'trig');
g = chebfun(@(x) 0.5*(1+sign(x)),dom,'splitting','on');
plot(f,'b',g,'r',LW,lw)

%%
% The product of $f$ and $g$ is given as
h = f.*g;
plot(h)

%%
% Note that $h$ is now represented as a chebfun.

%% 16.6 Trigonometric coefficients
% Another useful operation provided by the trigonometric technology is the 
% capability of computing the coefficients $a_k$ in the trigonometric series
% approximation $\mathcal{F}[u]$ discussed at the beginning of Section
% 16.1.  The function `trigcoeffs` provides this functionality.

%%
% As discussed in Section 16.1, the decay rate of the trigonometric
% coefficients $a_k$ is directly related to the decay of the underlying function.
% Typically, if $u$ and its periodic extension are twice continuously
% differentiable over $[-\pi,\pi]$ and $u'''(x)$ is piecewise continuous
% on $[-\pi,\pi]$ (or more specifically of bounded variation) then the
% trignometric coefficients of $u$ can be quickly computed by first constructing
% a trigfun representation of $u$, then calling the `trigcoeffs` function.
% Here is an example for a simple trigonometric polynomial:
u = chebfun(@(x) 1 - 4*cos(x) + 6*sin(2*x),dom,'trig');
a = trigcoeffs(u);
disp('Trig coeffs of 1 - 4*cos(x) + 6*sin(2*x):')
a

%%
% Note that `trigcoeffs` returns the coefficients from lowest degree
% to highest degree just as they appear in the definition for $q_N$ given
% in Section 16.1.  Also note that `trigcoeffs` by default returns the 
% coefficients in complex exponential form.  The equivalent coefficients in
% terms of cosines and sines can be obtained as:
[a,b] = trigcoeffs(u);
disp('Trig cos coeffs of 1 - 4*cos(x) + 6*sin(2*x)')
a
disp('Trig sin coeffs of 1 - 4*cos(x) + 6*sin(2*x)')
b

%%
% Note that `a` contains the constant term in the series as its first
% value followed by values of the coefficients for $\cos(x)$ and $\cos(2x)$, 
% while `b` starts with the coefficient for $\sin(x)$ followed by the 
% value of the coefficient for $\sin(2x)$.

%%
% The default behavior of `trigcoeffs` is to return all the trigonometric
% coefficients necessary to resolve the function to machine precision
% (assuming this number is less than 65537).  However, a specific number
% can be obtained with an additional input argument.  We illustrate this
% feature on the function $f(x) = 3/(5 - 4\cos(x))$, which is analytic in a 
% strip in the complex plane and has exact trigonometirc coefficients 
% (in complex exponential form) given by $a_k = 2^{-|k|}$:
numCoeffs = 11;
u = chebfun(@(x) 3./(5 - 4*cos(x)),dom,'trig');
a = trigcoeffs(u,numCoeffs);
disp('Trig coeffs of 3/(5-4cos(x)):')
a

%%
% We see that the computed results match the exact results to machine
% precision. All the coefficients can be visualized using the `plotcoeffs`
% function
plotcoeffs(u,'x-')

%%
% Here is an example for a less smooth function where we compute only the
% first 17 coefficients:
numCoeffs = 17;
u = chebfun(@(x) abs(sin(x)).^3,[-pi,pi],'trig');
a = trigcoeffs(u,numCoeffs);
disp('Trig coeffs of |sin(x)|^3')
a

%%
% We see that the coefficients decay much more slowly in this case
% in fact the number of terms required to resolve this function to machine
% precision is:
length(u)

%%
% The reason is that this function has only two continuous derivatives in
% $L^{2}[-\pi,\pi]$ and a piecewise continuous third derivative, so its
% Trigonometric coefficients decay as $O(|k|^{-4})$, as discussed in
% Section 16.1.  This decay rate can be seen by plotting the
% coefficients on a log-log scale, which can be easily done for the
% positive mode coefficients (i.e. $k>0$) using the `plotcoeffs` command:
plotcoeffs(u,'loglog'), ylim([1e-15 1]), hold on
k = [100 length(u)/2];
plot(k,10*k.^(-4),'k-',LW,lw)
text(500,50*(500)^(-4),'O(k^{-4})','FontSize',12)
hold off

%%
% The trigonometric coefficients of functions (and their periodic extensions) 
% with fewer than two continuous derivatives can also be computed.  However,
% the functions must first be constructed as a chebfun, i.e. using the default, 'non-periodic',
% option.  In this case the trigonometric coefficients are computed using the
% integral formulas defined at the beginning of section 16.1 (with the actual
% computations done with Chebfun's `sum` method) instead of the fast 
% Fourier Transform.
%
% The quintessential example of a non-smooth function is that of the square
% wave (or periodic extension of the step function), which can be defined
% by the `sign` function as
%
% $$ u(x) = \mbox{sign}(\sin(x)) $$ 
%
% This can be constructed as done in the previous section using Chebfun
% with 'splitting' turned on,
sq_wave = @(x) sign(sin((x)));
u = chebfun(sq_wave,dom,'splitting','on');

%%
% We can obtain the trigonometric coefficients of this function using again the
% `trigcoeffs` method.  Since the square wave is odd, it makes sense to
% only look at the sine coefficients in this case:
numCoeffs = 15;
[a,b] = trigcoeffs(u,numCoeffs);
disp('Fourier sine coeffs of unit step function:')
b

%%
% The exact values of the coefficients are 
%
% $$ b_k = \frac{4}{k\pi} ~~(k \mbox{ odd}) $$
%
% with $b_k = 0$ for $k$ even, $k\ge 1$.
% These values can be easily seen in the computed results:
disp('            k               pi/4*b_k')
disp([(1:7)' pi/4*real(b)])

%%
% The computed cosine coefficients are numerically zero, as expected:
norm(a,inf)

%%  16.7 Truncated trigonometric series approximations
% The `trigcoeffs` function can also be used in combination with the
% Chebfun constructor to generate truncated trigonometric series representations
% of functions of non-smooth functions.  Here's an example for the above square wave using 15
% trigonometric modes:
numModes = 15;
a = trigcoeffs(u,2*numModes+1);
u_trunc = chebfun(a,dom,'trig','coeffs');
plot(u,'k-',u_trunc,'b-',LW,lw)

%%
% This represents the best 15-mode trigonometric approximation to the square
% wave over $[-\pi,\pi]$ in the $L^2$ sense. The oscillations in the
% approximation are called the Gibbs' phenomenon.

%%
% To see the actual 'wave' it is useful to plot the approximation over 
% a larger interval, which can be done for $-4\pi \leq x \leq 4\pi$ as follows:
u = chebfun(sq_wave,[-4*pi,4*pi],'splitting','on');
u_trunc = chebfun(u_trunc,[-4*pi,4*pi],'trig');
plot(u,'k-',u_trunc,'b-',LW,lw)

%%
% Note that in this case we don't recompute the trigonometric coefficients over
% the larger domain, we simply construct `u_trunc` over the larger domain.

%%
% We conclude with one more classic example: the sawtooth wave.
% Again, we compute this over $[-\pi,\pi]$
% then expanded to a larger domain:
sawtooth = @(x) (mod(x+pi,2*pi))/(2*pi);
u = chebfun(sawtooth,dom,'splitting','on');
a = trigcoeffs(u,2*numModes+1);
u_trunc = chebfun(a,[-pi,pi],'trig','coeffs');

u = chebfun(sawtooth,[-4*pi,4*pi],'splitting','on');
u_trunc = chebfun(u_trunc,[-4*pi,4*pi],'trig');
plot(u,'k-',u_trunc,'b-',LW,lw)

%%
% To hear this wave using try `chebtune(u_trunc,6)`.

%% 16.8 Concluding remarks
% It should be noted that another MATLAB package based on trigonometric 
% interpolants called fourfun was developed independently
% by Kristyn Mcleod under the supervision of Rodrigo Platte [Mcleod 2014].
% This standalone MATLAB package has many common features with 
% the trigfun part of Chebfun, but does not have all the functionality that
% comes with the integrated Chebfun environment.
%
% [What's missing: trigbary, trigcf, and ODEs.  These perhaps could be separate
% guide chapters.]

%% 16.9  References
%
% [Canuto et al. 2006/7] C. Canuto, M. Y. Hussaini, A. Quarteroni and T. A.
% Zang, _Spectral Methods_, 2 vols., Springer, 2006 and 2007.
%
% [Henrici 1986] P. Henrici, _Applied and Computational Complex Analysis_,
% vol. 3, Wiley, 1986.
%
% [Hesthaven et. al. 2007] J. S. Hesthaven, S. Gottlieb, and D. Gottlieb, 
% _Spectral Methods for Time-Dependent Problems_, Cambridge U. Press, 2007.
%
% [Mcleod 2014] K. Mcleod, FOURFUN: A new system for automatic computations using
% Fourier expansions. Submitted, (2014).
%
% [Trefethen 2000] L. N. Trefethen, _Spectral Methods in MATLAB_, SIAM, 2000.
%
% [Van Loan 1992] C. Van Loan, _Computational Frameworks for the Fast
% Fourier Transform_, SIAM, 1992.
%
% [Zygmund 1959] A. Zygmund, _Trignometric Series_, Cambridge U. Press, 1959.
