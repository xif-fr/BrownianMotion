import numpy as np
from math import *
π = np.pi
import scipy.special as ss
import scipy.integrate as sint
import mpmath


def convolve_around_center (func1, func2, N1, Nout, Δx, x_center=0):
	u"""
	Convolve two functions func1 and func2, with func1 decreasing away from 0 (convolution kernel) :
	(func1*func2)(x) = ∫ dx1 func2(x-x1) func1(x1)
					 ≃ Δx ∑ func2(x-k⋅Δx) func1(k⋅Δx) from k=-N1 to +N1
	-> Only 2⋅N1+1 points of func1 are sampled around 0, while func2 is evaluated as needed (around x=x_center).
	The result is ( x, (func1*func2)(x) ) with X = [x_center-Nout⋅Δx, x_center+Nout⋅Δx].
	Typically, Nout << N1. Nout can be 0, and in that case, the result is simply (func1*func2)(x_center).
	"""
	# samples of func1
	X1 = np.linspace(-N1*Δx, +N1*Δx, 2*N1+1)
	Y1 = func1( X1 )
	# samples of func2
	X2 = x_center + np.linspace((-N1-Nout)*Δx, (+N1+Nout)*Δx, 2*(N1+Nout)+1)
	Y2 = func2( X2 )
	# output
	Conv_x = x_center + np.linspace(-Nout*Δx, +Nout*Δx, 2*Nout+1)
	Conv = np.zeros(2*Nout+1)
	for i in range(2*Nout+1):
		# pas optimal car ré-évaluation inutile de func2 :
		#  Y2 = func2( Conv_x[i] - X1 )
		#  Conv[i] = np.sum( Y1 * Y2 ) * Δx
		# mieux :
		Y2loc = Y2[i:i+2*N1+1]
		Conv[i] = np.sum( Y1 * Y2loc ) * Δx
	return Conv_x, Conv

def distr_x0_harmonic (x, σ):
	return np.exp(-(x/σ)**2/2)/sqrt(2*π)/σ

#----------------------------------------------------------------
# No resetting

def fpt_free_survival (L, t, D, σ):
	survdist_σ0 = lambda x, t: ss.erf( x/np.sqrt(4*D*t) ) * np.heaviside(x,0.5)
	if σ == 0:
		return survdist_σ0(L,t)
	else:
		assert np.isscalar(L)
		distr_x0 = lambda x: distr_x0_harmonic(x, σ)
		def ps (t):
			surv_f = lambda x: survdist_σ0(x, t)
			return convolve_around_center(distr_x0, surv_f, x_center=L, N1=1000+int(500*sqrt(4*D*np.max(t))/σ), Nout=0, Δx=0.01*σ)[1]
		return np.vectorize(ps)(t)

def fpt_free_distrib (t, x_targ):
	if σ == 0:
		return x_targ/(2*np.sqrt(π*D*t**3)) * np.exp(-x_targ**2/(4*D*t))
	else:
		pass

def fpt_2d_free_survival (R, t, D, Rtol, σ, regularize=True, split_domain=True):
	if σ == 0:
		a = Rtol/R
		c = R/np.sqrt(4*D*t)
		f = lambda x, a,c: np.exp(-x**2/(4*a**2*c**2)) / x * (ss.y0(x/a)*ss.j0(x)-ss.j0(x/a)*ss.y0(x)) / (ss.y0(x)**2+ss.j0(x)**2)
		if regularize:
			# regularization of the divergence of f at x=0 by substracting the leading-order term,
			# which is, amazingly, integrable analytically; this allows the integrator to better behave;
			# splitting the domain in two does improve the result a tiny bit;
			# (but this method seems to lead to a slight overestimation of the survival proba, if the langevin simulations are accurate)
			f_reg = lambda x, a,c: f(x,a,c) - 1/x * 2/π * log(1/a) / (1 + 4/π**2 * (np.euler_gamma+np.log(x/2))**2)
			if split_domain: ps0 = lambda a,c: 2*log(1/a) + 2/π * ( sint.quad(f_reg, 0, 1, args=(a,c), epsabs=1e-6, limit=1000)[0] + sint.quad(f_reg, 1, +np.inf, args=(a,c), epsabs=1e-5, limit=1000)[0] )
			else:            ps0 = lambda a,c: 2*log(1/a) + 2/π * ( sint.quad(f_reg, 0, +np.inf, args=(a,c), epsabs=1e-5, limit=1000)[0] )
		else:
			# splitting the domain in two (one near zero where there is a singularity, the other to infinity)
			# allows to use to integration methods, one on the finite domain which treats the singularity well
			# and the other which treats the rest of the infinite domain without singularity
			if split_domain: ps0 = lambda a,c: 2/π * ( sint.quad(f, 0, 0.1, args=(a,c), epsabs=1e-4, limit=1000)[0] + sint.quad(f, 0.1, +np.inf, args=(a,c), epsabs=1e-6, limit=1000)[0] )
			else:            ps0 = lambda a,c: 2/π * sint.quad(f, 0, +np.inf, args=(a,c), epsabs=1e-5, limit=1000)[0]
		return np.vectorize( lambda a,c: (ps0(a,c) if a < 0.999 else 0.) )(a,c)
	else:
		# just a convolution of a guassian with the σ=0 curve
		pass

#----------------------------------------------------------------
# Poissonian reset

def fpt_poisson_c (α, D, L):
	return sqrt(α/D)*L

def fpt_poisson_inverselapl (x, t, α, D, σ, fpt):
	mpmath.mp.dps = 30
	x = np.atleast_1d(x)
	t = np.atleast_1d(t)
	P = np.zeros((len(x),len(t)))

	sqrt2 = mpmath.sqrt(2)
	if fpt:
		ret_psr_lp = lambda psr,s: 1 - s*psr  # p(tf) = - d/dt psr
	else:
		ret_psr_lp = lambda psr,s: psr

	for i in range(len(x)):
		if σ == 0:
			def ps0_lp (κ, s):
				return (1 - mpmath.exp(-κ * x[i])) / s
		else:
			b = x[i] / σ
			def ps0_lp (κ, s):
				k = σ * κ
				return (1 - mpmath.exp(k**2/2)/2 * ( mpmath.exp(+κ*x[i]) * mpmath.erfc((b+k)/sqrt2)
				                                   + mpmath.exp(-κ*x[i]) * (1+mpmath.erf((b-k)/sqrt2)) ) ) / s
		def psr_lp (s):
			κ = mpmath.sqrt( (α+s) / D )
			ps0 = ps0_lp(κ, s=α+s)
			psr = ps0 / (1 - α*ps0)
			return ret_psr_lp(psr, s)

		for j in range(len(t)):
			if x[i] < 0:
				P[i,j] = 0
			else:
				P[i,j] = mpmath.invertlaplace(psr_lp, t[j], method='talbot', degree=20)
	return np.squeeze(P)

def fpt_poisson_survival (x, t, α, D, σ):
	return fpt_poisson_inverselapl(x, t, α, D, σ, False)

def fpt_poisson_distrib (x, t, α, D, σ):
	return fpt_poisson_inverselapl(x, t, α, D, σ, True)

def fpt_poisson_tau (b, c):
	if np.all(np.isinf(b)):
		return 4/c**2 * ( np.exp(c) - 1 )
	else:
		return 4/c**2 * ( (2*np.exp(-c**2/2/b**2)) / ( np.exp(c)*ss.erfc((c/b+b)/sqrt(2)) + np.exp(-c)*ss.erfc((c/b-b)/sqrt(2)) ) - 1 )

def fpt_2d_poisson_tau (b, c, a, do_warn_err=False):
	a = np.fmin(a, 1-1e-10)
	def func (a,b,c):
		if b > 18:
			if not np.isinf(b):
				print("warning : approximating b={:.3f} by b=inf".format(b)) 
			return ss.k0(a*c) / ss.k0(c) - 1
		else:
			# regularization, not needed :
			# f = lambda z, b,c: z * np.exp(-z**2/2) * ( ss.k0(c/b*z) * ss.i0(b*z) + np.log(z) )
			# d = -(a*b)**2/2
			# np.exp(-b**2/2) * sint.quad(f, a*b, max(10,2*b), args=(b,c), epsrel=1e-8)[0] - np.exp(d)*np.log(a*b) + ss.expi(d)/2
			f = lambda z, b,c: z * np.exp(-b**2/2-z**2/2) * ( ss.k0(c/b*z) * ss.i0(b*z) )
			I, Ierr = sint.quad(f, a*b, max(10,2*b), args=(b,c), epsrel=1e-8)
			x = ss.k0(a*c) / I - 1
			if do_warn_err:
				xp = ss.k0(a*c) / (I+Ierr) - 1
				xm = ss.k0(a*c) / (I-Ierr) - 1
				if abs((xp-xm)/x) > 1e-2:
					print("warning : rel. error can be >1% for b={:.3f}, c={:.3f}, a={:.3f}".format(b,c,a))          
			return x
	return 4/c**2 * np.vectorize(func)(a,b,c)

#----------------------------------------------------------------
# Periodical reset

def fpt_periodic_c (rT, D, L):
	return L/sqrt(4*D*rT)

def fpt_periodic_tau (b, c):
	if np.all(np.isinf(b)):
		return ( ss.erf(c) + 2*c*(np.exp(-c**2)/sqrt(π)-c*ss.erfc(c)) ) / ss.erfc(c) / c**2
	else:
		int_exp_erf = lambda v,b,c: sint.quad( lambda u, v,b,c: np.exp(-u**2/2) * ss.erf(c/np.sqrt(v)*np.abs(1-u/b)), -np.inf, +np.inf, args=(v,b,c), epsrel=1e-1 )[0]
		int_exp_erf = np.vectorize( int_exp_erf, excluded=(1,2) )
		int_v = np.vectorize( lambda b,c: sint.quad( int_exp_erf, 0, 1, args=(b,c), epsrel=1e-1 )[0] )
		int_exp_erfc = lambda b,c: sint.quad( lambda u, b,c: np.exp(-u**2/2) * ss.erfc(c*np.abs(1-u/b)), -np.inf, +np.inf, args=(b,c), epsrel=1e-3 )[0]
		int_exp_erfc = np.vectorize( int_exp_erfc )
		return int_v(b,c) / int_exp_erfc(b,c) / c**2

int_exp_erf = lambda b,c: sint.quad( lambda u, b,c: np.exp(-u**2/2) * ss.erf(c*np.abs(1-u/b)), -np.inf, +np.inf, args=(b,c), epsrel=1e-3 )[0]
int_exp_erf = np.vectorize( int_exp_erf )

def fpt_periodic_survival (t, rT, b, c):
	global int_exp_erf
	k = np.floor( t / rT )
	if np.all(np.isinf(b)):
		return ss.erf(c)**k * ss.erf(c*np.sqrt(rT/(t-k*rT)))
	else:
		int_exp_erf_kt = lambda b,c,k,t: sint.quad( lambda u, b,c,k,t: np.exp(-u**2/2) * ss.erf(c*np.abs(1-u/b)*np.sqrt(rT/(t-k*rT))), -np.inf, +np.inf, args=(b,c,k,t), epsrel=1e-3 )[0]
		int_exp_erf_kt = np.vectorize( int_exp_erf_kt )
		return (1/sqrt(2*π))**(k+1) * int_exp_erf(b,c)**k * int_exp_erf_kt(b,c,k,t)

def fpt_periodic_disrib (t, rT, b, c):
	global int_exp_erf
	k = np.floor( t / rT )
	prefact_dt_erf = lambda k,t: ((k+1)*rT-t)/sqrt(π*rT)*(t-k*rT)**(-3/2)
	if np.all(np.isinf(b)):
		return c * prefact_dt_erf(k,t) * ss.erf(c)**k * np.exp(-c**2*rT/(t-k*rT))
	else:
		int_exp_dterf = lambda b,c,k,t: sint.quad( lambda u, b,c,k,t: np.abs(1-u/b) * np.exp(-u**2/2 -c**2*(1-u/b)**2*rT/(t-k*rT)), -np.inf, +np.inf, args=(b,c,k,t), epsrel=1e-3 )[0]
		int_exp_dterf = np.vectorize( int_exp_dterf )
		return c * (1/sqrt(2*π))**(k+1) * int_exp_erf(b,c)**k * prefact_dt_erf(k,t) * int_exp_dterf(b,c,k,t)
