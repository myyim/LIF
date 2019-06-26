import scipy
import scipy.integrate
import numpy

def rate(mu=0.06, sig=0.2, V_r=0., V_th=1., t_m=20., t_r=2.):
        """"rate(mu=0.06, sig=0.2, V_r=0., V_th=1., t_m=20., t_r=2.) returns the firing rate of a LIF neuron receiving Gaussian white noise input. The output rate is in Hz. The temporal dimension for other parameters is in ms.
Note: You may need to test the convergence of scipy.integrate.quad(integrand1,0.,50.) by tuning the third (last) argument.
See Yim, Aertsen & Rotter (2013) PRE for details: http://pre.aps.org/abstract/PRE/v87/i3/e032710 """
        if sig == 0:
                if (V_th-mu*t_m)/(V_r-mu*t_m) <= 0:
                        return 0.
                else:
                        return 1000./(t_r-t_m*numpy.log((V_th-mu*t_m)/(V_r-mu*t_m)))
        else:
                y_r = (V_r-mu*t_m)/(sig*numpy.sqrt(t_m))
                y_th = (V_th-mu*t_m)/(sig*numpy.sqrt(t_m))
                def integr(u):
                        if u == 0:
                                return 2*(y_th-y_r)
                        else:
                                return (numpy.exp(-u*u+2*y_th*u)-numpy.exp(-u*u+2*y_r*u))/u
                return 1000./(scipy.integrate.quad(integr,0.,50.)[0]*t_m+t_r)
        # check convergence by tuning the third argument in scipy.integrate.quad
