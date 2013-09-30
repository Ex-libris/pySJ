"""
frequency to force conversion

"""

__all__ = ['ftf']


import numpy as np


def ftf(df, z, A, k, f0, f_scale):
    """
    Use the Sader-Jarvis inversion [JS] to obtain Force(z) from frequency(z)
    Inspired by MATLAB verion by [JW].
    
    Parameters
    ----------
    df : array_like 
        frequency in volts
    z : array_like
        input array length n
        tip position in meters in order from closest approach 
        to furthest from sample. 
        must have z[i] < z[i+1]
    A : float
        amplitude of oscillation in m
    k : float
        spring constant of cantilever
    f0 : float
        resonant frequency of cantilever used for experiment
    f_scale : float
        scaling factor in hz/volt
    
    Returns
    -------
    force : ndarray
            length m array of forces from the frequencies
    frequency : ndarray
                length m array of frequencies lined up with force
    zout : ndarray
        length m array of cantilever positions lined up with force
        
    
    Raises
    ------
    IndexError
    if length of frequency != z
    
    References
    ----------
    .. [JS] John E. Sader, Takayuki Uchihashi, Michael J. Higgins, Alan Farrell, 
    Yoshikazu Nakayama and Suzanne P. Jarvis 
    "Quantitative force measurements using frequency modulation atomic force
    microscopy - theoretical foundations" 
    Nanotechnology, 16 S94-S101 (2005)
    http://www.ampc.ms.unimelb.edu.au/afm/ref_files/Sader_NANO_2005.pdf
    
    .. [JW] Joachim Welker, Esther Illek, Franz J. Giessibl
    "Analysis of force-deconvolution methods in frequency-modulation 
    atomic force microscopy"
    Beilstein J. Nanotechnol. 3, 238 (2012)
    http://www.beilstein-journals.org/bjnano/content/pdf/2190-4286-3-27.pdf
    
    Example
    -------
    
    """
    
    #make input array in numpy array
    df = np.asarray(df)*f_scale
    z = np.asarray(z)
    force = np.zeros(np.size(df)-2)
    #calculate the simple derivative as (df[i+1]-df[i])/(z[i+1]-z[i])
    ddf_dz = np.diff(df)/np.diff(z)
    
    for i in range(np.size(z)-2):
        #integration range z_prime offset by one
        z_diff = z[i+1:-1] - z[i]
        #set up the integrand and then integrate with trapezoid
        integrand =\
        (1+((np.sqrt(A))/(8*np.sqrt(np.pi*(z_diff)))))*df[i+1:-1]-\
        (A**(3.0/2)/(np.sqrt(2*(z_diff))))*ddf_dz[i+1::]
        
        integral = np.trapz(integrand, z_diff+z[i])
        
        #add some correction factors after [JS]
        corr_1 = df[i]*np.diff(z)[i]
        corr_2 = 2*(np.sqrt(A)/(8*np.sqrt(np.pi)))*df[i]*np.sqrt(np.diff(z)[i])
        corr_3 = (-2)*((A**(3.0/2))/np.sqrt(2))\
        *ddf_dz[i]*np.sqrt(np.diff(z)[i])
        
        #make the total force and add it to the force array
        force[i] = ((2*k)/f0) *(corr_1+corr_2+corr_3+integral)

    zout = z[:np.size(force)]
    return zout, force
    
    
if __name__ == '__main__':
    print "this runs"