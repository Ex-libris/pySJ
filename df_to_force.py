__all__ = ["DFData"]

__version__ = "0.1"

import sys
import numpy as np
import datastrip
from scipy.signal import filtfilt, butter
import matplotlib.pyplot as plt



class DFData(object):
    """
    The basic class for a frequency shift curve and height data object
    :param z:
        numpy array of length n that contains tip-sample distance 
    :param df: 
        numpy array of length n that contains frequency shift
    :param A (optional):
        cantilever oscillation amplitude in m
    :param k (optional):
        cantilever stiffness in N/m
    :param f0 (optional):
        fundamental frequency of cantilever in Hz
    :param fs (optional):
        conversion from voltage to Hz in Hz/V      
              
                        
    """
    def __init__(self, z, df, 
                A = 2.5e-10, k = 1800, 
                f0 = 50000, fs = 15):
            self.z = z
            self.df = df*fs
            self.A = A
            self.k = k
            self.f0 = f0
            self.fs = fs  
            
    def filtfilt(self,order = 3, freq = 0.05):
        """
        a backwards and forwards filter from scipy added to this class
        uses a nth-order butterworth filter to provide a non-phase delayed 
        signal smoothing
        """
        b, a = butter(order, freq)
        # Use filtfilt to apply the filter.
        self.ff = filtfilt(b, a, self.df)
        return self.ff
        
    def sader_jarvis(self):
        """
        implement the Sader-Jarvis inversion
        
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
        """
            #make input array in numpy array
        self.force = np.zeros(np.size(self.sr)-2)
        #calculate the simple derivative as (df[i+1]-df[i])/(z[i+1]-z[i])
        dz = np.diff(self.sr)/np.diff(self.z)
    
        for i in range(np.size(self.z)-2):
            #integration range z_prime offset by one
            z_diff = z[i+1:-1] - z[i]
            #set up the integrand and then integrate with trapezoid
            integrand =\
            (1+((np.sqrt(self.A))/(8*np.sqrt(np.pi*(z_diff)))))*self.ff[i+1:-1]-\
        (self.A**(3.0/2)/(np.sqrt(2*(z_diff))))*dz[i+1::]
        
            integral = np.trapz(integrand, z_diff+self.z[i])
        
            #add some correction factors after [JS]
            corr_1 = self.sr[i]*np.diff(self.z)[i]
            corr_2 = 2*(np.sqrt(self.A)/(8*np.sqrt(np.pi)))*\
                    self.sr[i]*np.sqrt(np.diff(self.z)[i])
            corr_3 = (-2)*((self.A**(3.0/2))/np.sqrt(2))*\
            dz[i]*np.sqrt(np.diff(z)[i])
        
            #make the total force and add it to the force array
            self.force[i] = ((2*self.k)/self.f0) *(corr_1+corr_2+corr_3+integral)

        self.forcez = z[:np.size(self.force)]
        return self.forcez, self.force
        
class SmallMolecule(DFData):
    """
    inherits from DFData
    """        
    def remove_long_range(self, corner_hole):
        self.sr = self.ff-corner_hole.ff
        return self.sr
        
        
class CornerHole(DFData):
    """
    inherits from DFData
    """        
    def shift(self, shift=0.0):
        self.zs = self.z - shift
        return self.zs
        
#class DFDifference(object):
#    """
#    extension of DFDdata class that will allow one to subtract 
#    one :DFData: object from another
#    This is needed for removing the long-range forces in our experiments. 
#    """
            
if __name__ == '__main__':
    
    folder_root = "C:\\Users\\bendrevniok\\Dropbox\\DrevniokShareWithDrevniok\\\
benzene\\df_data\\30June2011\\"
    A = 2.7e-10
    k = 1.8e3
    f0 = 5.64e4
    fs = 15
    
    order = 3
    freq = 0.1
    
    z,df = datastrip.read_omicron_data(folder_root+"sm\\"+sys.argv[1])
    sm = SmallMolecule(z,df, A, k, f0, fs)
    sm.filtfilt(order,freq)
    #plt.plot(sm.z, sm.ff,sm.z,sm.df)
    
    z,df = datastrip.read_omicron_data(folder_root+"ch\\"+sys.argv[2])
    ch = CornerHole(z,df, A, k, f0, fs)
    ch.filtfilt(order,freq)
    
    sm.remove_long_range(ch)
    #plt.plot(sm.z, sm.sr)
    
    sm.sader_jarvis()
    plt.plot(sm.zforce, sm.force)
    
    
    plt.show()
    