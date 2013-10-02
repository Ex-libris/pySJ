__all__ = ["DFData", "SmallMolecule", "CornerHole",\
        "Adatom", "AverageCurve"]

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
        self.df = filtfilt(b, a, self.df)
        return self.df
    
    def make_plot(self):
        plt.plot(self.z,self.df)
         
        
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
        self.force = np.zeros(np.size(self.df)-2)
        #calculate the simple derivative as (df[i+1]-df[i])/(z[i+1]-z[i])
        dz = np.diff(self.df)/np.diff(self.z)
    
        for i in range(np.size(self.z)-2):
            #integration range z_prime offset by one
            z_diff = z[i+1:-1] - z[i]
            #set up the integrand and then integrate with trapezoid
            integrand =\
            (1+((np.sqrt(self.A))/(8*np.sqrt(np.pi*(z_diff)))))*self.df[i+1:-1]-\
        (self.A**(3.0/2)/(np.sqrt(2*(z_diff))))*dz[i+1:]
        
            integral = np.trapz(integrand, z_diff+self.z[i])
        
            #add some correction factors after [JS]
            corr_1 = self.df[i]*np.diff(self.z)[i]
            corr_2 = 2*(np.sqrt(self.A)/(8*np.sqrt(np.pi)))*\
                    self.df[i]*np.sqrt(np.diff(self.z)[i])
            corr_3 = (-2)*((self.A**(3.0/2))/np.sqrt(2))*\
            dz[i]*np.sqrt(np.diff(self.z)[i])
        
            #make the total force and add it to the force array
            self.force[i] = ((2*self.k)/self.f0) *(corr_1+corr_2+corr_3+integral)

        self.distance = z[:np.size(self.force)]
        return self.distance, self.force
        
    

        
class SmallMolecule(DFData):
    """
    inherits from DFData
    """        
    def remove_long_range(self, corner_hole):
        self.df = self.df-corner_hole.df
        return self.df
        
class Adatom(SmallMolecule):
    """
    inherits from SmallMolecule
    """
        
        
class CornerHole(DFData):
    """
    inherits from DFData
    """        
    def shift(self, shift=0.0):
        self.z = self.z - shift
        return self.z

class AverageCurve(object):
    """
    A class to hold averaged force curves.
    :param others:
        an array of DFData-type objects 
    """
    def __init__(self, others):
        self.force_curves = others
        print self.force_curves[0]
        
    
    def average(self):
        self.average = []
        #print "test"
        #print np.size(self.force_curves)
        
        for i in range(np.size(self.force_curves[0].df)):
            sum = 0
            for curve in self.force_curves:
                sum += curve.df[i]   
            self.average.append(sum/np.size(self.force_curves))
        return self.average  
          
    def make_plot(self):
        plt.plot(self.force_curves[0].z,self.average)

                    
        
            
if __name__ == '__main__':
    
    folder_root = "C:\\Users\\bendrevniok\\Dropbox\\"
    folder_root += "DrevniokShareWithDrevniok\\benzene\\df_data\\"
    folder_root += "30June2011\\" 
    #folder_root = ""
    A = 2.5e-10
    k = 1500
    f0 = 56400.0
    fs = 18.3
    order = 1
    freq = 0.1
    
    z,df = datastrip.read_omicron_data(folder_root+"sm\\"+sys.argv[1])
    sm = SmallMolecule(z,df, A, k, f0, fs)
    sm.filtfilt(order,freq)
    #plt.plot(sm.z, sm.ff,sm.z,sm.df)
    
    z,df = datastrip.read_omicron_data(folder_root+"ch\\"+sys.argv[2])
    ch = CornerHole(z,df, A, k, f0, fs)
    ch.filtfilt(order,freq)
      
    
    z,df = datastrip.read_omicron_data(folder_root+"ad\\"+sys.argv[3])
    ad =Adatom(z,df, A, k, f0, fs)
    ad.filtfilt(order,freq)
    
    
    sm.remove_long_range(ch)
    ad.remove_long_range(ch)
    #plt.plot(sm.z, sm.df)
    
    sm.sader_jarvis()
    ad.sader_jarvis()
    plt.plot(sm.distance, sm.force,ad.distance, ad.force)
    
    
    
    plt.show()
    