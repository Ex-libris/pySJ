# -*- coding: utf-8 -*-
__all__ = ["DFData", "SmallMolecule", "CornerHole", "Adatom", "AverageCurve", "ReadOmicronXY"]

__version__ = "0.1"
import sys
import numpy as np
from scipy.signal import filtfilt, butter
import matplotlib.pyplot as plt
import matplotlib as mpl



class DFData(object):
    """
    The basic class for a frequency shift curve and height data object
    :param z:
        numpy array of length n that contains tip-sample distance 
    :param df: 
        numpy array of length n that contains frequency shift
    :param a (optional):
        cantilever oscillation amplitude in m
    :param k (optional):
        cantilever stiffness in N/m
    :param f0 (optional):
        fundamental frequency of cantilever in Hz
    :param fs (optional):
        conversion from voltage to Hz in Hz/V
    """

    def __init__(self, distance, frequency_shift,
                 amplitude, spring_constant,
                 resonant_frequency, voltage_to_frequency,data_name, 
                 force_distance = None, integrated_force = None):
        self.z = distance
        self.df = frequency_shift * float(voltage_to_frequency)
        self.a = float(amplitude)
        self.k = float(spring_constant)
        self.f0 = float(resonant_frequency)
        self.distance = force_distance
        self.force = integrated_force
        self.name = data_name

    def filtfilt(self, filter_order=3, cutoff_freq=0.05):
        """
        a backwards and forwards filter from scipy added to this class
        uses a nth-order butterworth filter to provide a non-phase delayed
        signal smoothing
        """
        b, a = butter(filter_order, cutoff_freq)
        # Use filtfilt to apply the filter.
        self.df = filtfilt(b, a, self.df)
        return self.df

    def make_plot(self, xlim, ylim):
        plt.xlim((xlim[0],xlim[1]))
        plt.ylim((ylim[0],ylim[1]))
        plt.rc("axes", linewidth=3.0)
        plt.rc(("xtick","ytick","axes"), labelsize=20.0)

        if self.distance is not None:
            plt.ylabel('short-range force [nN]')
            plt.xlabel('distance [nm]')
            plt.plot(self.distance*1e9, self.force*1e9,
                    linewidth = 3.0, 
                    label = self.legend+" "+self.name+"f")
        else: 
            plt.ylabel('frequency shift [Hz]')
            plt.xlabel('distance [nm]')
            plt.plot(self.z*1e9, self.df-np.average(self.df[-50:]),
                    linewidth = 3.0,
                    label = self.legend+" "+self.name+"df")
            
    def sader_jarvis(self, difference = 5):
        """
        implement the Sader-Jarvis inversion
        
            References
        ----------
        .. [JS] John E. Sader, Takayuki Uchihashi, Michael J. Higgins,
        Alan Farrell, Yoshikazu Nakayama and Suzanne P. Jarvis
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
        #get the spacing in z
        h = self.z[1]-self.z[0]
        
        if difference == 5:
            #preallocate derivative array. It will be n-4 in length 
            #preallocate the force numpy array same length
            self.force = np.zeros(np.size(self.z) - 4)
            self.distance = self.z[2:-2]  
            derivative = np.zeros(np.size(self.z))
            weights = np.array([1.0/12, -2.0/3, 0, 2.0/3, -1.0/12])
                
            for i in range(2,np.size(self.df)-2):
                derivative[i] = (np.dot(weights, self.df[i-2:i+3])/h)  
       
        for i in range(2,np.size(self.z)-2):
            t = self.z[i+1:-2]
            df = self.df[i+1:-2]
            deriv = derivative[i+1:-2]
            integrand = ((1 + (self.a**(1.0/2))/
                    (8*((np.pi*(t-self.z[i]))**(1/2.0)))) * df
                    -(((self.a**(3.0/2))/(np.sqrt(2*(t-self.z[i]))))*deriv))

            
            integral = np.trapz(integrand,t)  
            #add some correction factors after [JS]
            corr_1 = self.df[i+2] * h
            corr_2 = (2 * (np.sqrt(self.a) / (8 * np.sqrt(np.pi))) 
                    * self.df[i+2] * np.sqrt(h))
            corr_3 = ((-2) * ((self.a ** (3.0 / 2)) / np.sqrt(2)) 
                    * derivative[i+2] * np.sqrt(h))
            self.force[i-2] = ((2 * self.k)/self.f0)*(
                    integral+corr_1+corr_2+corr_3)  
        return self.distance, self.force
        
    def normalize(self):
        """
        Normalize a force curve to the absolute value
        of the minimum value of a set of force curves    
        """
        self.force = self.force/np.abs(np.min(self.force))
    def shift_to_zero(self):
        """
        shift a curve to have 0Hz frequency shift 
        when it is far from the surface
        """
        self.df -= np.average(self.df[-50:])
        return self.df
        
class SmallMolecule(DFData):
    """
    inherits from DFData
    """
    def __init__(self, distance, frequency_shift,
                 amplitude, spring_constant,
                 resonant_frequency, voltage_to_frequency,data_name, 
                 force_distance = None, integrated_force = None,
                 legend_name = "sm"):
        DFData.__init__(self,distance, frequency_shift,
                        amplitude, spring_constant,
                        resonant_frequency, voltage_to_frequency,data_name, 
                        force_distance = None, integrated_force = None)
        self.legend = legend_name

    def remove_long_range(self, corner_hole):
        self.df -= corner_hole.df
        return self.df


class Adatom(DFData):
    """
    inherits from SmallMolecule
    """
    def __init__(self, distance, frequency_shift,
                 amplitude, spring_constant,
                 resonant_frequency, voltage_to_frequency,data_name, 
                 force_distance = None, integrated_force = None,
                 legend_name = "ad"):
                 
        DFData.__init__(self,distance, frequency_shift,
                        amplitude, spring_constant,
                        resonant_frequency, voltage_to_frequency,data_name, 
                        force_distance = None, integrated_force = None)
        self.legend = legend_name

    def remove_long_range(self, corner_hole):
        self.df -= corner_hole.df
        return self.df

class CornerHole(DFData):
    """
    inherits from DFData
    """
    def __init__(self, distance, frequency_shift,
                 amplitude, spring_constant,
                 resonant_frequency, voltage_to_frequency,data_name, 
                 force_distance = None, integrated_force = None,
                 legend_name = "ch"):
                 
        DFData.__init__(self,distance, frequency_shift,
                        amplitude, spring_constant,
                        resonant_frequency, voltage_to_frequency,data_name, 
                        force_distance = None, integrated_force = None)
                        
        self.legend = legend_name
        
    def shift(self, shift=0.0):
        self.z -= shift
        return self.z


class AverageCurve(object):
    """
    A class to hold averaged force curves.
    :param others:
        an array of DFData-type objects 
    """

    def __init__(self, others, average_curve = None):
        self.force_curves = others
        self.average_df = None
        self.average_f = None
        self.normalizer = 1.0
        self.ax = None
        self.fig = None

    def make_average_df(self):
        """
        Find the average force curve of multiple force curve objects
        """
        self.average_df = np.zeros(np.size(self.force_curves[0].z))
                
        for i in range(np.size(self.force_curves[0].df)):
            total = 0
            for curve in self.force_curves:
                total += curve.df[i]
            self.average_df[i] = total
        self.average_df = self.average_df/len(self.force_curves)            
        return self.average_df
        
    def make_average_f(self):
        """
        Find the average force curve of multiple force curve objects
        """
        self.average_f = np.zeros(np.size(self.force_curves[0].distance))
                
        for i in range(np.size(self.force_curves[0].distance)):
            total = 0
            for curve in self.force_curves:
                total += curve.force[i]
            self.average_f[i] = total
        self.average_f = self.average_f/len(self.force_curves)            
        return self.average_f
    
    
    def find_minimum(self):
        """
        find the minimum value of the set of force curves 
        """
        minima = np.zeros(len(self.force_curves))
        for i in range(len(self.force_curves)):
            minima[i] = np.min(self.force_curves[i].force)
        self.normalizer = np.abs(np.min(minima))
            
    def make_plot(self, key = 'f'):
        if key == 'f':
            plt.plot(self.force_curves[0].distance, 
                    self.average_f, label = "average force")
        else:
            plt.plot(self.force_curves[0].z, self.average_df,
                    label = "average df")
                    
    def plot_all_df(self,filename):
        """
        plot all plots contained in the average class
        That will make all the z,df or distance, f plots
        be plotted.
        
        """ 
        self.make_average_df()
        self.fig,self.ax = plt.subplots() 
        plt.xlim((np.min(self.force_curves[0].z)*1e9 - 0.06,
                np.max(self.force_curves[0].z)*1e9 + 0.05))
        plt.ylim((-60,2))
        plt.ylabel('frequency shift [Hz]')
        plt.xlabel('distance [nm]')
        size = len(self.force_curves)
        for i in range(len(self.force_curves)):
            if self.force_curves[i].legend == "sm":
                col = [1.0,0,0,0.9]
            elif self.force_curves[i].legend == "ad":
                col = [0.0,1.0,0.0,0.9]
            elif self.force_curves[i].legend == "ch":
                col = [0,0,1.0,0.9]
            rgb = [x * float(i)/(size+2) for x in col] 
            self.ax.plot(self.force_curves[i].z*1e9,self.force_curves[i].df,
                    'd',
                    markersize=6,
                    color = rgb,
                    label = self.force_curves[i].legend+" "
                    +self.force_curves[i].name+" df") 
        col[3] = 1.0
        self.ax.plot(self.force_curves[0].z*1e9, self.average_df,
                color = col,
                linewidth = 4.0,
                label = self.force_curves[i].legend +" average df",)
        self.ax.legend(loc = 'lower right')
       
        self.fig.savefig(filename+
                self.force_curves[i].legend+"\\"+self.force_curves[i].name
                +".pdf",
                transparent=True)
        self.fig.savefig(filename+
                self.force_curves[i].legend+"\\"+self.force_curves[i].name
                +".png",
                transparent=True)

class ReadOmicronXY(object):
    """
    Class to  hold omicron XY files from vernissage
    :param folder_root (optional):
        string for root folder
    :param filename (optional):
        string for the specific frequency-shift data that will be 
        examined
    """

    def __init__(self, folder_root="test\\", filename="test",
            xydata = None,z_in = None, z_out = None, 
            df_in = None, df_out = None):
        self.folder_root = folder_root
        self.filename = filename
        self.data = xydata
        self.z_approach = z_in
        self.z_retract = z_out
        self.df_approach = df_in
        self.df_retract = df_out
        self.name = filename[6:-4]

    def change_root(self, root):
        """
        change the object root folder name
    
        Parameters
        ----------
        root : string\of\path\to\data
        
        Returns
        -------
        self.folder_root : updated folder root
        """
        self.folder_root = root
        self.load_and_get()
        return self.folder_root

    def change_file(self, file):
        """
        change the object filename
    
        Parameters
        ----------
        root : string\of\path\to\data
        
        Returns
        -------
        self.folder_root : updated folder root
        """
        self.filename = file
        self.load_and_get()
        return self.filename

    def load_file(self):
        """
        load the file and create the numpy arrays for z and df
    
        Parameters
        ----------
        
        
        Returns
        -------
        self.z : numpy array
                tip-sample distance
        self.df : numpy array
                frequency shifts
        self.data : numpy array
                2n array of tip-sample distance and frequency shifts
        """
        self.data = np.loadtxt(self.folder_root + self.filename)
        self.get_column()

    def get_column(self):
        """
        get the distance and frequency shift from omicron data
        and put them into individual variables 
        Addtionally make all the data run from closest approach to furthest.
        That is, from interacting to not interacting.
        
        Parameters
        ----------
        
        
        Returns
        -------
        self.z_approach : numpy array
                       approach array of tip-sample  distance for force curve
                       organized from closest to furthest from sample
        self.z_retract : numpy array
                        retract array of tip-sample  distance for force curve
                        organized from closest to furthest from sample
        self.df_approach : numpy array
                        retract array of frequency shift for force curve
                        organized from closest to furthest from sample
        self.df_retract : numpy array
                        retract array of frequency shift for force curve
                        organized from closest to furthest from sample
            
        """
        self.z_approach = np.zeros(np.size(self.data[:, 0]))
        self.z_retract = np.zeros(np.size(self.data[:, 0]))
        self.df_approach = np.zeros(np.size(self.data[:, 0]))
        self.df_retract = np.zeros(np.size(self.data[:, 0]))

        if np.size(self.data[0:np.argmin(self.data[:, 0]), 0]) < (
            np.size(self.data[:, 0])):
            self.z_approach = self.data[np.size(self.data[:, 0]) / 2 - 1:0:-1, 0]
            self.z_retract = self.data[np.size(self.data[:, 0]) / 2:, 0]
            self.df_approach = self.data[np.size(self.data[:, 0]) / 2 - 1:0:-1, 1]
            self.df_retract = self.data[np.size(self.data[:, 0]) / 2:, 1]
        else:
            self.z_approach = self.data[-1::-1, 0]
            self.df_approach = self.data[-1::-1, 1]

        return self.z_approach, self.z_retract, self.df_approach, self.df_retract

    def load_and_get(self):
        """
        
        """
        self.load_file()
        self.get_column()


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

    data = ReadOmicronXY(folder_root, sys.argv[1])
    data.load_and_get()
    sm1 = SmallMolecule(data.z_approach, data.df_approach, A, k, f0, fs)
    data.change_file(sys.argv[2])
    sm1.filtfilt(order, freq)
    sm2 = SmallMolecule(data.z_approach, data.df_approach, A, k, f0, fs)
    sm2.filtfilt(order, freq)
    plt.plot(sm1.z, sm1.df, sm2.z, sm2.df)
    plt.show()
    