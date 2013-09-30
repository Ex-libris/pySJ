"""
read in frequency v distance files
"""

__all__ = ['read_omicron_data']

import numpy as np
import scipy.signal as signal

def read_omicron_data(file_name):
    """
    use numpy.loadtxt to create data that is acceptable for 
    frequency shift to force conversion
    
    Parameters
    ----------
    file_name : string 
        total path and filename of data 
        
    Returns
    -------
    df : array_like
        frequency shift data array
    
    z : array_like
        cantilever position data
        
    """
    data = np.loadtxt(file_name)
    #make the two arrays and reverse them 
    #if the curve had approach and retract
    #just take the approach for now
    df = data[np.argmin(data[:,0]):0:-1,1]
    z = data[np.argmin(data[:,0]):0:-1,0]
    return z,df
    
    
def simple_average(z,df, stride):
    """
    Calculate a simple average of the data a few blocks at a time
        
    Parameters
    ----------
    z : array_like
        length n
    df : array_like
        length n 
    stride : int
        size of averaging steps        
    
    Returns
    -------
    av_z : array_like
    
    av_df : array_like
    """
    av_z = np.zeros(np.size(z))
    av_df = np.zeros(np.size(z))
    for i in range(0,np.size(z)-stride,stride):
        av_z[i/5] = np.average(z[i:i+stride-1])
        av_df[i/5] = np.average(df[i:i+stride-1])
    av_z = np.trim_zeros(av_z)
    av_df = np.trim_zeros(av_df)
    return av_z, av_df
    


def filtfilt(df, order = 3):
    # Create an order 3 lowpass butterworth filter.
    b, a = signal.butter(order, 0.05)
    # Use filtfilt to apply the filter.
    return signal.filtfilt(b, a, df)
               
     
if __name__ == '__main__':
    folder_name = "C:\\Users\\bendrevniok\\Dropbox\\\
DrevniokShareWithDrevniok\\benzene\\df_data\\mol_small\\"
    file_name = "default_2011Jun29-110255_Custom\
-AFMHybrid_NonContact_QPlus_Twin--28_1-Aux2(Z).txt"
    file_name  = folder_name+file_name
    
    z,df = read_omicron_data(file_name)
    print z
        
        
    