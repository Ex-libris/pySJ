"""
read in frequency v distance files
"""

__all__ = ['read_omicron_data', 'get_meta_data']

import csv
import re
import numpy as np

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

def get_meta_data(file_name):
    meta_data = []
    for line in file_name:
        if line.startswith('#'):
            yield line
def get_meta_data2(lines):
    meta_data = []
    for line in lines:
        if line.startswith('#'):
            line = re.sub('#', '', line)
            
            meta_data.append(line.strip("\r\n"))
    return meta_data
            
     
if __name__ == '__main__':
    folder_name = "C:\\Users\\bendrevniok\\Dropbox\\\
DrevniokShareWithDrevniok\\benzene\\df_data\\30June2011\\"
    file_name = "sm\\df_90.txt"
    file_name  = folder_name+file_name
    

    z,df = read_omicron_data(file_name)
    #print z  
    with open(file_name,'rb') as f:
        lines = f.readlines()
        meta_data = get_meta_data2(lines)

    print meta_data[11]        
    
        
    