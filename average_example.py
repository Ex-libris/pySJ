import datastrip
import df_to_force
import sys
import numpy as np
from scipy.signal import filtfilt, butter
import matplotlib.pyplot as plt




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
    freq = 0.8
    
    
    sm = []
    files = sys.argv[1:]
    
    for i in range(np.size(files)):
        z,df = datastrip.read_omicron_data(folder_root+"sm\\"+files[i])
        sm.append(df_to_force.SmallMolecule(z,df, A, k, f0, fs))
        sm[i].filtfilt(order,freq)
        sm[i].make_plot()
    
    sma = df_to_force.AverageCurve(sm)
    sma.average()
    sma.make_plot()
    plt.show()
    
    
    