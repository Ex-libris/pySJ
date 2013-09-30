import sj
import datastrip
import matplotlib.pyplot as plt
import os
import sys


folder_root = "C:\\Users\\bendrevniok\\Dropbox\\DrevniokShareWithDrevniok\\\
benzene\\df_data\\30June2011\\"

def molecule(folder_root,argv):
    A = 2.5e-10
    k = 1500.0
    f0 = 56400.0
    fs = 15
    order = 1
    
    ch_z,ch_df = datastrip.read_omicron_data(folder_root+"ch\\"+argv[0])
    ch_ff = datastrip.filtfilt(ch_df,order)
    
    sm_z,sm_df = datastrip.read_omicron_data(folder_root+"sm\\"+argv[1])
    sm_ff = datastrip.filtfilt(sm_df,order)
    
    
    z,force = sj.ftf(sm_ff-ch_ff,ch_z,A,k,f0,fs)
    

    figa = plt.figure()
    ax = figa.add_subplot(111)
    ax.plot(z,force)
    
    figb = plt.figure()
    bx = figb.add_subplot(111)
    bx.plot(ch_z,ch_ff,sm_z, sm_ff) 
    
def adatom(folder_root,argv):
    A = 2.5e-10
    k = 1500.0
    f0 = 56400.0
    fs = 15
    order = 1
    
    ch_z,ch_df = datastrip.read_omicron_data(folder_root+"ch\\"+argv[0])
    ch_ff = datastrip.filtfilt(ch_df,order)
    
    sm_z,sm_df = datastrip.read_omicron_data(folder_root+"ad\\"+argv[1])
    sm_ff = datastrip.filtfilt(sm_df,order)
    
    
    z,force = sj.ftf(sm_ff-ch_ff,ch_z,A,k,f0,fs)
    

    figa = plt.figure()
    ax = figa.add_subplot(111)
    ax.plot(z,force)
    
    figb = plt.figure()
    bx = figb.add_subplot(111)
    bx.plot(ch_z,ch_ff,sm_z, sm_ff) 

    
if __name__ == '__main__':
    molecule(folder_root, sys.argv[1:])
    adatom(folder_root, (sys.argv[1],sys.argv[3]))
    plt.show()
    