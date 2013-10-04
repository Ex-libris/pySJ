import df_to_force
import sys
import numpy as np
import matplotlib.pyplot as plt

folder_root = "C:\\Users\\bendrevniok\\Dropbox\\"
folder_root += "DrevniokShareWithDrevniok\\benzene\\df_data\\"
folder_root += "30June2011\\"

A = 2.5e-10
k = 1500
f0 = 56400.0
fs = 18.3
order = 1
freq = 0.8

sm_files = ["df_90.txt", "df_91.txt","df_98.txt", "df_107.txt","df_108.txt","df_113.txt","df_114.txt" ]
ad_files = ["df_86.txt","df_87.txt","df_97.txt","df_109.txt","df_110.txt"]
ch_files = ["df_88.txt","df_89.txt","df_95.txt","df_96.txt","df_111.txt","df_112.txt"]

sm = []
ad = []
ch = []

for file_name in sm_files:
    data = ReadOmicronXY(folder_root, file_name)
    data.load_and_get()
    sm.append(SmallMolecule(data.z_approach, data.df_approach, A, k, f0, fs))
    #sm[i].filtfilt(order, freq)
    sm[i].make_plot()


#sma = df_to_force.AverageCurve(sm)
#sma.average()
#sma.make_plot()
#plt.show()
