import df_to_force
import sys
import numpy as np
import matplotlib.pyplot as plt

folder_root = "C:\\Users\\bendrevniok\\Dropbox\\"
folder_root += "DrevniokShareWithDrevniok\\benzene\\df_data\\"
folder_root += "30June2011\\"

A = 2.5e-10
k = 2800
f0 = 56400.0
fs = 18.3
order = 1
freq = 0.1

sm_files = ["df_90.txt", "df_91.txt","df_98.txt", "df_107.txt","df_108.txt","df_113.txt","df_114.txt" ]
ad_files = ["df_86.txt","df_87.txt","df_97.txt","df_109.txt","df_110.txt"]
ch_files = ["df_88.txt","df_89.txt","df_95.txt","df_96.txt","df_111.txt","df_112.txt"]

sm = []
ad = []
ch = []

for file_name in sm_files:
    data = df_to_force.ReadOmicronXY(folder_root, "sm\\"+file_name)
    data.load_and_get()
    sm.append(df_to_force.SmallMolecule(data.z_approach, data.df_approach, 
                                        A, k, f0, fs))
    #sm[i].filtfilt(order, freq)

for file_name in ad_files[3:]:
    data = df_to_force.ReadOmicronXY(folder_root, "ad\\"+file_name)
    data.load_and_get()
    ad.append(df_to_force.Adatom(data.z_approach, data.df_approach, 
                                A, k, f0, fs))
    #sm[i].filtfilt(order, freq)

for file_name in ch_files:
    data = df_to_force.ReadOmicronXY(folder_root, "ch\\"+file_name)
    data.load_and_get()
    ch.append(df_to_force.CornerHole(data.z_approach, data.df_approach, 
                                    A, k, f0, fs))
    
    
for i in range(len(sm)):
    if i < 3:
        sm[i].remove_long_range(ch[1])
    else: 
        sm[i].remove_long_range(ch[4])
    sm[i].filtfilt(order, freq)
    sm[i].sader_jarvis()
    #sm[i].make_plot("sm", sm_files[i][3:])
    sm[i].normalize()
    sm[i].make_plot("sm_norm", sm_files[i][3:])
#for i in range(len(ad)):
#    ad[i].remove_long_range(ch[4])
#    ad[i].filtfilt(order, freq)
#    print ad[i].k
#    ad[i].sader_jarvis()
#    ad[i].make_plot()
#    
#ad_av = df_to_force.AverageCurve(sm)
#ad_av.make_average_f()
#ad_av.make_plot()

sm_av = df_to_force.AverageCurve(sm)
sm_av.find_minimum()
sm_av.make_average_f()
sm_av.make_plot()



#for i in range(len(ad)):
#    ad[i].remove_long_range(ch[0])
#    ad[i].sader_jarvis()
#    ad[i].make_plot(ad[i].legend, ad_files[i][3:])
    
#for i in range(len(ch)):    
#    ch[i].make_plot(ch[i].legend, ch_files[i][3:])
#    


#sma = df_to_force.AverageCurve(sm)
#sma.average()
#sma.make_plot()
plt.legend()
plt.show()
