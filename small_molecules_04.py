import df_to_force
import sys
import numpy as np
import matplotlib.pyplot as plt

plt.close("all")
folder_list = ["C:\\Users\\bendrevniok\\Dropbox\\",
        "DrevniokShareWithDrevniok\\benzene\\df_data\\","30June2011\\"]
folder_root = folder_list[0]+folder_list[1]+folder_list[2]

A = 2.5e-10
k = 3500
f0 = 56400.0
fs = 18.3
order = 1
freq = 0.1

sm_files = ["df_90.txt", "df_91.txt", "df_107.txt","df_108.txt","df_113.txt"]
sm_files = ["df_90.txt", "df_91.txt"]
sm_files = [ "df_107.txt","df_108.txt"]
#sm_files = ["df_90.txt", "df_91.txt","df_98.txt", "df_107.txt",
        #"df_108.txt","df_113.txt","df_114.txt" ]
sm_files = ["df_98.txt", "df_107.txt",
        "df_108.txt","df_113.txt","df_114.txt"]


ad_files = ["df_97.txt","df_109.txt","df_110.txt"]
#ad_files = ["df_86.txt","df_87.txt","df_109.txt","df_110.txt"]
#ad_files = ["df_86.txt","df_87.txt"]


#ch_files = ["df_88.txt","df_89.txt","df_95.txt","df_96.txt","df_111.txt","df_112.txt"]
ch_files = ["df_95.txt","df_96.txt","df_111.txt","df_112.txt"]

sm = []
ad = []
ch = []

for file_name in sm_files:
    data = df_to_force.ReadOmicronXY(folder_root, "sm\\"+file_name)
    data.load_and_get()
    sm.append(df_to_force.SmallMolecule(data.z_approach, data.df_approach, 
                                        A, k, f0, fs,data.name))

for file_name in ad_files:
    data = df_to_force.ReadOmicronXY(folder_root, "ad\\"+file_name)
    data.load_and_get()
    ad.append(df_to_force.Adatom(data.z_approach, data.df_approach, 
                                A, k, f0, fs,data.name))

for file_name in ch_files:
    data = df_to_force.ReadOmicronXY(folder_root, "ch\\"+file_name)
    data.load_and_get()
    ch.append(df_to_force.CornerHole(data.z_approach, data.df_approach, 
                                    A, k, f0, fs,data.name))
    
    
for i in range(len(sm)):
    sm[i].shift_to_zero()
    
for i in range(len(ch)):
    ch[i].shift_to_zero()
            
for i in range(len(ad)):
    ad[i].shift_to_zero()


ad_av = df_to_force.AverageCurve(ad)
ch_av = df_to_force.AverageCurve(ch)
sm_av = df_to_force.AverageCurve(sm)


ad_av.make_average_df()
ch_av.make_average_df()
sm_av.make_average_df()


chav = df_to_force.CornerHole(ch_av.force_curves[0].z, ch_av.average_df, 
        A, k, f0,1.0,data.name)
smav = df_to_force.SmallMolecule(sm_av.force_curves[0].z, sm_av.average_df, 
        A, k, f0,1.0,data.name)
adav = df_to_force.Adatom(ad_av.force_curves[0].z, ad_av.average_df, 
        A, k, f0,1.0,data.name)
        
chav.filtfilt(order, freq)
smav.filtfilt(order, freq)
adav.filtfilt(order, freq)

smav.shift_to_zero()
adav.shift_to_zero()
chav.shift_to_zero()

smav.sader_jarvis()
adav.sader_jarvis()
chav.sader_jarvis()

#plt.plot((smav.distance)*1e9, (smav.force-chav.force)*1e9,
#                    linewidth = 3.0, 
#                    label = smav.legend+" "+smav.name+"f")
#                    
#plt.plot((adav.distance)*1e9, (adav.force-chav.force)*1e9,
#                    linewidth = 3.0, 
#                    label = adav.legend+" "+adav.name+"f")

#chav.make_plot([-0.2,1.65],[-3.5,0.25])
#adav.make_plot([-0.2,1.65],[-3.5,0.25])
smav.make_plot([-0.2,1.65],[-3.5,0.25])



plt.show()
