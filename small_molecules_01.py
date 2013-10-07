import df_to_force
import sys
import numpy as np
import matplotlib.pyplot as plt

plt.close("all")
folder_root = "C:\\Users\\bendrevniok\\Dropbox\\"
folder_root += "DrevniokShareWithDrevniok\\benzene\\df_data\\"
folder_root += "30June2011\\"

A = 2.5e-10
k = 2800
f0 = 56400.0
fs = 18.3
order = 1
freq = 0.2

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
    #sm[i].remove_long_range(ch[4])
    #sm[i].filtfilt(order, freq)
    #sm[i].sader_jarvis()
    #sm[i].make_plot("sm", sm_files[i][3:])
    #sm[i].normalize()
    #sm[i].make_plot("sm_norm", sm_files[i][3:])
    
for i in range(len(ch)):
    ch[i].shift_to_zero()
    #ch[i].remove_long_range(ch[4])
    #ch[i].filtfilt(order, freq)
    #ch[i].sader_jarvis()
    #ch[i].make_plot("ch", ch_files[i][3:])
    #ch[i].normalize()
    #ch[i].make_plot("sm_norm", sm_files[i][3:])
            
for i in range(len(ad)):
    ad[i].shift_to_zero()
    #ad[i].filtfilt(order, freq)
    #ad[i].sader_jarvis()
    #ad[i].make_plot("ad", ad_files[i][3:])
    #ad[i].normalize()
    #ad[i].make_plot("ad_norm", ad_files[i][3:])

#ad_av = df_to_force.AverageCurve(ad)
#ad_av.make_average_f()
#ad_av.make_plot()
ch_av = df_to_force.AverageCurve(ch)
ch_av.make_average_df()
#ch_av.plot_all_df()
ad_av = df_to_force.AverageCurve(ad)
ad_av.make_average_df()
#ad_av.plot_all_df()
sm_av = df_to_force.AverageCurve(sm)
sm_av.make_average_df()
#sm_av.plot_all_df()


chav = df_to_force.CornerHole(ch_av.force_curves[0].z, ch_av.average_df, 
        A, k, f0,1.0,data.name)
smav = df_to_force.SmallMolecule(sm_av.force_curves[0].z, sm_av.average_df, 
        A, k, f0,1.0,data.name)
adav = df_to_force.Adatom(ad_av.force_curves[0].z, ad_av.average_df, 
        A, k, f0,1.0,data.name)
        
chav.filtfilt(order, freq)
smav.filtfilt(order, freq)
adav.filtfilt(order, freq)
smav.remove_long_range(chav)
adav.remove_long_range(chav)
smav.sader_jarvis()
adav.sader_jarvis()
#chav.make_plot()
smav.make_plot()
adav.make_plot()
#sm_av.make_average_f()
#sm_av.make_plot()


plt.show()
