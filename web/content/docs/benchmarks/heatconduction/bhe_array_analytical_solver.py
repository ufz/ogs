# -*- coding: utf-8 -*- {}
import matplotlib.pyplot as plt
import numpy as np
import scipy
import os
import math
from mpmath import *
mp.dps = 25; mp.pretty = True

#part 1:gif_splitfile to txt

#新建文件夹txt，用于存储原文件的分解文件集。
os.mkdir("C://george//PhD//papers//BHE_sc//result//txt")
#os.chdir(path) 方法用于改变当前工作目录到指定的路径。此处创建路径到新建文件夹txt
os.chdir("C://george//PhD//papers//BHE_sc//result//txt")
txtfiles = os.getcwd() #获得当前运行脚本所在目录作为文件储存地址

# define the path of the data file
str_tec_file_path = "C://george//PhD//papers//BHE_sc//result//tec_0_5m.tec"

data_num =0
with open(str_tec_file_path, 'r') as f_in:
    i = -1 #此处写-1是因为要使输出文件名第一个必须是以0结尾，不然和part2的文件读入循环对不上。必须要对上是因为\n
           #part2数组值从0位开始储存的。
    for line in f_in:
        if '=' not in line:
            txtfiles.write(line)
            #读出单位文件中有多少行数据，存入变量data_num
            data_num = data_num + 1
        elif 'ZONE' in line: #从关键词ZONE开始读取行
            data_num = 0
            i = i + 1
            txtfiles = open('txtfile{}.txt'.format(i), 'w')#格式化命名文件名存储

#记录总共的分出文件数目，用于part 2
numfile = i + 1
#关闭并完整文件释放内存
txtfiles.close()


# plotting T_in T_out curves,建立作图文件夹及及改变工作路径
os.mkdir("C://george//PhD//papers//BHE_sc//result//png")
os.chdir("C://george//PhD//papers//BHE_sc//result//png")
pngfiles = os.getcwd()

#建立source term坐标矩阵，用以计算各热源对referece point的温度影响矩阵
po_x = np.array([40,40,40,40,40],dtype=float).reshape(-1,1) + np.arange(0,25,5)
po_y = np.arange(60,35,-5).reshape(-1,1) + np.zeros(5)

po_dist_to_referencepo = np.zeros([5,5])
Temp_po_to_referencepo = np.zeros([5,5])
#解析解1基础数据项,q是timestep变化量
q1 = -35
q2 = 35
q3 = 0

time_trans = 120*24*60*60 #3个月一输出数据
T0 = 10
T = T0


poro = 10e-20
lamda_f = 0.5984
density_f = 998.2032
cp_f = 4182.0
lamda_sp = 2.0
density_sp = 1950.0
cp_sp = 1500.0

alpha_f = lamda_f/(density_f*cp_f)
alpha_sp = lamda_sp/(density_sp*cp_sp)

lamda_med = (1 - poro)*lamda_sp + poro*lamda_f
alpha_med = (1 - poro)*alpha_sp + poro*alpha_f

#解析解2基础数据项：
qq = np.array([-35,0,0,-35,0,0,-35,0,0,-35,0]).reshape(1,-1)
qq_all = np.repeat(qq,data_num,axis=0)
expandedloads= np.array(qq)
numtimesteps = 11

numbhe = 25
ppo_x = np.array([40,40,40,40,40],dtype=float).reshape(-1,1) + np.arange(0,25,5)
ppo_y = np.arange(60,35,-5).reshape(-1,1) + np.zeros(5)
ppo_x_re = np.reshape(po_x, (-1,1))
ppo_y_re = np.reshape(po_y, (-1,1))

# 建立动态变量名
createVar = locals()

# 建立月份名字表，因为数据记录从1月1日有一组数据，所以多出一个月数据。
month = [" Temperature after 0 months ", " Temperature after 4 months", " Temperature after 8 months", " Temperature after 12 months", " Temperature after 16 months", " Temperature after 20 months",
          " Temperature after 24 months", " Temperature after 28 months", " Temperature after 32 months", " Temperature after 36 months"," Temperature after 40 months"," Temperature after 44 months"]

#解析解2项：
txtpath2 = f"C://george//PhD//papers//BHE_sc//result//txt//txtfile0.txt"#f格式化字符串
dist_arrayy = 0.0
with open(txtpath2, 'r') as f_in_txt:
    for line in f_in_txt:
        # split line
        data_line = line.split(" ") #maxsplit
        # time_double = 0.0
        dist_double2 = float(data_line[0])
        cordinate2 = math.sqrt(dist_double2**2/2)

        #  put into the list
        dist_arrayy= np.vstack((dist_arrayy, cordinate2))

    dist_arrayy_new = np.delete(dist_arrayy, 0, 0)

    numtemppoints = len(dist_arrayy_new)
    T2=np.zeros([numtemppoints,numtimesteps])

    coeff_all = np.zeros([numtemppoints,numtimesteps])

    for currstep in range(0,numtimesteps):
        Temp_po_to_referencepo= np.zeros([numtemppoints,numbhe])
        po_dist_to_referencepo= np.zeros([numtemppoints,numbhe])
        localcoeff_all= np.zeros([numtemppoints,1])
        localcoeff= np.zeros([numtemppoints,numbhe])
        localcoeff1= np.zeros([numtemppoints,numbhe])
        for i in range(0,numbhe):
            if(time_trans*(currstep+1)-time_trans*0>0):
                for j in range(0,numtemppoints):
                    po_dist_to_referencepo[j,i] = abs(ppo_x_re[i] - (dist_arrayy_new[j]+1))**2 + abs(ppo_y_re[i] - dist_arrayy_new[j])**2
                    exp = po_dist_to_referencepo[j,i]/(4*alpha_med*time_trans*(currstep+1))
                    n = e1(exp)
                    localcoeff[j,i] = 1/(4*math.pi*lamda_med)*n
            if(time_trans*(currstep+1)-time_trans*1>0):
                for j in range(0,numtemppoints):
                    po_dist_to_referencepo[j,i] = abs(ppo_x_re[i] - (dist_arrayy_new[j]+1))**2 + abs(ppo_y_re[i] - dist_arrayy_new[j])**2
                    exp1 = po_dist_to_referencepo[j,i]/(4*alpha_med*time_trans*currstep)
                    n1 = e1(exp1)
                    localcoeff[j,i] = localcoeff[j,i] - 1/(4*math.pi*lamda_med)*n1

        localcoeff_all= np.sum(localcoeff,axis=1).reshape(-1,1)
        coeff_all[:,1:]=coeff_all[:,:numtimesteps-1]
        coeff_all[:,:1]=localcoeff_all

    for currstep in range(0,numtimesteps):
        T2[:,currstep] = np.sum(coeff_all[:,numtimesteps-1-currstep:]*qq_all[:,:currstep+1],axis=1) +10



# loop over
for m in range(0 , numfile):
    #文件读取路径
    txtpath      = f"C://george//PhD//papers//BHE_sc//result//txt//txtfile{m}.txt"#f格式化字符串
    txtpath_ogs6 = f"C://george//PhD//papers//BHE_sc//result//txt_ogs6//txtfile{m}.txt"

    #动态变量名每组赋初值0.0
    createVar['dist_array'+ str(m)] = 0.0
    createVar['temperature_array_ogs5'+ str(m)] = 0.0
    createVar['temperature_array_ana'+ str(m)] = 0.0
    createVar['dist_array2'+ str(m)] = 0.0
    createVar['temperature_array_ana2'+ str(m)] = 0.0

    createVar['dist_array_ogs6'+ str(m)]             =  0.0
    createVar['temperature_array_ogs6'+ str(m)]      =  0.0

    #解析解1：
#    if m == 0 or m ==2 or m ==3 or m ==5 or m ==6 or m ==8 or m ==9 or m ==10:
#        q = q2
#    elif m == 1 or m ==4 or m ==7:
#        q = q1
    if m == 0:
        with open(txtpath, 'r') as f_in_txt:
            for line in f_in_txt:
                # split line
                data_line = line.split(" ") #maxsplit
                # time_double = 0.0
                dist_double = float(data_line[0])
                cordinate = math.sqrt(dist_double**2/2)
                temperature = float(data_line[1])

                #解析解数据项
                time = time_trans + 10e-6
                for i in range(0,5):
                    for j in range(0,5):
                        po_dist_to_referencepo[i,j] = abs(po_x[i,j]-(cordinate+1))**2 + abs( po_y[i,j]-cordinate)**2
                for i in range(0,5):
                    for j in range(0,5):
                        exp = po_dist_to_referencepo[i,j]/(4*alpha_med*time)
                        n1 = e1(exp)
                        Temp_po_to_referencepo[i,j] = q3/(4*math.pi*lamda_med)*n1

                T = np.sum(Temp_po_to_referencepo) + T0

                #  put into the list
                createVar['dist_array'+ str(m)] = np.vstack((createVar['dist_array'+ str(m)], cordinate))
                createVar['temperature_array_ogs5'+ str(m)] = np.vstack((createVar['temperature_array_ogs5'+ str(m)], temperature))
                createVar['temperature_array_ana'+ str(m)] = np.vstack((createVar['temperature_array_ana'+ str(m)], T))

        f_in_txt.close()
        # end of loop

        #去除21和22行的0.0值，因为此值被带入了dist_array的list第一行，可查看variable explorer temperature_array观察。
        #目的是去除做图时的第一个数据0值。
        dist_array_new = np.delete(createVar['dist_array'+ str(m)], 0, 0)
        temperature_array_ogs5_new = np.delete(createVar['temperature_array_ogs5'+ str(m)], 0, 0)
        temperature_array_ana_new = np.delete(createVar['temperature_array_ana'+ str(m)], 0, 0)

        with open(txtpath_ogs6, 'r') as f_in_txt:
            for line in f_in_txt:
                # split line
                data_line = line.split(" ") #maxsplit
                # time_double = 0.0
                dist_double = float(data_line[0])
                cordinate = math.sqrt(dist_double**2/2)
                temperature = float(data_line[1])
                createVar['dist_array_ogs6'+ str(m)] = np.vstack((createVar['dist_array_ogs6'+ str(m)], cordinate))
                createVar['temperature_array_ogs6'+ str(m)] = np.vstack((createVar['temperature_array_ogs6'+ str(m)], temperature))
        f_in_txt.close()
        dist_array_new_ogs6 = np.delete(createVar['dist_array_ogs6'+ str(m)], 0, 0)
        temperature_array_ogs6_new = np.delete(createVar['temperature_array_ogs6'+ str(m)], 0, 0)

        #plotting
        plt.figure()
        plt.plot(scipy.dot(dist_array_new,1.0),scipy.dot(temperature_array_ogs5_new,1.0),color = 'r',ls=':',lw=1, marker='^',markersize=1.5,  label= 'OGS5')
        plt.plot(scipy.dot(dist_array_new_ogs6,1.0),scipy.dot(temperature_array_ogs6_new,1.0),c='g', ls=':',lw=1, marker='o',markersize=1,  label= 'OGS6')
#        plt.plot(scipy.dot(dist_array_new,1.0),scipy.dot(temperature_array_ana_new,1.0), color = 'g',label= 'ana1' + months[m])
        plt.plot(scipy.dot(dist_array_new,1.0),scipy.dot(temperature_array_ana_new,1.0), "b",label= 'Analytical')
        plt.xlim([0,100])
        plt.ylim([-10,20])
        plt.ylabel('Temperature [$^\circ$C]')
        plt.xlabel('x [m]')
        plt.legend(loc='best',fontsize=8)
        plt.title(month[m],fontsize=12)
        # plt.show()
        plt.savefig('pngfile{}.png'.format(m),dpi = 300, transparent = False)


    else:
        #解析解数据项,实际是温度场斜三角矩阵相加
        T = np.zeros([12])
        with open(txtpath, 'r') as f_in_txt:
            for line in f_in_txt:
                # split line
                data_line = line.split(" ") #maxsplit
                # time_double = 0.0
                dist_double = float(data_line[0])
                cordinate = math.sqrt(dist_double**2/2)

                i1 = m +1
                for m1 in range(1, m + 1 ):
                    i1 = i1 - 1
                    #时间变量，timestep温度随时间叠加矩阵每行重新设为0
                    time = time_trans * i1 + 10e-6
                    if  m1 == 2 or m1 == 5 or m1 == 8 or m1 == 11:
                        q = q2
                    elif m1 == 1 or m1 ==4 or m1 ==7 or m1 ==10:
                        q = q1
                    elif m1 == 3 or m1 == 6 or m1 == 9:
                        q = q3

                    for i in range(0,5):
                        for j in range(0,5):
                            po_dist_to_referencepo[i,j] = abs(po_x[i,j]-(cordinate+1))**2 + abs( po_y[i,j]-cordinate)**2
                    for i in range(0,5):
                        for j in range(0,5):
                            exp = po_dist_to_referencepo[i,j]/(4*alpha_med*time)
                            n1 = e1(exp)
                            Temp_po_to_referencepo[i,j] = q/(4*math.pi*lamda_med)*n1

                    T[m1] = np.sum(Temp_po_to_referencepo)

                T_sum = np.sum(T) + T0

                createVar['temperature_array_ana'+ str(m)] = np.vstack((createVar['temperature_array_ana'+ str(m)], T_sum))

        f_in_txt.close()


        #ogs5项：
        with open(txtpath, 'r') as f_in_txt:
            for line in f_in_txt:
                # split line
                data_line = line.split(" ") #maxsplit
                # time_double = 0.0
                dist_double = float(data_line[0])
                cordinate = math.sqrt(dist_double**2/2)
                temperature = float(data_line[1])

                #  put into the list
                createVar['dist_array'+ str(m)] = np.vstack((createVar['dist_array'+ str(m)], cordinate))
                createVar['temperature_array_ogs5'+ str(m)] = np.vstack((createVar['temperature_array_ogs5'+ str(m)], temperature))

        f_in_txt.close()
        #去除21和22行的0.0值，因为此值被带入了dist_array的list第一行，可查看variable explorer temperature_array观察。
        #目的是去除做图时的第一个数据0值。
        dist_array_new = np.delete(createVar['dist_array'+ str(m)], 0, 0)
        temperature_array_ogs5_new = np.delete(createVar['temperature_array_ogs5'+ str(m)], 0, 0)
        temperature_array_ana_new = np.delete(createVar['temperature_array_ana'+ str(m)], 0, 0)


        #ogs6项
        with open(txtpath_ogs6, 'r') as f_in_txt:
            for line in f_in_txt:
                # split line
                data_line = line.split(" ") #maxsplit
                # time_double = 0.0
                dist_double = float(data_line[0])
                cordinate = math.sqrt(dist_double**2/2)
                temperature = float(data_line[1])

                #  put into the list
                createVar['dist_array_ogs6'+ str(m)] = np.vstack((createVar['dist_array_ogs6'+ str(m)], cordinate))
                createVar['temperature_array_ogs6'+ str(m)] = np.vstack((createVar['temperature_array_ogs6'+ str(m)], temperature))

        f_in_txt.close()
        dist_array_new_ogs6 = np.delete(createVar['dist_array_ogs6'+ str(m)], 0, 0)
        temperature_array_ogs6_new = np.delete(createVar['temperature_array_ogs6'+ str(m)], 0, 0)

        #plotting
        plt.figure()
        plt.plot(scipy.dot(dist_array_new,1.0),scipy.dot(temperature_array_ogs5_new,1.0),color = 'r',ls=':',lw=1, marker='^',markersize=1.5,  label= 'OGS5')
        plt.plot(scipy.dot(dist_array_new_ogs6,1.0),scipy.dot(temperature_array_ogs6_new,1.0),c='g', ls=':',lw=1, marker='o',markersize=1,  label= 'OGS6')

#        plt.plot(scipy.dot(dist_array_new,1.0),scipy.dot(temperature_array_ana_new,1.0), color = 'g',label= 'ana' + month[m])
        plt.plot(scipy.dot(dist_array_new,1.0),T2[:,m-1],"b",label= 'Analytical')
        plt.xlim([0,100])
        plt.ylim([-10,20])
        plt.ylabel('Temperature [$^\circ$C]')
        plt.xlabel('x [m]')
        plt.legend(loc='best',fontsize=8)
        plt.title(month[m],fontsize=12)
        # plt.show()
        plt.savefig('pngfile{}.png'.format(m), dpi = 300, transparent = False)

        # end of loop
txtpath.close()
pngfiles.close()