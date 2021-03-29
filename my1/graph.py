#!/usr/bin/env python
# -*- coding: utf-8 -*-
# author: yu.yang1@yitu-inc.com

# https://blog.csdn.net/u012328159/article/details/79240652
# https://blog.csdn.net/Nyte2018/article/details/88835890

import matplotlib
 
matplotlib.use('TkAgg')   #或者PDF, SVG或PS
 
import matplotlib.pyplot as plt

import json
import os

import sys

plt.rcParams['font.sans-serif']=['SimHei'] #用来正常显示中文标签
plt.rcParams['axes.unicode_minus'] = False #用来正常显示负号

#plt.figure(1) # 创建图表1
#plt.figure(2) # 创建图表2
#x1 = plt.subplot(211) # 在图表2中创建子图1
#x2 = plt.subplot(212) # 在图表2中创建子图2

# https://blog.csdn.net/qq_32867925/article/details/103270531?utm_medium=distribute.pc_relevant.none-task-blog-BlogCommendFromMachineLearnPai2-2.control&dist_request_id=&depth_1-utm_source=distribute.pc_relevant.none-task-blog-BlogCommendFromMachineLearnPai2-2.control

file_name = sys.argv[1]
jpg_file = file_name.replace(".json", ".jpg")

dump_data = json.load(open(file_name))
with open(file_name, "r") as ouf:
    '''
    xt = [i for i in range(0, len(dump_data["ABC"]['time']))]
    plt.subplot(122)
    #plt.figure(1) # 创建图表1
    if "ABC" in dump_data:
        time1 = dump_data["ABC"]['time']
        #time1 = [x/10 for x in time1]
        plt.plot(xt,time1,'-',color = 'r',label="ABC")
    if "EABC" in dump_data:
        time2 = dump_data["EABC"]['time']
        #time2 = [x/10 for x in time1]
        plt.plot(xt,time2,'-',color = 'g',label="EABC")
    if "LSABC" in dump_data:
        time3 = dump_data["LSABC"]['time']
        #time3 = [x/10 for x in time1]
        plt.plot(xt,time3,'-.',color = 'y',label="LSABC")
    if "DEABC" in dump_data:
        time4 = dump_data["DEABC"]['time']
        #time4 = [x/10 for x in time1]
        plt.plot(xt,time4,':',color = 'b',label="DEABC")
    if "HSABC" in dump_data:
        time5 = dump_data["HSABC"]['time']
        #time5 = [x/10 for x in time1]
        plt.plot(xt,time5,':',color = 'k',label="HSABC")
    if "IEABC" in dump_data:
        time6 = dump_data["IEABC"]['time']
        #time6 = [x/10 for x in time1]
        plt.plot(xt,time6,'--',color = 'c',label="IEABC")
    if "SFLA" in dump_data:
        time7 = dump_data["SFLA"]['time']
        #time7 = [x/10 for x in time1]
        plt.plot(xt,time7,'--',color = 'm',label="SFLA")
    plt.xlabel(u"迭代次数") #横坐标名字
    plt.ylabel(u"累计时间(微妙)")  #纵坐标名字
    plt.legend(loc = "best")    #图例
    '''

    xf = [i for i in range(0, len(dump_data["ABC"]['fitness']))]
    plt.subplot(111)
    #plt.figure(2)# 创建图表1
    if "ABC" in dump_data:
        fitness1 = dump_data["ABC"]['fitness']
        plt.plot(xf,fitness1,',',color = 'k',label="ABC")

    if "EABC" in dump_data:
        fitness2 = dump_data["EABC"]['fitness']
        plt.plot(xf,fitness2,'-.',color = 'g',label="EABC")

    '''
    if "IEABC" in dump_data:
        fitness6 = dump_data["IEABC"]['fitness']
        plt.plot(xf,fitness6,'--',color = 'c',label="IEABC")
    '''

    if "LSABC" in dump_data:
        fitness3 = dump_data["LSABC"]['fitness']
        plt.plot(xf,fitness3,'--',color = 'y',label="LSABC")

    if "DEABC" in dump_data:
        fitness4 = dump_data["DEABC"]['fitness']
        plt.plot(xf,fitness4,':',color = 'b',label="DEABC")

    if "HSABC" in dump_data:
        fitness5 = dump_data["HSABC"]['fitness']
        plt.plot(xf,fitness5,':',color = 'r',label="HSABC")

    if "SFLA" in dump_data:
        fitness7 = dump_data["SFLA"]['fitness']
        plt.plot(xf,fitness7,'-',color = 'm',label="SFLA")

    #解决中文显示问题
    #plt.rcParams['font.sans-serif']=['FangSong']
    #plt.rcParams['axes.unicode_minus'] = False
    plt.xlabel(u"迭代次数") #横坐标名字
    plt.ylabel(u"适应度值")  #纵坐标名字
    plt.legend(loc = "best")    #图例

    plt.savefig(jpg_file)

    #plt.show()
