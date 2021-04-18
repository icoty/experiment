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

plt.figure(1) # 创建图表1
plt.figure(2) # 创建图表2
#x1 = plt.subplot(211) # 在图表2中创建子图1
#x2 = plt.subplot(212) # 在图表2中创建子图2

# https://blog.csdn.net/qq_32867925/article/details/103270531?utm_medium=distribute.pc_relevant.none-task-blog-BlogCommendFromMachineLearnPai2-2.control&dist_request_id=&depth_1-utm_source=distribute.pc_relevant.none-task-blog-BlogCommendFromMachineLearnPai2-2.control

file_name = sys.argv[1]
jpg_file = file_name.replace(".json", ".jpg")

dump_data = json.load(open(file_name))
with open(file_name, "r") as ouf:
    #time1 = [x/1000 for x in dump_data["ABC"]['time']]
    #time2 = [x/1000 for x in dump_data["EABC"]['time']]
    time3 = [x/1000 for x in dump_data["HABC"]['time']]
    time4 = [x/1000 for x in dump_data["OABC"]['time']]
    time5 = [x/1000 for x in dump_data["mABC"]['time']]
    time6 = [x/1000 for x in dump_data["MeanABC"]['time']]
    time7 = [x/1000 for x in dump_data["NABC"]['time']]

    #fitness1 = dump_data["ABC"]['fitness']
    #fitness2 = dump_data["EABC"]['fitness']
    fitness3 = dump_data["HABC"]['fitness']
    fitness4 = dump_data["OABC"]['fitness']
    fitness5 = dump_data["mABC"]['fitness']
    fitness6 = dump_data["MeanABC"]['fitness']
    fitness7 = dump_data["NABC"]['fitness']

    xt = [i for i in range(0, len(time3))]

    plt.figure(1) # 创建图表1
    #plt.plot(xt,time1,',',color = 'k',label="ABC")
    #plt.plot(xt,time2,'-.',color = 'g',label="EABC")
    plt.plot(xt,time3,'--',color = 'r',label="HABC")
    plt.plot(xt,time4,',',color = 'b',label="OABC")
    plt.plot(xt,time5,':',color = 'm',label="mABC")
    plt.plot(xt,time6,'--',color = 'c',label="MeanABC")
    plt.plot(xt,time7,'-.',color = 'g',label="NABC")
    if "SFLA" in dump_data:
        time8 = [x/1000 for x in dump_data["SFLA"]['time']]
        #plt.plot(xt,time8,':',color = 'y',label="SFLA")

    plt.xlabel(u"迭代次数") #横坐标名字
    plt.ylabel(u"累计时间(毫妙)")  #纵坐标名字
    plt.legend(loc = "best")    #图例

    xf = [i for i in range(0, len(fitness3))]
    plt.figure(2) # 创建图表2
    #plt.plot(xf,fitness1,',',color = 'k',label="ABC")
    #plt.plot(xf,fitness2,'-.',color = 'g',label="EABC")
    plt.plot(xf,fitness3,'--',color = 'r',label="HABC")
    plt.plot(xf,fitness4,',',color = 'b',label="OABC")
    plt.plot(xf,fitness5,':',color = 'm',label="mABC")
    plt.plot(xf,fitness6,'--',color = 'c',label="MeanABC")
    plt.plot(xf,fitness7,'-.',color = 'g',label="NABC")
    if "SFLA" in dump_data:
        fitness8 = dump_data["SFLA"]['fitness']
        plt.plot(xf,fitness8,':',color = 'y',label="SFLA")

    plt.xlabel(u"迭代次数") #横坐标名字
    plt.ylabel(u"适应度值")  #纵坐标名字
    plt.legend(loc = "best")    #图例

    plt.show()
