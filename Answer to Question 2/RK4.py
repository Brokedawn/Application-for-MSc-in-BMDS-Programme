import numpy as np
import math
import matplotlib.pyplot as plt
import pandas as pd
#dE/dt=-k1[E][S]+k2[ES]+k3[ES]
#dS/dt=k2[ES]-k1[E][S]
#dES/dt=- dE/dt=k1[E][S]-k2[ES]-k3[ES]
#dP/dt=k3[ES]
#k1=100/µM/min, k2=600/min, k3=150/min
#y1:Et0=1µM, y2:St0=10µM,y3:ESt0=0, y4:Pt0=0
def f1(x, y1, y2,y3):
    df = -100*y1*y2+600*y3+150*y3
    return df
 
def f2(x, y1, y2,y3):
    df = 600*y3-100*y1*y2
    return df
 
def f3(x,y1,y2,y3):
    df=100*y1*y2-600*y3-150*y3
    return df
 
def f4(x,y1,y2,y3) :
    df=150*y3
    return df

def RK4(x, y1, y2,y3,y4, h):
    """
    :param x: Initial value of X
    :param y1: Initial value of y1
    :param y2: Initial value of y2
    :param y3: Initial value of y3
    :param y4: Initial value of y4
    :param h: time step
    :return: New iterative solution
    """
    xarray, y1array, y2array,y3array,y4array= [], [], [], [],[]
    while x <= 1:
        xarray.append(x)
        y1array.append(y1)
        y2array.append(y2)
        y3array.append(y3)
        y4array.append(y4)
        x += h
 
        K_1 = f1(x, y1, y2,y3)
        L_1 = f2(x, y1, y2,y3)
        M_1 = f3(x, y1, y2,y3)
        N_1 = f4(x,y1,y2,y3)
        K_2 = f1(x + h / 2, y1 + h / 2 * K_1, y2 + h / 2 * L_1 , y3 + h/2 * M_1)
        L_2 = f2(x + h / 2, y1 + h / 2 * K_1, y2 + h / 2 * L_1 , y3 + h/2 * M_1)
        M_2 = f3(x + h / 2, y1 + h / 2 * K_1, y2 + h / 2 * L_1, y3 + h / 2 * M_1)
        N_2 = f4(x + h / 2, y1 + h / 2 * K_1, y2 + h / 2 * L_1, y3 + h / 2 * M_1)
        K_3 = f1(x + h / 2, y1 + h / 2 * K_2, y2 + h / 2 * L_2 , y3 + h/2 * M_2)
        L_3 = f2(x + h / 2, y1 + h / 2 * K_2, y2 + h / 2 * L_2 , y3 + h/2 * M_2)
        M_3 = f3(x + h / 2, y1 + h / 2 * K_2, y2 + h / 2 * L_2, y3 + h / 2 * M_2)
        N_3 = f4(x + h / 2, y1 + h / 2 * K_2, y2 + h / 2 * L_2, y3 + h / 2 * M_2)
        K_4 = f1(x + h, y1 + h * K_3, y2 + h * L_3, y3 + h * M_3)
        L_4 = f2(x + h, y1 + h * K_3, y2 + h * L_3, y3 + h * M_3)
        M_4 = f3(x + h, y1 + h * K_3, y2 + h * L_3, y3 + h * M_3)
        N_4 = f4(x + h, y1 + h * K_3, y2 + h * L_3, y3 + h * M_3)
        y1 = y1 + (K_1 + 2 * K_2 + 2 * K_3 + K_4) * h / 6
        y2 = y2 + (L_1 + 2 * L_2 + 2 * L_3 + L_4) * h / 6
        y3 = y3 + (M_1 + 2 * M_2 + 2 * M_3 + M_4) * h / 6
        y4 = y4 + (N_1 + 2 * N_2 + 2 * N_3 + N_4) * h / 6
    return xarray, y1array, y2array,y3array,y4array
 
def getrateofchange(list1):
    fund_arr = np.array(list1)
    array = np.diff(fund_arr) / fund_arr[:-1]
   # plt.plot(array)
    return array.tolist()
    
def main():
    #2.2
    xarray, y1array, y2array,y3array,y4array = RK4(0,1,10,0,0,0.001)
    plt.figure('Runge Kutta numerical results')
    plt.subplot(221)
    
    plt.scatter(xarray, y1array, label='E_scatter', s=1, c='#DC143C', alpha=0.6)
    
    plt.legend()
    plt.subplot(222)
    #plt.plot(xarray, y2array, label='y2_runge_kutta')
    plt.scatter(xarray, y2array, label='S_scatter', s=1, c='#DC143C', alpha=0.6)
    #plt.y2label('x')
    plt.legend()
    plt.subplot(223)
    #plt.plot(xarray, y3array, label='y3_runge_kutta')
    plt.scatter(xarray, y3array, label='ES_scatter', s=1, c='#DC143C', alpha=0.6)
    #plt.y3label('x')
    plt.legend()
    plt.subplot(224)
    plt.scatter(xarray, y4array, label='P_scatter', s=1, c='#DC143C', alpha=0.6)
    plt.legend()
    plt.show()
    
    #2.3
    y2array.remove(y2array[0])
    v=getrateofchange(y4array)
    
    plt.plot(y2array,v)
    
    
if __name__ == '__main__':
    main()