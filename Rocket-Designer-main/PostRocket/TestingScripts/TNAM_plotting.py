#Author: Alex Post, Jasper Stedman
#Date: 07/28/2024

import numpy as np
from numpy import genfromtxt
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, AutoMinorLocator
from scipy.optimize import curve_fit

#Extract the data from the following csv files
errorlist = np.ndarray([500,1])
errorlist[:,0] = np.logspace(-6,0,500) # errorlist needs to be an ndarray for the plotting functions to work properly

MSE_alt_avg = genfromtxt('MSE_alt_avg.csv', delimiter=',')
MSE_vel_avg = genfromtxt('MSE_vel_avg.csv', delimiter=',')
MSE_acc_avg = genfromtxt('MSE_acc_avg.csv', delimiter=',')
Runtime_avg = genfromtxt('Runtime_avg.csv', delimiter=',')

MSE_alt_dev = genfromtxt('MSE_alt_dev.csv', delimiter=',')
MSE_vel_dev = genfromtxt('MSE_alt_dev.csv', delimiter=',')
MSE_acc_dev = genfromtxt('MSE_acc_dev.csv', delimiter=',')
Runtime_dev = genfromtxt('Runtime_dev.csv', delimiter=',')

for j in range(3):
    for i,point in enumerate(MSE_acc_avg[:,j]):
        if point > 100 * np.mean(MSE_acc_avg[:,j]):
         MSE_acc_avg[i,j] = MSE_acc_avg[i-1,j]
         MSE_vel_avg[i,j] = MSE_vel_avg[i-1,j]
         MSE_alt_avg[i,j] = MSE_alt_avg[i-1,j]
         Runtime_avg[i,j] = Runtime_avg[i-1,j]

         MSE_acc_dev[i,j] = MSE_acc_dev[i-1,j]
         MSE_vel_dev[i,j] = MSE_vel_dev[i-1,j]
         MSE_alt_dev[i,j] = MSE_alt_dev[i-1,j]
         Runtime_dev[i,j] = Runtime_dev[i-1,j]

RMSE_vel_avg = np.sqrt(MSE_vel_avg)
RMSE_alt_avg = np.sqrt(MSE_alt_avg)
RMSE_acc_avg = np.sqrt(MSE_acc_avg)
RMSE_alt_dev = np.sqrt(MSE_alt_dev)
RMSE_vel_dev = np.sqrt(MSE_vel_dev)
RMSE_acc_dev = np.sqrt(MSE_acc_dev)


plotColors = ['orange', 'green', 'blue']
plotLabels = ['RK45' , 'RK23' , 'LSODA']

# def logLine(x,a,n):
#    y = a*x**n
#    return y
    
# def logSin(x,a,b,c,d):
#    y = 10**( a*np.sin(b*np.log10(x) + c) + d )
#    return y
#         # fitData = curve_fit(logLine, xData, yData[:,i])
#         # ax.plot(xData,logLine(xData,*fitData[0]), color= colors[i], label=labels[i], alpha = 0.9)

# # Plotting

def logLogPolyfit(xData,yData):
   
    xLog = np.log10(xData)
    yLog = np.log10(yData)

    fitCoeffs = np.polyfit(xLog,yLog,2)
    yFitLog = np.poly1d(fitCoeffs)
    yFit = 10**yFitLog(xLog)
    
    return yFit

def logLogScatter(xData,yData,yDev,colors,labels):
   
    fig, ax = plt.subplots()
    m0 = np.shape(xData)[-1]
    m1 = np.shape(yData)[-1]

    bothPlotting = m0 == m1

    for j in range(m1):
        if bothPlotting:
            i = j
        else:
            i = 0

        yFit = logLogPolyfit(xData[:,i],yData[:,j])
        ax.scatter(xData[:,i],yData[:,j], 3, colors[j], label=labels[j] + ' data', alpha = 0.4) # Plot Data Scatter Plot
        ax.plot(xData[:,i],yFit, color= colors[j], label=labels[j] + ' fit', alpha = 0.9) # Logarithmic polynomial fit bc curve_fit with ax^n wasn't working
        ax.legend(fontsize=12, loc='best')

        ax.fill_between(xData[:,i], yFit + 2*yDev[:,j], yFit - 2*yDev[:,j], color=colors[j], alpha=0.25) # Standard Deviation Plotting

        ax.set_xscale("log", base=10)
        ax.set_yscale("log", base=10)
        ax.grid(True,'major')
        ax.grid(True,'minor', alpha= 0.25)
   
    return fig, ax

def logLinScatter(xData,yData,yDev,colors,labels,*args):
   
    fig, ax = plt.subplots()
    m0 = np.shape(xData)[-1]
    m1 = np.shape(yData)[-1]

    bothPlotting = m0 == m1

    for j in range(m1):
        if bothPlotting:
            i = j
        else:
            i = 0

        yFit = logLogPolyfit(xData[:,i],yData[:,j])
        ax.scatter(xData[:,i],yData[:,j], 3, colors[j], label=labels[j] + ' data', alpha = 0.4) # Plot Data Scatter Plot
        ax.plot(xData[:,i],yFit, color= colors[j], label=labels[j] + ' fit', alpha = 0.9) # Logarithmic polynomial fit bc curve_fit with ax^n wasn't working
        ax.legend(fontsize=12, loc='best')

        ax.fill_between(xData[:,i], yFit + 2*yDev[:,j], yFit - 2*yDev[:,j], color=colors[j], alpha=0.25) # Standard Deviation Plotting

        ax.set_xscale("log", base=10)
        ax.grid(True,'major')
        ax.grid(True,'minor', alpha= 0.25)
   
    return fig, ax



# MSE Alt
fig0, ax0 = logLogScatter(errorlist, RMSE_alt_avg, RMSE_alt_dev, plotColors, plotLabels)
ax0.set_title('Root Mean Square Error Altitude', fontsize=20)
ax0.set_xlabel('Tolerance', fontsize=15)
ax0.set_ylabel('Root Mean Square Error (m)', fontsize=15)


# #MSE Vel
fig1, ax1 = logLogScatter(errorlist, RMSE_vel_avg, RMSE_vel_dev, plotColors, plotLabels)

ax1.set_title('Root Mean Square Error Velocity', fontsize=20)
ax1.set_xlabel('Tolerance', fontsize=15)
ax1.set_ylabel('Root Mean Square Error (m/s)', fontsize=15)



#MSE Acc

fig2, ax2 = logLogScatter(errorlist, RMSE_acc_avg, RMSE_acc_dev, plotColors, plotLabels)
ax2.set_title('Root Mean Square Error Acceleration', fontsize=20)
ax2.set_xlabel('Tolerance', fontsize=15)
ax2.set_ylabel('Root Mean Square Error (m/s^2)', fontsize=15)

# #Runtime with standard deviation

fig3, ax3 = logLinScatter(errorlist, Runtime_avg, Runtime_dev, plotColors, plotLabels)
ax3.set_title('Average Runtime', fontsize=20)
ax3.set_xlabel('Tolerance', fontsize=15)
ax3.set_ylabel('Runtime (s)', fontsize=15)

#MSE (x) Runtime (y)
RMSE_alt_avg_sort = RMSE_alt_avg
Runtime_avg_sort = Runtime_avg
Runtime_dev_sort = Runtime_dev
for i in range(3):
    sortedInds = np.argsort(RMSE_alt_avg[:,i])
    RMSE_alt_avg_sort[:,i] = RMSE_alt_avg[:,i][sortedInds]
    Runtime_avg_sort[:,i] = Runtime_avg[:,i][sortedInds]
    Runtime_dev_sort[:,i] = Runtime_dev[:,i][sortedInds]

fig4, ax4 = logLinScatter(RMSE_alt_avg_sort, Runtime_avg_sort, Runtime_dev_sort, plotColors, plotLabels)
ax4.set_title('Runtime vs Root Mean Square Error', fontsize=20)
ax4.set_xlabel('Root Mean Square Error (m)', fontsize=15)
ax4.set_ylabel('Runtime (s)', fontsize=15)


plt.show()


