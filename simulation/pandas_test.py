from __future__ import print_function
import pandas
import matplotlib.pyplot as plt
# from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import numpy as np

if __name__ == '__main__':
    TB_simulations =  []
    for i in range(1,6):
        TB_simulations.append('E:/PYTHON_PROJECTS/TB_simulation/simulation_statistics_'+str(i)+'0.xlsx')
    # output = file("E:/PYTHON_PROJECTS/TB_simulation/pearson_correlation.txt",'w')
    fig = plt.figure(figsize=(40, 40))
    ax = []
    for i in range (0,5):
        # ax.append(plt.subplot(521+2*i))
        ax.append(fig.add_axes([0.1, 1-0.155*(i+1), 0.7, 0.14]))
    # ymajorLocator = MultipleLocator(10)
    for i in range (0,4):
        data1 = pandas.read_excel(TB_simulations[i])['resistant_mean']
        data2 = pandas.read_excel(TB_simulations[i])['resistant_ratio']
        data3 = pandas.read_excel(TB_simulations[i])['resistant_std']
        plt.sca(ax[i])
        data1.plot(kind='line',color='g')
        data2.plot(kind='line',color='r')
        data3.plot(kind='line',color='b')
        # ax[i].yaxis.set_major_locator(ymajorLocator)
        ax[i].set_ylabel('timespan = '+str(10+10*i),fontsize=10)
        ax[i].yaxis.set_ticks(np.arange(0, 51, 10))
        ax[i].yaxis.grid(True, which='major')
        ax[i].xaxis.set_visible(False)

    parameters = list(pandas.read_excel(TB_simulations[4])['parameters'])
    for i in range (0,len(parameters)):
        parameters[i] = parameters[i].encode("ascii")[3:]
    data1 = pandas.read_excel(TB_simulations[4])['resistant_mean']
    data2 = pandas.read_excel(TB_simulations[4])['resistant_ratio']
    data3 = pandas.read_excel(TB_simulations[4])['resistant_std']
    plt.sca(ax[4])
    data1.plot(kind='line', color='g')
    data2.plot(kind='line', color='r')
    data3.plot(kind='line', color='b')
    ax[4].set_ylabel('timespan = ' + str(50),fontsize=10)
    ax[4].yaxis.grid(True, which='major')
    ax[4].yaxis.set_ticks(np.arange(0, 51, 10))
    ax[4].xaxis.set_ticks(np.arange(0, 661, 10))
    ax[4].set_xticklabels(parameters[0::10])
    plt.setp(ax[4].get_xticklabels(), rotation=90, fontsize=8)
    ax[4].set_xlabel('parameters')
    handles, labels = ax[0].get_legend_handles_labels()
    fig.legend(handles, labels,loc='upper right')
    plt.show()
    # ss = str(data.corr())
    # s1 = '\t'.join(filter(lambda x: x, ss.split(' ')))
    # output.write(s1)
    # output.close()