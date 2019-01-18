import numpy as np
import scipy as scp

import matplotlib.pyplot as plt
from scipy.interpolate import interp1d



"""
class that loads data and works out the equivilant width of the line profile
"""
class equivWidth:
    def __init__(self,input_file,x,y):
        self.input_file = input_file
        self.x = x
        self.y = y

    """
    Loads data from specified data input_file
    updates to implement - different file format handling
    """
    def load_data(self):
        data = np.loadtxt(self.input_file)
        shape = np.shape(data)
        print("Input array size: {}".format(np.shape(data)))
        return data

    def plot_profile(self,data,func):
        plt.plot(data[:,self.x],data[:,self.y],'bx',markersize=5,label="Data")
        min = np.amin(data[:,self.x])
        max = np.amax(data[:,self.x])
        x = np.linspace(min,max, num=10*np.shape(data)[0], endpoint=True)
        plt.plot(x,func(x),'k',linewidth=1,label="Cubic fit")
        plt.legend()
        plt.show()

    def interpCurve(self,data):
        func = interp1d(data[:,self.x], data[:,self.y], kind='cubic')
        return func

    def find_line_limits(self,data,func):
        dif = 0.001
        rows = np.shape(data)[0]
        av_len = int(0.05 * rows)
        upper_bound = rows - 1
        lower_bound = 0
        run_av = 0.0
        temp_av = 0.0

        for i in range(av_len):
            run_av += data[i,self.y]
        run_av /= av_len

        for i in range(av_len,int(rows/2.0),1):
                temp_av = run_av
                run_av += data[i,self.y]/av_len
                run_av -= data[i-av_len,self.y]/av_len
                if(abs(run_av-temp_av) >= dif):
                    lower_bound = i
                    break

        run_av = 0.0
        for i in range(rows-1,av_len,-1):
            run_av += data[i,self.y]
        run_av /= av_len

        for i in range((rows-1-av_len),int(rows/2.0),-1):
                temp_av = run_av
                run_av += data[i,self.y]/av_len
                run_av -= data[i-av_len,self.y]/av_len
                if(abs(run_av-temp_av) >= dif):
                    upper_bound = i
                    break
        return(lower_bound,upper_bound)

    def find_baseline(self,data,func):
        rows = np.shape(data)[0]
        bounds = equivWidth.find_line_limits(self,data,func)
        length =  bounds[0] + (rows - bounds[1])

        continuum = np.zeros((length,2))
        count = 0
        for i in range(0,bounds[0],1):
            continuum[count,0] = data[i,self.x]
            continuum[count,1] = data[i,self.y]
            count += 1
        for i in range(bounds[1],rows,1):
            continuum[count,0] = data[i,self.x]
            continuum[count,1] = data[i,self.y]
            count += 1

        order = 2
        coef = np.polyfit(continuum[:,0],continuum[:,1],order)

        data_cont = np.zeros((rows,2))
        data_cont[:,0] = data[:,self.x]

        x=np.linspace(-500,500,100,)

        # for i in range(1):
        #     cont = 0.0
        #     for j in range(order,-1,-1):
        #         cont += coef[j]*np.power(data_cont[i,0],j)
        #         print(data_cont[i,0],cont)
        #     data_cont[i,1] = data[i,1] - cont

        #plt.plot(continuum[:,0],continuum[:,1],'r*') #debugging
        plt.plot(data_cont[:,0],data_cont[:,1])
        plt.show()


    def calc_width(self):
        data = equivWidth.load_data(self)
        func = equivWidth.interpCurve(self,data)
        equivWidth.find_baseline(self,data,func)
        #equivWidth.plot_profile(self,data,func)


a = equivWidth("halpha_035_12818.dat",0,1).calc_width()
