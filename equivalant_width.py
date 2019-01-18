import numpy as np
import scipy as scp

import matplotlib.pyplot as plt
from scipy.interpolate import interp1d



"""
class that loads data and works out the equivilant width of the line profile
"""
class equivWidth:

    def __init__(self,input_file,x,y,order,wavelength):
        self.input_file = input_file
        self.x = x
        self.y = y
        self.order = order
        self.wavelength = wavelength

    """
    Loads data from specified data input_file
    updates to implement - different file format handling
    """
    def load_data(self):
        data = np.loadtxt(self.input_file)
        shape = np.shape(data)
        print("Input array size: {}".format(np.shape(data)))
        return data


    def vel_2_wavelength(self,vel): #lamda is the central wavelength
        c = 299792458.0
        ratio = (c + vel) / (c - vel)
        return self.wavelength + (self.wavelength / ratio)

    def plot_profile(self,data,func,poly,width,continuum):
        plt.plot(data[:,self.x],data[:,self.y],'bx',markersize=5,label="Data")
        min = np.amin(data[:,self.x])
        max = np.amax(data[:,self.x])
        x = np.linspace(min,max, num=10*np.shape(data)[0], endpoint=True)
        plt.plot(x,func(x),'k',linewidth=1,label="Cubic fit")
        plt.plot(x,poly(x),'r',linewidth=1,label="continuum fit")
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
        for i in range(0,bounds[0],1):          #write continuum data to array for polyfit
            continuum[count,0] = data[i,self.x]
            continuum[count,1] = data[i,self.y]
            count += 1
        for i in range(bounds[1],rows,1):
            continuum[count,0] = data[i,self.x]
            continuum[count,1] = data[i,self.y]
            count += 1

        coef = np.polyfit(continuum[:,0],continuum[:,1],self.order)
        poly = np.poly1d(coef)
        data_cont = np.zeros((rows,2))
        data_cont[:,0] = data[:,self.x]

        return poly

    def trapezium(self,func,lower,upper):
        n = int(abs(upper-lower)) * 10
        diff = abs(upper-lower)/n

        y_vals= np.zeros((n))
        x_vals = np.linspace(lower,upper,n,endpoint=False)
        y_vals = func(x_vals)
        area = 0.0
        for i in range(1,n-1,1):
            area += (2.0 * y_vals[i])
        area += (y_vals[0] + y_vals[-1])
        area *= (0.5 * diff)

        return area

    def calc_width(self):
        data = equivWidth.load_data(self)
        func = equivWidth.interpCurve(self,data)
        poly = equivWidth.find_baseline(self,data,func)

        x_min = np.amin(data[:,self.x])
        x_max = np.amax(data[:,self.x])
        prof_area = equivWidth.trapezium(self,func,x_min,x_max)

        prof_area -= equivWidth.trapezium(self,poly,x_min,x_max)
        continuum = np.mean(poly(data[:,0]))
        width = prof_area / continuum

        print("width ",width)
        print("prof_area ", prof_area)
        equivWidth.plot_profile(self,data,func,poly,width,continuum)
        return width

    def test_data(self):
        count = -50
        y=np.zeros((100,2))
        for i in range(100):
            if(count<=-5 or count >5):
                y[i,0] = count
                y[i,1] = 2
                count += 1
            else:
                y[i,0] = count
                y[i,1] = 12.
                count += 1
            np.savetxt('test.out',y)





x_col = 0
y_col = 1
order = 6
wavelength = 128.18
a = equivWidth("halpha_035_12818.dat",x_col,y_col,order,wavelength).calc_width()
