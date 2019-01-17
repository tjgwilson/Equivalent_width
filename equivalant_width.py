import numpy as np
import scipy as scp
import matplotlib.pyplot as plt



"""
class that loads data and works out the equivilant width of the line profile
"""
class equivWidth:
    from scipy.interpolate import interp1d
    def __init__(self,input_file):
        self.input_file = input_file

    """
    Loads data from specified data input_file
    updates to implement - different file format handling
    """
    def load_data(self):
        data = np.loadtxt(self.input_file)
        print("shape of data is {}".format(np.shape(data)))
        return data

    def plot_profile(self):
        data = equivWidth.load_data(self)
        plt.plot(data[:,0],data[:,1])
        plt.show()

    def interpCurve(self):
        func = interp1d(x, y, kind='cubic')
        return func


a=equivWidth("halpha_035_12818.dat")
a.plot_profile()
