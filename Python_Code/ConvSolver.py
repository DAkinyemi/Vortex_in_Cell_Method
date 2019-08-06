import numpy as np
import copy
from scipy import signal

class ConvSolver:
    def __init__(self):
        self.data = []

    def kernel(self, x):
        if x == 0:
            return 0
        else:
            return -1/(2*np.pi)*np.log(x)

    # Assumes grid is a 1D grid and the domain is square
    def __init__(self, grid, fullGrid = True):
        self._m_dx = np.abs(grid[0] - grid[1])
        self.fullGridBool = fullGrid
        if(not self.fullGridBool):
            shape = grid.shape
            numPoints = shape[0]
            self._m_M = int(np.log2(numPoints))
        else:
            shape = grid.shape
            numPoints = shape[0]
            self._m_M = int(np.log2(numPoints - 1))
        self.grid = grid
        xCen = grid[int(len(grid)/2)]
        yCen = xCen

        #This pads a sub-maximal domain with zeros at the boundary to pass to the convolve function
        gridSize = 2**self._m_M
        self.KernelArray = np.zeros((gridSize+1, gridSize+1))
        for indexX, X in np.ndenumerate(grid):
            for indexY, Y in np.ndenumerate(grid):
                self.KernelArray[indexX, indexY] = self.kernel((X - 0.5)**2 + (Y - 0.5)**2)

    def solve(self, RHSArray):

        if (not self.fullGridBool):
            #This pads a sub-maximal domain with zeros at the boundary to pass to the convolve function
            gridSize = 2**self._m_M
            inputArray = np.zeros((gridSize+1, gridSize+1))
            inputArray[0:-1, 0:-1] = copy.deepcopy(RHSArray)
        else:
            #Not necessary to do a deep copy here but I want to be transparent about data being overwritten
            inputArray = copy.deepcopy(RHSArray);

        soln = signal.convolve2d( self.KernelArray, inputArray, mode='same', boundary = 'fill', fillvalue = '0.0')
        scale = 1/(len(self.grid)**2)
        soln *= scale

        if (not self.fullGridBool):
            outputArray = soln[0:-1, 0:-1]
        else:
            outputArray = soln;

        return outputArray
