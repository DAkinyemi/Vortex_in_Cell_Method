#Importing required libraries
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
import random
import copy
np.set_printoptions(suppress=True)
# import sys
# sys.path.insert(0, "./../../vortexincell/src/wrappedCode/")
import ConvSolver
import importlib
importlib.reload(ConvSolver)

#Vortex Particle-in-Cell Method
class VPIC:
    #Takes points and breaks into x, y, vorticity strength, dx and dy
    def __init__(self,points,dx,dy):
        self.points = points
        self.dx = dx
        self.dy = dy
        self.x = points[:,0]
        self.y = points[:,1]
        self.size = points[:,2]

    #Interpolates the points to 4 corners indexed
    def interpolate(self):
        interpolated_points = []
        for i in range(0,len(self.points)):
            stuff = np.zeros((4,2))
            stuff[0][0] = math.floor(self.x[i]/self.dx) #bottom left x
            stuff[0][1] = math.floor(self.y[i]/self.dy) #bottom left y
            stuff[1][0] = math.floor(self.x[i]/self.dx)  + 1 #bottom right x
            stuff[1][1] = math.floor(self.y[i]/self.dy) #bottom right y
            stuff[2][0] = math.floor(self.x[i]/self.dx) + 1 #top right x
            stuff[2][1] = math.floor(self.y[i]/self.dy) + 1 #top right y
            stuff[3][0] = math.floor(self.x[i]/self.dx) #top left x
            stuff[3][1] = math.floor(self.y[i]/self.dy) + 1 #top left y
            interpolated_points.append(stuff)
        interpolated_points = np.concatenate(interpolated_points)
        return interpolated_points

    #Deposits the weights relative to the distance from points
    #Where xi and yi are the interpolated x and y coordinates
    def weight_deposition(self,xi,yi):
        array_of_newsizes = []
        for i in range(0,len(self.points)):
            for k in range((4*i),((4*i)+4)):
                volfact = self.size[i]#/(self.dx*self.dy)
                if self.y[i] - (self.dy * yi[k]) > 0:
                    yweights = 1 - ((self.y[i] - (self.dy * yi[k]))/self.dy)
                else:
                    yweights = 1 - abs((self.y[i] - (self.dy * yi[k]))/self.dy)
                if self.x[i] - (self.dx * xi[k]) > 0:
                    xweights = 1 - ((self.x[i] - (self.dx * xi[k]))/self.dx)
                else:
                    xweights = 1 - abs((self.x[i] - (self.dx * xi[k]))/self.dx)
                totweight = volfact * (yweights) * (xweights)
                array_of_newsizes.append([(self.dx * xi[k]),(self.dy * yi[k]), totweight])
        array_of_newsizes = np.asarray(array_of_newsizes)
        return (array_of_newsizes)

    #Finds the places where two arrays overlap in both x and y
    def findliketerms(self,list1,list2):
        liketerms = []
        for i in range(0,len(list1)):
            for k in range(i+1,len(list2)):
                if list1[i][0] == list2[k][0] and list1[i][1] == list2[k][1]:
                    liketerms.append([i,k])
        liketerms = np.asarray(liketerms)
        liketerms = np.resize(liketerms,(len(liketerms),2))
        return liketerms

    def findliketerms2(self,list1,list2):
        liketerms = []
        for i in range(0,len(list1)):
            for k in range(0,len(list2)):
                if list1[i][0] == list2[k][0] and list1[i][1] == list2[k][1]:
                    liketerms.append([i,k])
        liketerms = np.asarray(liketerms)
        liketerms = np.resize(liketerms,(len(liketerms),2))
        return liketerms

    #Combines the vorticcity strengths of the same x and y coordinate
    def combineliketerms(self, array_of_newsizes, liketerms, interpolated_points):
        for k in range(0,len(liketerms)):
            array_of_newsizes[liketerms[k][0]] =+ array_of_newsizes[liketerms[k][0]] + array_of_newsizes[liketerms[k][1]]
        array_of_newsizes = np.delete(array_of_newsizes, ((liketerms[:,1]).tolist()), axis = 0)
        interpolated_points = np.delete(interpolated_points, ((liketerms[:,1]).tolist()), axis = 0)
        data = np.zeros((len(interpolated_points),3))
        data[:,0] = interpolated_points[:,0] * self.dx
        data[:,1] = interpolated_points[:,1] * self.dy
        data[:,2] = array_of_newsizes
        if np.any(data) == 0:
            return array_of_newsizes
        else:
            return data

    #creates an x-grid & runs imported Poisson Solver
    #RHS = Right Hand Side
    def RHSstart(self, xlengthstart, xlengthstop, combinedliketerms):
        x_grid = np.linspace(xlengthstart/self.dx, xlengthstop/self.dx, (N+1))
        x_grid = x_grid[:-1]
        #print(x_grid*self.dx)
        gridSize = int(len(x_grid))
        RHSArray = np.zeros((int(gridSize), int(gridSize)))

        for item in combinedliketerms:
            #print(item)
            try:
                RHSArray[int(item[0]/self.dx), int(item[1]/self.dx)] = item[2]
            except IndexError:
                print("Error: Particle was deleted because it was indexed outside of accepted range; results are inconclusive.")
                pass
        #print(RHSArray)
        PS = ConvSolver.ConvSolver(x_grid*self.dx, False)
        solution = PS.solve(RHSArray)
        return solution

    #Takes the solution from the poisson solver and run potential difference
    def potentialcalculator(self, points2):
        potentials = []
        for i in range(1,len(points2)-1):
            for k in range(1,len(points2)-1):
                ycentraldiff = ((points2[i, k+1]) - (points2[i, k-1])) / (2 * self.dx)
                xcentraldiff = ((points2[i+1, k]) - (points2[i-1, k])) / (2 * self.dy)
                potentials.append(([self.dx * i, self.dy * k , ycentraldiff, -xcentraldiff]))
        potentials = np.asarray(potentials)
        potentials = np.resize(potentials,(len(potentials),4))
        #plt.quiver(potentials[:,0], potentials[:,1], potentials[:,2], potentials[:,3])
        return potentials

    #Field Interpolation calculation
    def weight_deposition2(self,x,y,xi,yi,xweight,yweight):
        if yi - y > 0:
            yweights = 1.0 - ((yi - y)/self.dy)
        else:
            yweights = 1.0 - abs((yi - y)/self.dy)
        if xi - x > 0:
            xweights = 1.0 - ((xi - x)/self.dx)
        else:
            xweights = 1.0 - abs((xi - x)/self.dx)
        xvelocity = xweight * (yweights) * (xweights)
        yvelocity = yweight * (yweights) * (xweights)
        return (xvelocity, yvelocity)

    #Looks for the interpolated velocity for each correspoiding point to interpolate back
    def interpolation_backtopoints(self,newdata):
        velocities_on_grid = []
        for i in range(0,len(self.points)):
            for k in range(0, len(newdata)):
                #print(newdata[k][0], math.floor(self.points[i][0]/self.dx))
                if (newdata[k][0] == math.floor(self.points[i][0]/self.dx)*self.dx and newdata[k][1] == math.floor(self.points[i][1]/self.dy)*self.dy):
                    #print(i,k,"t1")
                    velocities_on_grid.append([i,k,self.weight_deposition2(self.x[i],self.y[i],newdata[k][0], newdata[k][1], newdata[k][2], newdata[k][3])])
                if (newdata[k][0] == (math.floor(self.points[i][0]/self.dx) + 1)*self.dx and newdata[k][1] == math.floor(self.points[i][1]/self.dy)*self.dy):
                    #print(i,k,"t2")
                    velocities_on_grid.append([i,k,self.weight_deposition2(self.x[i],self.y[i],newdata[k][0], newdata[k][1], newdata[k][2], newdata[k][3])])
                if (newdata[k][0] == (math.floor(self.points[i][0]/self.dx) +1)*self.dx and newdata[k][1] == (math.floor(self.points[i][1]/self.dy) + 1)*self.dy):
                    #print(i,k,"t3")
                    velocities_on_grid.append([i,k,self.weight_deposition2(self.x[i],self.y[i],newdata[k][0], newdata[k][1], newdata[k][2], newdata[k][3])])
                if (newdata[k][0] == math.floor(self.points[i][0]/self.dx)*self.dx and newdata[k][1] == (math.floor(self.points[i][1]/self.dy) + 1)*self.dy):
                    #print(i,k,"t4")
                    velocities_on_grid.append([i,k,self.weight_deposition2(self.x[i],self.y[i],newdata[k][0], newdata[k][1], newdata[k][2], newdata[k][3])])
        velocities_on_grid = np.array(velocities_on_grid)
        #velocities_on_grid = np.resize(velocities_on_grid,(len(velocities_on_grid),4))
        return velocities_on_grid

    #Combines velocities based on the points they belong to
    def combining_velocities(self,velocities_on_grid):
        empty = np.zeros((len(self.points),3))
        for i in range(0,len(velocities_on_grid)):
            for k in range(0,len(velocities_on_grid)):
                #print(i, velocities_on_grid[k][0])
                if i == velocities_on_grid[k][0]:
                    #print(i, velocities_on_grid[k])
                    empty[i][0] = i
                    empty[i][1] =+ empty[i][1] + velocities_on_grid[k][2][0]
                    empty[i][2] =+ empty[i][2] + velocities_on_grid[k][2][1]
        return(empty[:,1:3])

#INITIAL CONDITIONS
#X and Y grid size
xlengthstop = 1
ylengthstop = 1

xlengthstart = 0
ylengthstart = 0

xlength = xlengthstop - xlengthstart
ylength = ylengthstop - ylengthstart

startvalue = 0
particles = 3

#Grid criteria
N = 2**5
dx = (xlengthstop - xlengthstart)/N
dy = (ylengthstop - ylengthstart)/N

weightmin = -1
weightmax = 1

#CREATES Random points to run the test on for if you want to run this method on random points
points = np.zeros((particles,3))
points[:, 0] = np.array(xlengthstop - xlengthstart)*np.random.random_sample(particles,) + xlengthstart
points[:, 1] = np.array(ylengthstop - ylengthstart)*np.random.random_sample(particles,) + ylengthstart
points[:, 2] = (weightmax-weightmin)*np.random.random_sample(particles)+weightmin

#Function that runs everything
def k2(points,dx,dy):
    initial_particles_with_sizes = VPIC(points,dx,dy)
    interpalations = initial_particles_with_sizes.interpolate()
    depositions = initial_particles_with_sizes.weight_deposition(interpalations[:,0], interpalations[:,1])
    liketerms = initial_particles_with_sizes.findliketerms(interpalations,interpalations)
    updated_data = initial_particles_with_sizes.combineliketerms(depositions[:,2],liketerms,interpalations)
    RHS_calculations = initial_particles_with_sizes.RHSstart(xlengthstart, xlengthstop, updated_data)
    potentials = initial_particles_with_sizes.potentialcalculator(RHS_calculations)
    likterms2 = initial_particles_with_sizes.findliketerms2(updated_data[:,0:2],potentials)
    new_data = potentials[likterms2[:,1].tolist()]
    velocities_on_grid = initial_particles_with_sizes.interpolation_backtopoints(new_data)
    return initial_particles_with_sizes.combining_velocities(velocities_on_grid)

#RK2 Solver
def RK2(points,time,N):
    np.set_printoptions(suppress=False)
    h = 140*0.025/N
    numTimeSteps = int(time/h)
    everything = []
    for i in range(0,numTimeSteps):
        pointsnewx = []
        pointsnewy = []
        tpoints = copy.deepcopy(points)
        #print(tpoints,"test")
        F = h * k2(points,dx,dy)
        #print(F,"F")
        for i in range(0,len(points)):
            k1x = F[i][0]
            k1y = F[i][1]
            tpoints[i][0] = tpoints[i][0] + (.5 * k1x)
            tpoints[i][1] = tpoints[i][1] + (.5 * k1y)
        F2 = h * k2(tpoints, dx, dy)
        #print(F2, "F2")
        for i in range(0,len(tpoints)):
            k2x = F2[i][0]
            k2y = F2[i][1]
            pointsnewx.append(points[i][0] + k2x)
            pointsnewy.append(points[i][1] + k2y)
        pointsnewx = np.asarray(pointsnewx)
        pointsnewy = np.asarray(pointsnewy)
        pointsnewx = np.resize(pointsnewx,(len(pointsnewx),1))
        pointsnewy = np.resize(pointsnewy,(len(pointsnewy),1))
        points[:,0] = pointsnewx[:,0]
        points[:,1] = pointsnewy[:,0]
        everything.append(copy.deepcopy(points))
    everything = np.asarray(everything)
    everything = np.concatenate(everything)
    return everything


#Test 1
figure(num=1, figsize=(6, 6), dpi=80, facecolor='w', edgecolor='k')
points_test = np.array([[0.5, 0.5, (1/(dx*dy))]])
test = RK2(points_test, 10, N)
plt.scatter(test[:,0], test[:,1])
plt.xlabel('X Position')
plt.ylabel('Y Position')
plt.title('Particle Position as a function of Time (Test 1)')

#Test 2
figure(num=2, figsize=(6, 6), dpi=80, facecolor='w', edgecolor='k')
points_test2 = np.array([[0.25, 0.5, 0],[0.5,0.5,1/(dx*dy)]])
test2 = RK2(points_test2, 10, N)
plt.scatter(test2[:,0], test2[:,1])
plt.xlabel('X Position')
plt.ylabel('Y Position')
plt.title('Particle Position as a function of Time (Test 2)')

#Test 3
figure(num=3, figsize=(6, 6), dpi=80, facecolor='w', edgecolor='k')
points_test3 = np.array([[0.5, 0.5,1/(dx*dy)],[0.75,0.5,1/(dx*dy)]])
test3 = RK2(points_test3, 10, N*2)
plt.scatter(test3[:,0], test3[:,1])
plt.xlabel('X Position')
plt.ylabel('Y Position')
plt.title('Particle Position as a function of Time (Test 3)')

# Test 4A Set up
np.set_printoptions(suppress=False)
Np = int(1/(dx/2))
hp = 1/Np
total_data = []
for i in range(0,Np):
    for k in range(0,Np):
        if abs(np.sqrt((i*dx - 0.5)**2 + (k*dy - 0.375)**2))<= .12 or abs(np.sqrt((i*dx - 0.5)**2 + (k*dy - 0.625)**2)) <= .12:
            total_data.append([dx *i,dy * k,((hp**2)/(dx*dy))])
        else:
            total_data.append([dx * i, dy * k,0])
total_data = np.array(total_data)
new_data = []
for item in total_data:
    if item[2] != 0:
        new_data.append([item[0],item[1],item[2]])
new_data = np.array(new_data)
figure(num=4, figsize=(6, 6), dpi=80, facecolor='w', edgecolor='k')
plt.scatter(new_data[:,0], new_data[:,1])
plt.xlabel('X Position')
plt.ylabel('Y Position')
plt.title('Inital Particle Position (Test 4A)')

#Test 4A
figure(num=5, figsize=(6, 6), dpi=80, facecolor='w', edgecolor='k')
test4 = RK2(new_data, 12.5, N)
plt.scatter(test4[:,0], test4[:,1],0.1)
plt.xlabel('X Position')
plt.ylabel('Y Position')
plt.title('Particle Position as a function of Time (Test 4A)')

# Test 4B Set up
np.set_printoptions(suppress=True)
Np = int(1/(dx/2))
hp = 1/Np
rbdry = 1/4
total_data = []
for i in range(0,Np):
    for k in range(0,Np):
        if abs(np.sqrt(((i*dx) - 0.5)**2 + ((k*dy) - 0.5)**2)) <= rbdry:
            total_data.append([dx *i,dy * k,((rbdry - abs((i*dx - 0.5)**2 + (k*dy - 0.5)**2))**7)*100000])
        else:
            total_data.append([dx * i, dy * k,0])
total_data = np.array(total_data)

new_data2 = []
for item in total_data:
    if item[2] != 0:
        new_data2.append([item[0],item[1],item[2]])
new_data2 = np.array(new_data2)
figure(num=6, figsize=(6, 6), dpi=80, facecolor='w', edgecolor='k')
plt.scatter(new_data2[:,0], new_data2[:,1] , new_data2[:,2])
plt.xlabel('X Position')
plt.ylabel('Y Position')
plt.title('Inital Particle Position (Test 4B)')

#Test 4B
figure(num=7, figsize=(6, 6), dpi=80, facecolor='w', edgecolor='k')
test5 = RK2(new_data2, 12.5, N)
plt.scatter(test5[:,0], test5[:,1],0.1)
plt.xlabel('X Position')
plt.ylabel('Y Position')
plt.title('Particle Position as a function of Time (Test 4B)')
plt.show()
