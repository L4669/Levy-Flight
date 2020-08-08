#! /usr/bin/python3.6

import math
import numpy as np
import random
import matplotlib.pyplot as plt
from matplotlib import animation
from IPython.display import HTML


def proximity_sensor(Tx, Ty, xi, yi, l_min):
    l_prox = 0
    tx_ctf = 0
    ty_ctf = 0
    ctf = 0
    for tx, ty in zip(Tx, Ty):
        l_prox = ((tx-xi)**2 + (ty-yi)**2)**0.5
        if l_prox <= l_min:
            ctf = 1
            tx_ctf = tx
            ty_ctf = ty
            break
    return ctf, tx_ctf, ty_ctf


def proximity_sensor_path(Tx, Ty, xi, yi, xf, yf, l_min):
    tx_ctf = 0
    ty_ctf = 0
    x_path = 0
    y_path = 0
    ctf = 0
    a = -(yf-yi)/(xf-xi)
    b = 1
    c = xi*(yf-yi)/(xf-xi) - yi
    xmax = max(xi, xf)
    xmin = min(xi, xf)
    ymax = max(yi, yf)
    ymin = min(yi, yf)
    for tx, ty in zip(Tx, Ty):
        l_prox = abs(a*tx + b*ty + c)/(a**2 + b**2)**0.5
        x_path = (b*(b*tx - a*ty) - a*c)/(a**2 + b**2)
        y_path = (a*(-b*tx + a*ty) - b*c)/(a**2 + b**2)
        if (l_prox <= l_min) and (x_path <= xmax) and \
                (x_path >= xmin) and (y_path <= ymax) and (y_path >= ymin):
            ctf = 1
            tx_ctf = tx
            ty_ctf = ty
            break
    return ctf, x_path, y_path, tx_ctf, ty_ctf


def levy(mu, l_min, A_ll, A_ul):
    r = random.random()
    theta = random.uniform(A_ll, A_ul)
    l = l_min*(1-r)**(1/(1-mu))
    return l, theta

def brownian(l_ll, l_ul, A_ll, A_ul):
    theta = random.uniform(A_ll, A_ul)
    l  = random.uniform(l_ll, l_ul)
    return l, theta


def flight_simulation(pdf, N, R, C, l_min, mu, Tx, Ty, xi, yi, predator_range):
    # capture the flag
    ctf = 0 

    # initialization of temp variables
    l = []
    X = [xi]
    Y = [yi]
    n = 1
    l_prox = 0
    
    # simulation loop
    while (1):
        # check if the prey lies within the proximity when the predator is
        # yet to start
        ctf, Tx_ctf, Ty_ctf = proximity_sensor(Tx, Ty, xi, yi, l_min)
        
        # condition for stopping the simulation
        if (ctf == 1):
            X.append(Tx_ctf)
            Y.append(Ty_ctf)
            n += 1
            print("Success with n =", n)
            print("Captured Prey at:", Tx_ctf, Ty_ctf)
            break
        
        if (n == N):
            print("Simulation stopped as max. n reached")
            break
        
        # generating parameters, length & direction, as per Levy Distribution
        if pdf == "levy":
            l_var, theta = levy(mu, l_min, 0, 2*math.pi)
        elif pdf == "brownian":
            l_var, theta = brownian(l_min, 50, 0, 2*math.pi)
        else:
            print("PDF not supported")
            return -1

        xf = xi + l_var*math.cos(theta)
        yf = yi + l_var*math.sin(theta)
        
        # distance of the next point from origin
        l_origin = ((xf-int(R/2))**2 + (yf-int(C/2))**2)**0.5
      
        # if the next point on the path lies outside the 
        # range of predator, it need to be handled separately
        if (l_origin <= predator_range):
            ctf, x_path, y_path, Tx_ctf, Ty_ctf = \
                    proximity_sensor_path(Tx, Ty, xi, yi, xf, yf, l_min)
            if (ctf == 1):
                X.append(x_path)
                Y.append(y_path)
                X.append(Tx_ctf)
                Y.append(Ty_ctf)
                n += 2
                print("Success with n =", n)
                print("Captured Prey at:", Tx_ctf, Ty_ctf)
                break
            X.append(xf)
            Y.append(yf)
            xi = xf
            yi = yf
            n += 1
        else:
            xf = predator_range*math.cos(theta) + int(R/2)
            yf = predator_range*math.sin(theta) + int(C/2)
            ctf, x_path, y_path, Tx_ctf, Ty_ctf = \
                    proximity_sensor_path(Tx, Ty, xi, yi, xf, yf, l_min)
            if (ctf == 1):
                X.append(x_path)
                Y.append(y_path)
                X.append(Tx_ctf)
                Y.append(Ty_ctf)
                n += 2
                print("Success with n =", n)
                print("Captured Prey at:", Tx_ctf, Ty_ctf)
                break
                
            X.append(xf)
            Y.append(yf)
            xi = xf
            yi = yf
            
            # theta to be selected so that the predator remains inside the circle
            if (theta >= 0) and (theta <= 0.5*math.pi):
                A_ll = math.pi - (0.5*math.pi-theta)
                A_ul = 1.5*math.pi + theta
            elif (theta >= 0.5*math.pi) and (theta <= math.pi):
                A_ll = 1.5*math.pi - (1*math.pi-theta)
                A_ul = 2*math.pi + theta - 0.5*math.pi
            elif (theta >= math.pi) and (theta <= 1.5*math.pi):
                A_ll = 0*math.pi - (1.5*math.pi-theta)
                A_ul = 0.5*math.pi + theta - math.pi
            else:
                A_ll = 0.5*math.pi - (2*math.pi-theta)
                A_ul = math.pi + theta - 1.5*math.pi
            
            if pdf == "levy":  
                l_var, theta = levy(mu, l_min, A_ll, A_ul)
            elif pdf == "brownian":
                l_var, theta = brownian(0, 10, A_ll, A_ul)
            else:
                print("PDF not supported")
                return -1
            xf = xi + l_var*math.cos(theta)
            yf = yi + l_var*math.sin(theta)
            
            ctf, x_path, y_path, Tx_ctf, Ty_ctf = \
                    proximity_sensor_path(Tx, Ty, xi, yi, xf, yf, l_min)
            if (ctf == 1):
                X.append(x_path)
                Y.append(y_path)
                X.append(Tx_ctf)
                Y.append(Ty_ctf)
                n += 2
                print("Success with n =", n)
                print("Captured Prey at:", Tx_ctf, Ty_ctf)
                break
                
            X.append(xf)
            Y.append(yf)
            xi = xf
            yi = yf
            n += 2
    
    return X, Y, Tx_ctf, Ty_ctf
  
if __name__ == "__main__":
    # Intial Conditions
    l_min = 10
    mu = 2.5
    N = 10000

    # grid size
    R = 1000
    C = 1000

    # coordinate of prey on the grid
    num_prey = 5
    Tx = []
    Ty = []

    for i in range(num_prey):
        Tx.append(random.randint(0, R))
        Ty.append(random.randint(0, C))
    
    # initial location of predator
    xi = int(R/2)
    yi = int(C/2)
    
    # circle of some radius from initial point
    # here it is assumed as diagonal of the square grid defined above
    predator_range = (R**2 + C**2)**0.5/2 
    
    # range for plotting purpose
    xmax = xi + predator_range
    ymax = yi + predator_range
    xmin = xi - predator_range
    ymin = yi - predator_range

    pdf = "levy"  
    X, Y, Tx_ctf, Ty_ctf = flight_simulation(pdf, N, R, C, l_min, mu, Tx, Ty,\
            xi, yi, predator_range)
    
    pdf = "brownian"
    X1, Y1, Tx_ctf1, Ty_ctf1 = flight_simulation(pdf, N, R, C, l_min, mu, Tx, Ty,\
            xi, yi, predator_range)
    
    # plot routing
    fig, ax = plt.subplots()
    fig.set_dpi(100)
    fig.set_size_inches(7, 7)
    plt.grid(linewidth=0.5)

    ax = plt.axes(xlim=(xmin, xmax), ylim=(ymin, ymax))
    plt.grid(linewidth=0.5)
    plt.gca().set_aspect("equal")

    plt.plot(X, Y, 'b')
    plt.plot(X1, Y1, 'y')

    plt.plot(int(R/2), int(C/2), '*')
    ax.add_patch(plt.Circle((int(R/2), int(R/2)), \
            predator_range, edgecolor='orange', facecolor='white'))
    for tx, ty in zip(Tx, Ty):
        ax.add_patch(plt.Circle((tx, ty), l_min, color='black'))
        
    captured_target = plt.Circle((Tx_ctf, Ty_ctf), l_min, color='r')
    ax.add_patch(captured_target)
    captured_target_1 = plt.Circle((Tx_ctf1, Ty_ctf1), l_min, color='r')
    ax.add_patch(captured_target_1)
    plt.show()

