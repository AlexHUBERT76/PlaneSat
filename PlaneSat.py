import numpy as np
from numpy import linalg as LA
import random as rnd
from math import *


def rand(fCA=0):
    # module returning random values as a tuple (x,x,z)
    # one can tune the values here to make an artificial plan
    x = rnd.choice((-1, 1)) * rnd.random() * 100
    y = rnd.choice((-1, 1)) * rnd.random() * 100
    z = rnd.choice((-1, 1)) * rnd.random() * 100 * fCA  # fCA = c/a

    return np.array([x, y, z])


def Centroid(r):
    # derivation of centroid r0
    N = len(r)
    r0 = np.divide(np.sum(r, axis=0), N)
    return r0


def PlaneFitting(r, D=3):
    # module taking r (np.array) and dimensionality D of problem then returning the best fitting orthogonal eigen vector of an hypothetic plane as a vector tuple (x,y,z)

    N = len(r)
    r0 = Centroid(r)
    T0 = np.zeros((3, 3))

    i = 0
    while i < N:
        # loop calculating the moment of inertia tensor T0
        rs = np.subtract(r[i], r0)
        sq = np.dot(rs, rs)
        id = np.multiply(sq, np.identity(3))
        id = np.subtract(id, np.outer(rs, np.transpose(rs)))
        T0 = np.add(T0, id)

        i += 1

    # returning the eigen values (w) and eigen vectors (v) of tensor T0
    w, v = LA.eig(T0)
    v = np.transpose(v)

    # returning index of highest eigen value
    i = np.argmax(w)

    if D == 2:
        w[i] = 0
        # returning index of second highest eigen value
        i = np.argmax(w)
        return v[i]

    if D == 1:
        w[i] = 0
        # returning index of second third eigen value
        i = np.argmax(w)
        return v[i]

    return v[i]




def RndValArr(N=1000, fCA=0):
    # module returning random values as N tuples (x,x,z)
    tab = np.empty([N, 3])
    for i in range(0, N):
        val = rand(fCA)
        tab[i] = val
        i += 1
    return tab

def RndRotMat():
    # module returning random 3D rotation matrix
    a = rnd.choice((-1, 1)) * rnd.random() * 2 * np.pi
    b = rnd.choice((-1, 1)) * rnd.random() * 2 * np.pi
    c = rnd.choice((-1, 1)) * rnd.random() * 2 * np.pi

    Rotx = ([[1,0,0],[0,np.cos(a),-np.sin(a)],[0,np.sin(a),np.cos(a)]])
    Roty = ([[np.cos(b),0,np.sin(b)],[0,1,0],[-np.sin(b),0,np.cos(b)]])
    Rotz = ([[np.cos(c),-np.sin(c),0],[np.sin(c),np.cos(c),0],[0,0,1]])

    matrix = np.matmul(Rotx,Roty)
    matrix = np.matmul(matrix,Rotz)

    return matrix


def GalCoord(n):
    # module returning a cartesian coordinate vector n(x,y,z) into a galactic coordinate vector n(l,b) (in deg°)
    x = n[0]
    y = n[1]
    z = n[2]

    if x != 0:
        l = np.arctan(y / x)
        l = np.rad2deg(l)
    else:
        l = 0.0
    if x != 0 and y != 0:
        b = np.arctan(z / sqrt(x * x + y * y))
        b = np.rad2deg(b)
    else:
        b = 90.0
    return np.array([l, b])


def GalPlot(r, n = ([0,0,0]), size = 100, title = ""):  # plot galactic plan given galaxies positions r(x,y,z) and normal vector n(x,y,z)

    from mpl_toolkits import mplot3d
    import matplotlib.pyplot as plt

    r0 = Centroid(r)

    fig = plt.figure()
    ax = plt.axes(projection="3d")

    z_points = []
    x_points = []
    y_points = []
    i = 0
    while i < len(r):
        z_points.append(r[i, (2)])
        x_points.append(r[i, (0)])
        y_points.append(r[i, (1)])
        i = i + 1

    ax.scatter3D(x_points, y_points, z_points, c=z_points, cmap='hsv')  # plotting galaxies position

    n = np.multiply(n, size)

    plt.plot([r0[0], r0[0] + n[0]], [r0[1], r0[1] + n[1]], [r0[2], r0[2] + n[2]])  # plotting normal

    plt.title(title)

    plt.show()

    return ()


def main():

    Ngal = 10
    Nmin = 0
    Nmax = 30

    stat = 10000
    CA = 0.5
    result = np.zeros([Nmax - Nmin + 1, 3])
    tracker = 0

    for Galpol in range(Nmin, Nmax + 1):
        i = Galpol + Ngal
        Nn = []
        Nnproj = []
        rproj = np.zeros([i, 3])


        for k in range(0, stat):
            rndproj = rnd.choice((0, 1, 2))
            matrot = RndRotMat()

            normal = [0, 0, 1]
            normalext = np.zeros([3, 3])
            normalext[0] = normal
            normalext = np.matmul(normalext, matrot)
            normal = normalext[0]


            #r = RndValArr(i, CA)
            r = np.concatenate((RndValArr(Ngal,CA),RndValArr(Galpol,1)))

            for j in range(0, i-1):
                riext = np.zeros([3, 3])
                riext[0] = r[j]
                riext = np.matmul(riext,matrot)
                r[j] = riext[0]



            projvector = np.zeros([1, 3])
            if rndproj == 0:
                projvector = ([[0,1,1]])
            elif rndproj == 1:
                projvector = ([[1, 0, 1]])
            elif rndproj == 2:
                projvector = ([[1,1,0]])

            rproj = np.multiply(r, projvector)

            n2D = PlaneFitting(rproj, 2)
            normalproj = np.multiply(normal, projvector)



            dot = np.dot(n2D,normal)
            dotproj = np.dot(normalproj, n2D)


            delta = np.divide(dot,LA.norm(normal)*LA.norm(n2D))
            delta = np.arccos(delta)
            delta = np.rad2deg(delta)

            deltaproj = np.divide(dotproj,LA.norm(normalproj)*LA.norm(n2D))
            deltaproj = np.arccos(deltaproj)
            deltaproj = np.rad2deg(deltaproj)



            if delta > 90:
                delta = 180-delta

            if deltaproj > 90:
                deltaproj = 180-deltaproj

            Nnproj.append(deltaproj)
            Nn.append(delta)
            tracker = tracker + 1
            if tracker % 1000 == 0:
                print("Progress :", tracker * 100.0 / (stat * (Nmax - Nmin)), "%")

        mean_angle = np.mean(Nn)
        mean_angleproj = np.mean(Nnproj)
        N = Galpol
        result[Galpol] = [N, mean_angle, mean_angleproj]



    import matplotlib.pyplot as plt

    xresult = []
    yresult = []
    zresult = []

    for i in range(0, len(result)):
        xresult.append(result[i, (0)])
        yresult.append(result[i, (1)])
        zresult.append(result[i, (2)])

    plt.plot(xresult, yresult, "r", label = "Error with real normal")
    plt.plot(xresult, zresult, "b", label = "Error with projected normal")

    plt.legend()

    plt.xlabel("Satellites number polution N")
    plt.ylabel("Absolute error in degree")

    title = "Mean error on PlaneFitting2D over " + str(stat) + " trials per case and with C/A = " + str(CA)
    plt.title(title)
    plt.show()


if __name__ == '__main__':
    main()



