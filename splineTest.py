import numpy as np
import matplotlib.pyplot as plt

def  plot_rows(points, ax):
    
    numrows = 5
    numcols = 7

    delr = 1
    delc = 1.5

    radius = 0.25

    # fig, ax = plt.subplots()

    for r in range(numrows):
        for c in range(numcols):
            circle = plt.Circle((r*delr,c*delc),radius)
            ax.add_patch(circle)

    plt.scatter(points[:,0], points[:,1],c='red')


def calc_spline(points, N=20, b=0, tau=0, f=1):
    virt_start = endpoint_extrapolation(points[:2,:],-f)
    virt_end = endpoint_extrapolation(points[-2:,:],f)

    points = np.concatenate((virt_start,points,virt_end))

    fullSpline = []
    for i in range(points.shape[0]-3):
        splineSeg = KB_spline(points[i:i+4,:], N, b, tau)
        fullSpline.append(splineSeg)

    return np.vstack(fullSpline)
        


def endpoint_extrapolation(pts,factor=1):
    """ Computes a vitrual point for the spline endpoints based on linear extrapolation

    Args:
        pts: numpy array of 2 points, the second of which is the endpoint
        factor: relative length along line of the new point away from the other 2. 1 means the old endpont will be at the midpoint

    Returns:
        endpt: the new virtual endpoint [x,y]
    """

    x = pts[0,0] + (factor+1)*(pts[1,0]-pts[0,0])
    y = pts[0,1] + (factor+1)*(pts[1,1]-pts[0,1])

    return np.column_stack((x,y))


def KB_spline(pts, N=20, b=0, tau=0):
    """ Computes the modified Kochanek–Bartels spline segment between the points. b = tau = 0 for Catmull-Rom spline

    Args:
        pts: numpy array of 4 points
        N: num of interpolation points of spline
        b: Kochanek–Bartels bias term
        tau: Kochanek–Bartels tension term

    Returns:
        CR_spline: 
    """
    M = np.array([[0, 1, 0, 0],
                    [-0.5*(1-tau)*(1+b), (1-tau)*b, 0.5*(1-tau)*(1-b), 0],
                    [(1-tau)*(1+b), 0.5*(1-tau)*(1-3*b)-3, -(1-tau)+3, -0.5*(1-tau)*(1-b)],
                    [-0.5*(1-tau)*(1+b), -0.5*(1-tau)*(1-b)+2, 0.5*(1-tau)*(1+b)-2, 0.5*(1-tau)*(1-b)]])

    C_x = M @ pts[:,0]
    C_y = M @ pts[:,1]

    t = np.linspace(0,1,N)
    T = np.array([np.ones(N), t, t**2, t**3])
    X = T.T @ C_x
    Y = T.T @ C_y

    CR_spline = np.column_stack((X,Y))
    return CR_spline

def two_pt_spline(locs,heads,f, N):
    """ Computes a cubic Hermite spline segment between 2 points given their headings

    Args:
        locs: numpy array of 2 point locations: [x,y]
        heads: numpy array of 2 point headings: [dx, dy]
        f: scale factor multiplied to heading directions
        N: num of interpolation points of spline

    Returns:
        spline
    """
    M = np.array([[1, 0, 0, 0],
                  [0, 0, 1, 0],
                  [-3, 3, 1, 2],
                  [2, -2, -1, -1]])

    C_x = M @ [locs[:,0], heads[:,0]]
    C_y = M @ [locs[:,1], heads[:,1]]
    
    t = np.linspace(0,1,N)
    T = np.array([np.ones(N), t, t**2, t**3])
    X = T.T @ C_x
    Y = T.T @ C_y

    spline = np.column_stack((X,Y))
    return spline

def main():
    
    # points = np.array([[0,-1], #starting point
    #                 [0.5, 0],
    #                 [0.5, 9],
    #                 [1.5, 9],
    #                 [1.5, 0],
    #                 [2.5,0]])
    
    points = np.array([[0,0], #starting point
                [4, 5],
                [0, 10],
                [-4, 15],
                [0, 20]])
    

    fig, ax = plt.subplots()
    # ax.set_xlim((-2,7))
    # ax.set_ylim((-2,10))
    plt.axis('equal')

    # plot_rows(points,ax)

    fullSpline = calc_spline(points, 20, b=0, tau=0, f=0.5)

    # fig, ax = plt.subplots()
    
    plt.scatter(fullSpline[:,0], fullSpline[:,1],c='green',s=2) #intrerpolated spline points
    plt.scatter(points[:,0], points[:,1],c='red',s=6) #waypoints
    plt.show()


if __name__ == "__main__":
    main()