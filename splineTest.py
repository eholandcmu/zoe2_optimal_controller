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


def calc_spline(points, del_d=0.5, b=0, tau=0, f=1):
    virt_start = endpoint_extrapolation(points[:2,:],-f)
    virt_end = endpoint_extrapolation(points[-2:,:],f)

    points = np.concatenate((virt_start,points,virt_end))

    fullSpline = []
    fullSplineHeadings = []
    for i in range(points.shape[0]-3):
        splineSeg, headingSeg = KB_spline(points[i:i+4,:], del_d, b, tau)
        fullSpline.append(splineSeg[:-1,:]) #don't include the last points because they will be the first points in the next segment
        fullSplineHeadings.append(headingSeg[:-1])
        # fullSplineHeadings = [fullSplineHeadings, headingSeg]

    fullSpline.append(splineSeg[-1,:]) #add the last points to the very end
    fullSplineHeadings.append([headingSeg[-1]])

    return np.vstack(fullSpline), np.concatenate(fullSplineHeadings)
        


def endpoint_extrapolation(pts,factor=1):
    """ Computes a virtual point for the spline endpoints based on linear extrapolation

    Args:
        pts: numpy array of 2 points, the second of which is the endpoint
        factor: relative length along line of the new point away from the other 2. 1 means the old endpont will be at the midpoint

    Returns:
        endpt: the new virtual endpoint [x,y]
    """

    x = pts[0,0] + (factor+1)*(pts[1,0]-pts[0,0])
    y = pts[0,1] + (factor+1)*(pts[1,1]-pts[0,1])

    return np.column_stack((x,y))


def KB_spline(pts, del_d=0.5, b=0, tau=0):
    """ Computes the modified Kochanek–Bartels spline segment between the points. b = tau = 0 for Catmull-Rom spline

    Args:
        pts: numpy array of 4 points
        N: num of interpolation points of spline
        del_d: the desired distance between spline points
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

    # del_d = 0.2 #m
    N_int = 1000
    t_int = np.linspace(0,1,N_int)
    T_int = np.array([np.ones(N_int), 2*t_int, 3*(t_int**2)])
    X_p = T_int.T @ C_x[1:]
    Y_p = T_int.T @ C_y[1:]
    L = np.trapz(np.sqrt(X_p**2+Y_p**2),t_int)

    N = np.round(L/del_d).astype(int)

    print(f'Segment length: {L}')
    print(f'Num pts in Segment: {N}')

    t = np.linspace(0,1,N)
    T = np.array([np.ones(N), t, t**2, t**3])
    X = T.T @ C_x
    Y = T.T @ C_y

    CR_spline = np.column_stack((X,Y))

    T_p = np.array([np.ones(N), 2*t, 3*(t**2)])

    X_p = T_p.T @ C_x[1:]
    Y_p = T_p.T @ C_y[1:]

    headings = np.arctan2(Y_p,X_p)  #*180/np.pi

    return CR_spline, headings

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

def generate_trajectory(points, params):
    
    v_b = params["v_b"] #nominal body speed of the rover (m/s)
    dt = params["dt"] #time step (s)
    r = params["r"] #rover wheel radius (m)
    

    del_d = v_b*dt  #the desired distance between spline points

    fullSpline, fullHeadings = calc_spline(points, del_d, b=0, tau=0, f=0.5)

    N = len(fullHeadings)

    Xref = [[fullSpline[i,0],fullSpline[i,1],fullHeadings[i],0,0] for i in range(N)]

    Uref = [[v_b/r * np.ones(4)] for _ in range(N)]

    return Xref, Uref

def main():
    
    # points = np.array([[0,-1], #starting point
    #                 [0.5, 0],
    #                 [0.5, 9],
    #                 [1.5, 9],
    #                 [1.5, 0],
    #                 [2.5,0]])
    
    points = np.array([[0,0], #starting point
                [4, 5],
                # [0, 10],
                [-4, 15],
                [0, 20]])
    

    fig, ax = plt.subplots()
    # ax.set_xlim((-2,7))
    # ax.set_ylim((-2,10))
    ax.axis('equal')

    # plot_rows(points,ax)
    del_d = 0.25  #the desired distance between spline points.  Should be set to v_b*dt

    fullSpline, fullHeadings = calc_spline(points, del_d, b=0, tau=0, f=0.5)

    # fig, ax = plt.subplots()
    
    ax.scatter(fullSpline[:,0], fullSpline[:,1],c='green',s=2) #interpolated spline points
    ax.scatter(points[:,0], points[:,1],c='red',s=6) #waypoints

    fig,ax2 = plt.subplots()

    ax2.plot(fullHeadings, marker='.')

    plt.show()

    #### sample usage: ####
    # points = np.array([[0,0], #starting point
    #         [4, 5],
    #         [-4, 15],
    #         [0, 20]])
    # params = {"v_b":10,
    #           "dt":4,
    #           "r":.375}
    
    # Xref,Uref = generate_trajectory(points, params) 

if __name__ == "__main__":
    main()