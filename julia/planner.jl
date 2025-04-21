
function calc_spline(points::Vector{Vector{Float64}}, del_d::Float64=0.5, b::Float64=0.0, tau::Float64=0.0, f::Float64=1.0)
    virt_start = endpoint_extrapolation(points[1:2], -f)
    virt_end = endpoint_extrapolation(points[end-1:end], f)

    points_ext = vcat([virt_start], points, [virt_end])

    # println(points_ext)

    fullSpline = Vector{Vector{Float64}}()
    fullSplineHeadings = Vector{Float64}()
    splineSeg = Vector{Vector{Float64}}()
    headingSeg = Vector{Float64}()


    for i in 1:(length(points_ext)-3)
        splineSeg, headingSeg = KB_spline(points_ext[i:i+3], del_d, b, tau)
        # println(splineSeg)
        # println(typeof(splineSeg))
        # println(splineSeg[1:end-1])

        append!(fullSpline, splineSeg[1:end-1])
        append!(fullSplineHeadings, headingSeg[1:end-1])
        # println(fullSpline)
        # println(splineSeg[1:end-1, :])
        # println("----------")
    end

    
    # println(splineSeg[end])
    
    # println(headingSeg[end])

    push!(fullSpline, splineSeg[end])
    push!(fullSplineHeadings, headingSeg[end])

    # println(fullSpline)
    # println(fullSplineHeadings)

    return fullSpline, fullSplineHeadings#vcat(fullSpline...), vcat(fullSplineHeadings...)
end


function endpoint_extrapolation(pts::Vector{Vector{Float64}}, factor::Float64=1.0)
    """ Computes a virtual point for the spline endpoints based on linear extrapolation

    Args:
        pts: numpy array of 2 points, the second of which is the endpoint
        factor: relative length along line of the new point away from the other 2. 1 means the old endpont will be at the midpoint

    Returns:
        endpt: the new virtual endpoint [x,y]
    """

    x = pts[1][1] + (factor + 1) * (pts[2][1] - pts[1][1])
    y = pts[1][2] + (factor + 1) * (pts[2][2] - pts[1][2])

    return [x, y]
end

function KB_spline(pts::Vector{Vector{Float64}}, del_d::Float64=0.5, b::Float64=0.0, tau::Float64=0.0)
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
    M = [
        0 1 0 0;
        -0.5*(1-tau)*(1+b) (1-tau)*b 0.5*(1-tau)*(1-b) 0;
        (1-tau)*(1+b) 0.5*(1-tau)*(1-3*b)-3 -(1-tau)+3 -0.5*(1-tau)*(1-b);
        -0.5*(1-tau)*(1+b) -0.5*(1-tau)*(1-b)+2 0.5*(1-tau)*(1+b)-2 0.5*(1-tau)*(1-b)
    ]

    # println(pts)
    pts = transpose(reduce(hcat, pts))  # gives a matrix, each vector is a row
    # println(pts)

    Cx = M * pts[:,1]
    Cy = M * pts[:,2]

    N_int = 1000
    t_int = range(0, 1, length=N_int)
    T_int = [ones(N_int) 2*t_int 3*t_int.^2]

    Xp = T_int * Cx[2:end]
    Yp = T_int * Cy[2:end]

    L = trapz(t_int, sqrt.(Xp.^2 .+ Yp.^2))
    N = max(Int(round(L / del_d)), 2)

    println("Segment length: $L")
    println("Num pts in Segment: $N")

    t = range(0, 1, length=N)
    T = [ones(N) t t.^2 t.^3]

    X = T * Cx
    Y = T * Cy

    CR_spline = hcat(X, Y)

    Tp = [ones(N) 2 .* t 3 .* t.^2]
    Xp = Tp * Cx[2:end]
    Yp = Tp * Cy[2:end]

    headings = atan.(Yp, Xp)

    CR_spline = [CR_spline[i, :] for i in axes(CR_spline, 1)]

    # println(CR_spline)
    # println(typeof(CR_spline))

    return CR_spline, headings

end
    
function generate_trajectory(points::Vector{Vector{Float64}}, params::NamedTuple)
    v_b =params.v_b
    dt = params.dt
    r = params.r

    del_d = v_b * dt
    fullSpline, fullHeadings = calc_spline(points, del_d, 0.0, 0.0, 0.5)

    N = length(fullHeadings)

    Xref = [[fullSpline[i][1], fullSpline[i][2], fullHeadings[i], 0.0, 0.0] for i in 1:N]
    Uref = [fill(v_b/r, 4) for _ in 1:N]

    return Xref, Uref
end

function main()

    # points = [
    #     0.0 0.0;
    #     4.0 5.0;
    #     -4.0 15.0;
    #     0.0 20.0
    # ]
    points = [
        [0.0, 0.0],
        [4.0, 5.0],
        [-4.0, 15.0],
        [0.0, 20.0]
    ]
    del_d = 0.25
    fullSpline, fullHeadings = calc_spline(points, del_d, 0.0, 0.0, 10)

    # Plot spline
    scatter(fullSpline[:][1], fullSpline[:][2], markersize=2, color=:green, label="Spline")
    scatter!(points[:][1], points[:][2], markersize=4, color=:red, label="Waypoints", legend=:bottomright)
    axis_equal = true
    gui()

    # Plot headings
    plot(fullHeadings, marker=:circle, label="Heading")
    gui()

    #### sample usage: ####
    # points = np.array([[0,0], #starting point
    #         [4, 5],
    #         [-4, 15],
    #         [0, 20]])
    # params = {"v_b":10,
    #           "dt":4,
    #           "r":.375}

    # Xref,Uref = generate_trajectory(points, params) 
end

# main()