function zoe2_dynamics(model::NamedTuple, x, u; debug=false)
    """
    Compute the dynamics of the rover using a kinematic model.
    The rover is modeled as a 4-wheeled vehicle with front and rear steering.

    Args:
        The state vector is defined as:
        x = [x_b, y_b, ψ, θ_f, θ_r]
        where:
            x_b: x position of the rover body (world frame)
            y_b: y position of the rover body (world frame)
            ψ: heading of the rover body (world frame)
            θ_f: front steering angle (body frame)
            θ_r: rear steering angle (body frame)

        The control inputs are the angular velocities of the wheels:
        u = [ω_fl, ω_fr, ω_rl, ω_rr]
        where:
            ω_fl: front left wheel angular velocity
            ω_fr: front right wheel angular velocity
            ω_rl: rear left wheel angular velocity
            ω_rr: rear right wheel angular velocity

    Returns:
        The function returns the state derivative vector:
        x_dot = [x_dot, y_dot, ψ_dot, θ_f_dot, θ_r_dot]
        where:
            x_dot: derivative of x position
            y_dot: derivative of y position
            ψ_dot: derivative of heading
            θ_f_dot: derivative of front steering angle
            θ_r_dot: derivative of rear steering angle
    """

    # Extract state variables
    x_b   = x[1]         # x position of the rover body (world frame)
    y_b   = x[2]         # y position of the rover body (world frame)
    ψ     = x[3]         # heading of the rover body (world frame)
    θ_f   = x[4]         # front steering angle (body frame)
    θ_r   = x[5]         # rear steering angle (body frame)

    # Extract model parameters
    L = model.L          # Robot length (steering joint to steering joint)
    B = model.B          # Axle width
    r = model.r          # Wheel radius

    # Extract control inputs (wheel angular velocities)
    ω_fl = u[1]          # Front left wheel angular velocity
    ω_fr = u[2]          # Front right wheel angular velocity
    ω_rl = u[3]          # Rear left wheel angular velocity
    ω_rr = u[4]          # Rear right wheel angular velocity

    # Compute the turning radius of the rover body using an average of the front and rear steering
    ϵ = 1e-6  # Small value to avoid division by zero
    R_b = L/4 * (1/(tan(θ_f) + ϵ) + 1/(tan(θ_r) + ϵ))
    if debug
        println("R_b: ", R_b)
    end

    ϕ_f = (π/2) - θ_f
    ϕ_r = (π/2) - θ_r
    R_b = L/4 * (tan(ϕ_f) + tan(ϕ_r))
    if debug
        println("Reparameterized R_b: ", R_b)
    end

    # Compute individual turning radii for each wheel using trigonometric relations
    R_fl = (L - B*sin(θ_f)) / (2*sin(θ_f) + ϵ)  # Front left wheel radius
    R_fr = (L + B*sin(θ_f)) / (2*sin(θ_f) + ϵ)  # Front right wheel radius
    R_rl = (L - B*sin(θ_r)) / (2*sin(θ_r) + ϵ)  # Rear left wheel radius
    R_rr = (L + B*sin(θ_r)) / (2*sin(θ_r) + ϵ)  # Rear right wheel radius
    if debug
        println("R_fl: ", R_fl)
        println("R_fr: ", R_fr)
        println("R_rl: ", R_rl)
        println("R_rr: ", R_rr)
    end

    # Compute individual turning radii for each wheel using trigonometric relations
    R_fl = (L - B*cos(ϕ_f)) / (2*cos(ϕ_f))  # Front left wheel radius
    R_fr = (L + B*cos(ϕ_f)) / (2*cos(ϕ_f))  # Front right wheel radius
    R_rl = (L - B*cos(ϕ_r)) / (2*cos(ϕ_r))  # Rear left wheel radius
    R_rr = (L + B*cos(ϕ_r)) / (2*cos(ϕ_r))  # Rear right wheel radius
    if debug
        println("Reparameterized R_fl: ", R_fl)
        println("Reparameterized R_fr: ", R_fr)
        println("Reparameterized R_rl: ", R_rl)
        println("Reparameterized R_rr: ", R_rr)
    end

    # Compute the forward velocity v_b as an average contribution from all 4 wheels
    # Each wheel’s contribution is based on the relationship v = r*ω and the local turning radius
    v_b = (L*r)/(8) * (1/(tan(θ_f) + ϵ) + 1/(tan(θ_r) + ϵ)) * (
        (ω_fl*sin(θ_f))/(L - B*sin(θ_f) + ϵ) +
        (ω_fr*sin(θ_f))/(L + B*sin(θ_f) + ϵ) +
        (ω_rl*sin(θ_r))/(L - B*sin(θ_r) + ϵ) +
        (ω_rr*sin(θ_r))/(L + B*sin(θ_r) + ϵ)
    )
    if debug
        println("v_b: ", v_b)
    end

    # Compute the forward velocity v_b as an average contribution from all 4 wheels
    # Each wheel’s contribution is based on the relationship v = r*ω and the local turning radius
    v_b = (L*r)/(8) * (tan(ϕ_f) + tan(ϕ_r)) * (
        (ω_fl*cos(ϕ_f))/(L - B*cos(ϕ_f)) +
        (ω_fr*cos(ϕ_f))/(L + B*cos(ϕ_f)) +
        (ω_rl*cos(ϕ_r))/(L - B*cos(ϕ_r)) +
        (ω_rr*cos(ϕ_r))/(L + B*cos(ϕ_r))
    )
    if debug
        println("Reparameterized v_b: ", v_b)
    end


    # Compute the rover's heading rate
    # The heading rate is v_b divided by the body turning radius
    ψ_dot = 4*v_b / (L * (1/(tan(θ_f) + ϵ) + 1/(tan(θ_r) + ϵ)))
    if debug
        println("ψ_dot: ", ψ_dot)
    end

    # Compute the rover's heading rate
    # The heading rate is v_b divided by the body turning radius
    ψ_dot = 4*v_b / (L * (tan(ϕ_f) + tan(ϕ_r)))
    if debug
        println("Reparameterized ψ_dot: ", ψ_dot)
    end

    # Compute the steering angle rates from the difference in wheel speeds on each axle
    θ_f_dot = (r / B) * (ω_fr - ω_fl)
    θ_r_dot = (r / B) * (ω_rl - ω_rr)
    if debug
        println("θ_f_dot: ", θ_f_dot)
        println("θ_r_dot: ", θ_r_dot)
    end

    # Compute the steering angle rates from the difference in wheel speeds on each axle
    ϕ_f_dot = (r / B) * (ω_fr - ω_fl)
    ϕ_r_dot = (r / B) * (ω_rl - ω_rr)
    if debug
        println("Reparameterized ϕ_f_dot: ", ϕ_f_dot)
        println("Reparameterized ϕ_r_dot: ", ϕ_r_dot)
    end

    θ_f_dot = -ϕ_f_dot
    θ_r_dot = -ϕ_r_dot
    if debug
        println("Reparameterized θ_f_dot: ", θ_f_dot)
        println("Reparameterized θ_r_dot: ", θ_r_dot)
    end

    # Compute the derivative of the rover's position (y is forward motion)
    x_dot = v_b * -sin(ψ) # have to flip to get positive x on the right
    y_dot = v_b * cos(ψ)
    if debug
        println("x_dot: ", x_dot)
        println("y_dot: ", y_dot)
    end

    # Robot Orientation:

    # y |
    #   |
    #   |
    #   |_______ x

    #  ψ = 0                  ψ > 0
    # |-----|                \-----\
    #    |                       \
    #    |                        \
    #    |                         \
    # |-----|                    \-----\

    # Steering:

    # θ_f, θ_r = 0          θ_f, θ_r > 0
    #   |-----|               \-----\
    #      |                     |
    #      |                     |
    #      |                     |
    #   |-----|               /-----/


    # Assemble the state derivative vector
    x_dot = [
        x_dot;
        y_dot;
        ψ_dot;
        θ_f_dot;
        θ_r_dot;
    ]

    return x_dot
end

function zero_torque!(torques::AbstractVector, t, state::MechanismState)
    torques .= 0
end

function animate_zoe2(Xsim, dt; heading_offset=π/2, Xref=nothing)
    """
    Animate the rover along the trajectory defined by Xsim.
    
    Args:
        Xsim is a vector of vectors, where each vector is a state/configuration:
        [x_b, y_b, ψ, θ_f, θ_r].

        dt is the time step between successive states in the animation.
    """
    
    # Define the path to your rover's URDF file
    urdf_path = joinpath(@__DIR__, "zoe2.urdf")
    
    # Create the MeshCat visualizer
    vis = Visualizer()
    
    # Parse the URDF file to create the robot model
    robot_obj = parse_urdf(urdf_path)
    
    # Create a mechanism visualizer to visualize the robot
    mvis = MechanismVisualizer(robot_obj, URDFVisuals(urdf_path), vis[:robot])

    # Create and set the initial state (all zeros)
    init_state = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    set_configuration!(mvis, init_state)

    # Create the state (with a little motion) including the velocity
    # (have to give a little so that the robot keeps updating)
    init_velocity = [1e-4, 1e-4, 1e-4, 1e-4, 1e-4, 1e-4, 1e-4]
    state = MechanismState(robot_obj, init_state, init_velocity)
    # println("dt: ", dt)
    # println("Xsim length: ", size(Xsim))

    # Simulate the rover up to the length of Xsim with delta time dt
    # t = collect(0:dt:(length(Xsim)-1)*dt)
    # q = [zeros(7) for _ in 1:length(Xsim)]
    t, q, v = simulate(state, (length(Xsim)-1) * dt, zero_torque!, Δt=dt)
    # println("q length: ", size(q))
    println("t length: ", length(t))

    # Add dots along the trajectory
    for k in 1:length(Xsim)
        pt = Xsim[k]
        dot = Sphere(Point3f(0, 0, 0.5), 0.05)
        node = vis[Symbol("dot_$k")]
        mat = MeshPhongMaterial(color = RGBA(0, 0.5, 1.0, 1.0))
        setobject!(node, dot, mat) 
        settransform!(node, Translation(pt[1], pt[2], 0.0))
    end

    # If a reference trajectory is provided, add it to the Animation
    if Xref !== nothing
        for k in 1:length(Xref)
            pt = Xref[k]
            dot = Sphere(Point3f(0, 0, 0.5), 0.05)
            node = vis[Symbol("dot_ref_$k")]
            mat = MeshPhongMaterial(color = RGBA(1.0, 0.5, 0.0, 1.0))
            setobject!(node, dot, mat) 
            settransform!(node, Translation(pt[1], pt[2], 0.0))
        end
    end
    
    # Overwrite the simulation with our own joint states
    for k in 1:length(Xsim)
        q[k][2] = Xsim[k][4]   # front axle steering
        q[k][3] = -Xsim[k][5]  # back axle steering (flipped)
    end
    # println("q[1]: ", q[1])

    # Create a new animation object with the time vector and joint states
    # fps = 1 / dt
    anim = MeshCat.Animation(mvis, t, q; fps=Int.(1/dt))
    # anim = MeshCat.Animation(mvis, t, q; fps=60)

    # Now animate the rover: update both the global transform and joint configuration per frame.
    for k = 1:length(Xsim)
        atframe(anim, k) do
            # Extract the state variables from the current configuration
            x_b = Xsim[k][1]
            y_b = Xsim[k][2]
            ψ   = Xsim[k][3] - heading_offset
            θ_f = Xsim[k][4]
            θ_r = Xsim[k][5]
            
            # --- Update Global Transform ---
            # Build a rotation about the z-axis:
            R = [ cos(ψ)  -sin(ψ)  0.0;
                  sin(ψ)   cos(ψ)  0.0;
                  0.0      0.0     1.0 ]
            # Build the translation
            T = Translation(x_b, y_b, 0.0)
            # Compose the overall transform
            TR = compose(T, LinearMap(R))
            settransform!(vis[:robot], TR)
        end
    end
    
    # Pass the animation to the visualizer.
    setanimation!(mvis, anim)
    
    # Render the visualizer (for example, in a Jupyter notebook).
    return render(vis)
end



