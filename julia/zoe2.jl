function zoe2_dynamics(model::NamedTuple, x, u)
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

    # Compute individual turning radii for each wheel using trigonometric relations
    R_fl = (L - B*sin(θ_f)) / (2*sin(θ_f) + ϵ)  # Front left wheel radius
    R_fr = (L + B*sin(θ_f)) / (2*sin(θ_f) + ϵ)  # Front right wheel radius
    R_rl = (L - B*sin(θ_r)) / (2*sin(θ_r) + ϵ)  # Rear left wheel radius
    R_rr = (L + B*sin(θ_r)) / (2*sin(θ_r) + ϵ)  # Rear right wheel radius

    # Compute the forward velocity v_b as an average contribution from all 4 wheels
    # Each wheel’s contribution is based on the relationship v = r*ω and the local turning radius
    v_b = (L*r)/(8) * (1/(tan(θ_f) + ϵ) + 1/(tan(θ_r) + ϵ)) * (
        (ω_fl*sin(θ_f))/(L - B*sin(θ_f) + ϵ) +
        (ω_fr*sin(θ_f))/(L + B*sin(θ_f) + ϵ) +
        (ω_rl*sin(θ_r))/(L - B*sin(θ_r) + ϵ) +
        (ω_rr*sin(θ_r))/(L + B*sin(θ_r) + ϵ)
    )

    # Compute the rover's heading rate
    # The heading rate is v_b divided by the body turning radius
    ψ_dot = 4*v_b / (L * (1/(tan(θ_f) + ϵ) + 1/(tan(θ_r) + ϵ)))

    # Compute the steering angle rates from the difference in wheel speeds on each axle
    θ_f_dot = (r / B) * (ω_fr - ω_fl)
    θ_r_dot = (r / B) * (ω_rl - ω_rr)

    # Assemble the state derivative vector
    x_dot = [
        v_b * cos(ψ);
        v_b * sin(ψ);
        ψ_dot;
        θ_f_dot;
        θ_r_dot;
    ]

    return x_dot
end

function animate_zoe2(Xsim, dt)
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
    
    # Build the animation from the joint states
    anim = Animation(vis, fps=floor(Int, 1/dt))
    
    # Set the position and orientation of the rover for each frame
    for k = 1:length(Xsim)
        atframe(anim, k) do

            # Extract the state variables from the current configuration
            x_b = Xsim[k][1]
            y_b = Xsim[k][2]
            ψ   = Xsim[k][3]
            θ_f = Xsim[k][4]
            θ_r = Xsim[k][5]

            # Set the transformations for the rover's position, heading, and steering angles
            T = Translation(x_b, y_b, 0.0)
            settransform!(vis[:robot], T)

            # Handle the heading rotation
            # The rover's heading is defined by the yaw angle ψ
            # We need to rotate the rover's body around the Z-axis by ψ
            # The rotation matrix for a yaw rotation is:
            # R = [cos(ψ) -sin(ψ) 0;
            #      sin(ψ)  cos(ψ) 0;
            #      0        0      1]
            # We can use the LinearMap function to create this rotation
            # R = LinearMap([cos(ψ), -sin(ψ), 0.0;
            #                sin(ψ),  cos(ψ), 0.0;
            #                0.0,     0.0,    1.0])
            # settransform!(vis[:robot], compose(T, R))

            # Build the full 7-DOF configuration vector.
            # Here we assume:
            #   Joint 1 (axle_roll_back_joint): 0.0 (unused)
            #   Joint 2 (axle_yaw_front_joint): θ_f
            #   Joint 3 (axle_yaw_back_joint):  θ_r
            #   Joints 4-7 (wheel joints): 0.0 (no wheel rotation simulated)
            q = [0.0, θ_f, θ_r, 0.0, 0.0, 0.0, 0.0]
            set_configuration!(mvis, q)
        end
    end
    
    # Pass the animation to the visualizer.
    setanimation!(vis, anim)
    
    # Render the visualizer (for example, in a Jupyter notebook).
    return render(vis)
end



