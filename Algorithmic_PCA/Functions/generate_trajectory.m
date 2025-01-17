% generate trajectory
function traj = generate_trajectory(env, n_steps)
    traj = HasselmoForage(env, n_steps);
    traj = round(traj);
end