function [xn,yn]=simulate_next_step(x0,dt,fv,fv2,cam_params,qc,MASK)
    % This function simulates a time step in the simulator.
    % Inputs:
    % x0   :  initial state of simulator;
    % dt   :  time for propagation [dt];
    % fv   :  triangulation object of detailed model;
    % fv2  :  triangulation object of reduced model;
    % npix :  number of pixels on vertical side of sensor;
    % Hf   :  horizontal size of sensor [m];
    % Vf   :  vertical size of sensor [m];
    % foc  :  focal length for camera [m];
    % b    :  baseline between cameras along y [m];
    %
    % Outputs:
    % xn   : vector of new state after dt;
    % yn   : vector of size [nh,m] of m measurements of size nh each;

    %% simulate next state:
    % simulates process to obtain next state after dt:
        xn=sim_states(x0,dt);   

    %% simulate measures at new state:
    % simulates acquired measures after dt:
        [yn]=sim_measure(xn,fv,fv2,cam_params,qc,MASK);

end