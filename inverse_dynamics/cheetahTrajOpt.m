clear;
clc;

tic

data = load('traj_opt.mat'); % load data from liam's traj opt
q_t_meas = data.x(66:100, [22, 36]); % tail's angles - flick starts on frame 66
dq_t_meas = data.dx(66:100, [22, 36]); % tail's angular velocity
q_t_meas_th = q_t_meas(:, 1); % tail's pitch angle measurements
q_t_meas_psi = q_t_meas(:, 2); % tail's yaw angle measurements
q0 = [zeros(3, 1); q_t_meas(1, :)'];
dq0 = [zeros(3, 1); dq_t_meas(1, :)'];

ngrid = 11;
N = length(q_t_meas); % number of samples
fps = 120;
simTime = N/fps;
time_tau = 0:(simTime)/(ngrid-1):simTime; % times at which torques are sampled
time_q = 0:(simTime)/(N-1):simTime; % times at which measurements are sampled
%% Model parameters - Here you put your model stuff.
A = [eye(3), zeros(3, 2)];
v = [10; 2; 0]; % forward velocity obtained from cheetah's initial head velocity
Cd = 1.2;
tau_max = 30;
tau0 = zeros(2*ngrid, 1);

% run the optimization
[tau,fmin,exitflag,output] = trajOpt(tau0, q0, dq0, q_t_meas_th, q_t_meas_psi,...
                                     A, v, Cd, tau_max,...
                                     ngrid, simTime, time_tau, time_q);

tau_th = tau(1:ngrid); 
tau_psi = tau(ngrid+1:end); % use end instead of 2*ngrid?

toc

delete(gcp('nocreate')); % close parallel pool - so that simulink doesn't give a warning for each iteration (makes opt very slow)

% omega_b = [dphi_b*cos(psi_b) + dth_b*cos(phi_b)*sin(psi_b)
%            dth_b*cos(phi_b)*cos(psi_b) - dphi_b*sin(psi_b)
%                       dpsi_b - dth_b*sin(phi_b)           ];

% ang_imp = 

function [xmin,fmin,exitflag,output] = trajOpt(tau0, q0, dq0, q_t_meas_th, q_t_meas_psi,...
                                                  A, v, Cd, tau_max,...
                                                  ngrid, simTime, time_tau, time_q)

mdl = 'cheetahSim';
xLast = []; % Last place compute all was called
simOut = []; % ODE solution structure

%% Define CONSTRAINTS and Bounds. The only variables you will need to optimise are the points for the input pitch and yaw (roll) torques. 
% These values are normalised from -1 to 1 and scaled in the simulation
UB = [ones(ngrid,1); ones(ngrid,1)];
LB = -1*UB;

%% global search - I did this but I coded it for you to just be a normal optimization

opts = optimoptions('fmincon', 'Algorithm', 'sqp',...
                    'display', 'iter',...
                    'UseParallel', true);

prob = createOptimProblem('fmincon', 'x0' ,tau0,...
                          'objective', @objfun,...
                          'lb' ,LB,'ub', UB,...
                          'options',opts);

%% Run Global Search Algorithm

[xmin,fmin,exitflag,output] = fmincon(prob);

    % Cost function
    function cost = objfun(x)
        
        if ~isequal(x,xLast) % Check if computation is necessary

            % unpack the vector of input torques and the time vector            
            tau_th = x(1:ngrid);
            tau_psi = x(ngrid+1:end); % use end instead of 2*ngrid?
            
            % Weird thing in simulink where you need to declare a global
            % variable again
            simTime; time_tau; time_q; q_t_meas_th; q_t_meas_psi;
            A; q0; dq0; v; Cd; tau_max;
            
            % Simulate the model
            simOut = sim(mdl,'SimulationMode','Accelerator','SrcWorkspace', 'current'); %'parent'? or 'current'?
            
            xLast = x;
        end
        cost = simOut.costs(end);
    end
end