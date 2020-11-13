clear;
clc;

mdl = 'tailSim';
load_system(mdl);
set_param([mdl, '/Plotting & To Workspace'], 'Commented', 'on');
save_system;
close_system;

data_path = 'data\2019_03_09 lily flick';
data = load([data_path, '/traj_opt.mat']); % load data from liam's traj opt
% indices = 1:16;
q_meas = data.x(:, [22, 36, 23, 37]); % tail's angles - flick starts on frame 66
dq_meas = data.dx(:, [22, 36, 23, 37]); % tail's angular velocity
angles_meas___I = data.x(:, [18, 4, 32, 19, 5, 33, 20, 21, 7, 35]);
v_meas__I = data.dx(:, [1, 2, 3]); % forward velocity obtained from cheetah's initial head velocity
v = zeros(size(v_meas__I)); % in torso's frame
for i = 1:length(v_meas__I)
    a = angles_meas___I(i, :);
    R_Imeas__head = RotZ(a(3))*RotX(a(2))*RotY(a(1));
    R_head__neck = RotZ(a(6))*RotX(a(5))*RotY(a(4))*R_Imeas__head;
    R_neck__backTorso = RotZ(a(10))*RotX(a(9))*RotY(a(8))*RotY(a(7))*R_head__neck;
    v(i, :) = (R_neck__backTorso * v_meas__I(i, :)')';
end
accel = data.ddx(:, [1, 2, 3]);

ngrid =11;
Ts = 1/90; % sample time
if contains(data_path, '2019')
    Ts = 1/120;
end
N = length(q_meas); % number of samples

ave_accel = (v_meas__I(end) - v_meas__I(1))/(Ts*N);

tau0 = zeros(4*ngrid, 1);
q0 = q_meas(1, :)';
dq0 = dq_meas(1, :)';

simTime = N*Ts;
time_tau = simTime/(ngrid-1)*(0:ngrid-1)'; % times at which torques are sampled
time = simTime/(N-1)*(0:N-1)'; % times at which measurements are sampled
%% Model parameters - Here you put your model stuff.
Cd = 1.2;
tau_max = 30;

% run the optimization
tic
[tau,fmin,exitflag,output] = trajOpt(mdl, tau0, q_meas, dq_meas,...
                                     q0, dq0, v, Cd, tau_max,...
                                     ngrid, simTime, time_tau, time);
toc

tau1_th = tau(1:ngrid);
tau1_psi = tau(ngrid+1:2*ngrid);
tau2_th = tau(2*ngrid+1:3*ngrid);
tau2_psi = tau(3*ngrid+1:end);
                                                                                                
load_system(mdl);
set_param([mdl, '/Plotting & To Workspace'], 'Commented', 'off');
save_system;
close_system;

delete(gcp('nocreate')); % close parallel pool - so that simulink doesn't give a warning for each iteration (makes opt very slow)

plotTailOnly = false;
tailAnim;
if ~plotTailOnly
    cheetahAnim;
end

function [xmin,fmin,exitflag,output] = trajOpt(mdl, tau0, q_meas, dq_meas,...
                                               q0, dq0, v, Cd, tau_max,...
                                               ngrid, simTime, time_tau, time)
    %% Initialize vars
    xLast = []; % Last place compute all was called
    simOut = [];
    
    th1_meas  = q_meas(:, 1);
    psi1_meas = q_meas(:, 2);
    th2_meas  = q_meas(:, 3);
    psi2_meas = q_meas(:, 4);
    
    dth1_meas  = dq_meas(:, 1);
    dpsi1_meas = dq_meas(:, 2);
    dth2_meas  = dq_meas(:, 3);
    dpsi2_meas = dq_meas(:, 4);
    
    vx = v(:, 1);
    vy = v(:, 2);
    vz = v(:, 3);
    %% Bounds on torques
    % These values are normalised from -1 to 1 and scaled in the simulation
    UB = ones(4*ngrid,1);
    LB = -1*UB;

    %% Optimization setup
    opts = optimoptions('fmincon', 'Algorithm', 'sqp',...
                        'MaxFunctionEvaluations', 10000, 'MaxIterations', 5000,...
                        'display', 'iter', 'UseParallel', true);

    prob = createOptimProblem('fmincon', 'x0', tau0, 'objective', @objfun,...
                              'lb', LB, 'ub', UB, 'options', opts);

    %% Optimize!
    [xmin,fmin,exitflag,output] = fmincon(prob);

    %% Function definitions
    function cost = objfun(x)
        if ~isequal(x,xLast) % Check if computation is necessary
            % unpack the vector of input torques           
            tau1_th = x(1:ngrid);
            tau1_psi = x(ngrid+1:2*ngrid);
            tau2_th = x(2*ngrid+1:3*ngrid);
            tau2_psi = x(3*ngrid+1:end);
            
            % Simulate the model
            simOut = sim(mdl,'SimulationMode','Accelerator','SrcWorkspace', 'current');
            
            xLast = x;
        end
        cost = simOut.costs(end);
    end
end
function rot = RotX(th)
    rot = [[1     0       0   ]
           [0  cos(th) sin(th)]
           [0 -sin(th) cos(th)]];
end
% function rot = RotY(th) % standard rotation about y
%     rot = [[ cos(th) 0 sin(th)]
%            [    0    1    0   ]
%            [-sin(th) 0 cos(th)]];
% end
function rot = RotY(th) % Liam's rotation about y
    rot = [[ cos(th) 0 -sin(th)]
           [    0    1    0   ]
           [sin(th) 0 cos(th)]];
end
function rot = RotZ(th)
    rot = [[ cos(th) sin(th) 0]
           [-sin(th) cos(th) 0]
           [    0       0    1]];
end