clear;
clc;
tic
if gcp('nocreate')==0
    parpool;
end

ngrid = 11;
simTime = 0.3;
tau0 = zeros(2*ngrid, 1);

% load_system('cheetahSim4trajOpt');
% set_param('cheetahSim4trajOpt/Plotting&ToWorkspace','commented','on');
% close_system('cheetahSim4trajOpt', 1, 'OverwriteIfChangedOnDisk', true);
% set_param('cheetahSim4trajOpt/Scopes','commented','on');

[tau,fmin,exitflag,output] = trajOpt(tau0);

tau_th = tau(1:ngrid);
tau_psi = tau(ngrid+1:2*ngrid); % use end instead of 2*ngrid?
timeVec = 0:(simTime)/(ngrid-1):simTime;

% delete(gcp('nocreate'));
toc

function [xming,fming,exitflag,output] = trajOpt(tau0)

mdl = 'cheetahSim4trajOpt';
xLast = []; % Last place compute all was called
simOut = []; % ODE solution structure
ngrid = 11; % Grid size of control. You can probably make this '11'
simTime = 0.3;  %Check the simulation to see how fast the flick occurs. I'm guessing it's probably much faster than this.
timeVec = 0:(simTime)/(ngrid-1):simTime;
%% Model parameters - Here you put your model stuff.
q0  = [-0.16; 0.3; 0.27; -1; 0.2];
dq0 = [-1; -2.5; -1.8; -6; -11];
% q0  = [0; pi/18; 0; 0; -pi/4];
% dq0 = [0; 0; 0; 0; 0];
v   = [10; 2; 0];
Cd  = 1.2;
tau_max = 30;  % Check my Aerodynamics paper for this number.
%% Define CONSTRAINTS and Bounds. The only variables you will need to optimise are the points for the input pitch and yaw (roll) torques. 
% These values are normalised from -1 to 1 and scaled in the simulation
% Aineq = []; Bineq = []; Aeq = []; Beq = [];
UB = [ones(ngrid,1); ones(ngrid,1)];
LB = -1*UB;

%% global search - I did this but I coded it for you to just be a normal optimization

opts = optimoptions('fmincon', 'Algorithm', 'sqp',...
                    'display', 'iter',...
                    'UseParallel', true);

prob = createOptimProblem('fmincon', 'x0' ,tau0,...
                          'objective', @objfun,'nonlcon', @constr,...
                          'lb' ,LB,'ub', UB,...
                          'options',opts);

%% Run Global Search Algorithm

%[xming,fming,flagg,outptg,manyminsg] = run(gs,problem);
[xming,fming,exitflag,output] = fmincon(prob);

    % Cost function
    function cost = objfun(x)
        
        if ~isequal(x,xLast) % Check if computation is necessary

            % Pack the vector of input torques and the time vector            
            tau_th = x(1:ngrid);
            tau_psi = x(ngrid+1:end); % use end instead of 2*ngrid?
            
            % Weird thing in simulink where you need to declare a global
            % variable again
            q0; dq0; v; Cd; tau_max;
            
            % Simulate the model
            simOut = sim(mdl,'SimulationMode','Accelerator','SrcWorkspace', 'current'); %'parent'? or 'current'?
            
            xLast = x;
        end
        cost = simOut.costs(end);
%         cost = cost(end);
    end

    function [cineq,ceq] = constr(x)
                
        if ~isequal(x,xLast) % Check if computation is necessary
            
            simTime;
            
            tau_th = x(1:ngrid);
            tau_psi = x(ngrid+1:end);
            
            % Weird thing in simulink where you need to declare a global
            % variable again
            q0; dq0; v; Cd; tau_max;
            
            % Simulate the model
            simOut = sim(mdl,'SimulationMode','Accelerator', 'on','SrcWorkspace', 'current'); %'parent'? or 'current'?
        
            xLast = x;
        end
        
        q  = simOut.q;

        % Get angles as a vector.
%         th_b = q(:,1);
        phi_b = q(:,2);
%         psi_b = q(:,3);
        th_t = q(:,4);
        psi_t = q(:,5);

        % Inequality Constraints - body roll to above 90 or less than -90, you
        % should put it as an inequlaity
        cineq = [abs(phi_b)-pi/2; -th_t-pi/3; th_t-3*pi/4; abs(psi_t)-pi/2.1];
        
        % Equality constraints - final value of angles must be (0,0,0)
%         ceq = [th_b(end); phi_b(end); psi_b(end)];
        ceq = [];
    end
end