L = 0.75;

p1_end__1 = [-L/2; 0; 0];

pF1_1__1  = p1_end__1/4;
pF1_2__1  = p1_end__1 - p1_end__1/4;
pF2_1__2  = p1_end__1/4;
pF2_2__2  = p1_end__1 - p1_end__1/4;
pF = [pF1_1__1, pF1_2__1, pF2_1__2, pF2_2__2];

%% Uncomment this if you want to run the sim without optimisation data

% ngrid = 11;
% N = 100; % number of samples
% Ts = 1/120; % sample time
% simTime = N*Ts;
% time_tau = simTime/(ngrid-1)*(0:ngrid-1)'; % times at which torques are sampled
% time = simTime/(N-1)*(0:N-1)'; % times at which measurements are sampled
% 
% q_meas  = zeros(N, 4); % tail's pitch angle measurements
% dq_meas  = zeros(N, 4);
% 
% q0  = zeros(4, 1);
% dq0 = zeros(4, 1);
% 
% v   = [0; 10; 0]; % forward velocity obtained from cheetah's initial head velocity
% Cd  = 1.2;
% tau_max = 30;
% 
% tau1_th  = zeros(ngrid, 1); 
% tau1_psi = zeros(ngrid, 1);
% tau2_th  = zeros(ngrid, 1)/2; 
% tau2_psi = zeros(ngrid, 1);

%% run simulink model
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

tailSim
out = sim('tailSim');

if plotTailOnly
    
    tail_lines = out.tail_lines;
    Fd_lines = out.F_lines;
    t = out.tout;
    N = length(t);

    % Create figure
    figure(1);
    title('Tail', 'FontWeight', 'bold');
    xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)');
    axis('equal'); grid on;
    time_txt =  text(0, 0, 1, 't = 0 s');
    
    % tail colors
    goldenrod = [218, 165, 32]/255; % RGB normalised

    hold on
    % tail lines setup
    l1 = line('LineWidth', 8, 'Color', goldenrod, 'Marker', 'o', 'MarkerSize', 3, 'MarkerFaceColor', goldenrod);
    l2 = line('LineWidth', 8, 'Color', goldenrod, 'MarkerFaceColor', goldenrod); % 'Marker', 'o', 'MarkerSize', 3
    l = [l1, l2];

    % drag force setup
    f1_1 = quiver3(0, 0, 0, 0, 0, 0, 'LineWidth', 1, 'LineStyle', ':',  'Color', 'r');
    f1_2 = quiver3(0, 0, 0, 0, 0, 0, 'LineWidth', 1, 'LineStyle', ':',  'Color', 'r');
    f2_1 = quiver3(0, 0, 0, 0, 0, 0, 'LineWidth', 1, 'LineStyle', ':',  'Color', 'r');
    f2_2 = quiver3(0, 0, 0, 0, 0, 0, 'LineWidth', 1, 'LineStyle', ':',  'Color', 'r');
    f = [f1_1, f1_2, f2_1, f2_2];
    hold off

    idx = 0;
    while 1
        idx = mod(idx, N) + 1;

        tails = tail_lines(:,:,:,idx);
        for i = 1:2
            l(i).XData = tails(1,:,i);
            l(i).YData = tails(2,:,i);
            l(i).ZData = tails(3,:,i);
        end

        forces = Fd_lines(:,:,:,idx);
        for i = 1:4
            f(i).XData = forces(1,1,i);
            f(i).YData = forces(2,1,i);
            f(i).ZData = forces(3,1,i);
            f(i).UData = forces(1,2,i);
            f(i).VData = forces(2,2,i);
            f(i).WData = forces(3,2,i);
        end

        axis([-1 1, -1 1, -1 1]);
        time_txt.String = ['t = ', num2str(t(idx), '%2.2f'), ' s'];
        drawnow;
    %     pause(0.2);
    end
end