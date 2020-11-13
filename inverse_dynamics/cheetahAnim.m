%% initalise positions
% box & tail dimensions
l_b = 0.25; h_b = 0.3; b_b = 0.2;
l_t = 0.75; De = 0.08;
m_b = 15;

% Body's COM & normal force positions wrt body frame
p_s = [0; 0; 0.2 + h_b];
p_tb = [-l_b; 0; 0]; % rel to spine
p1_end__1 = [-l_t/2; 0; 0];

pF1_1__1  = p1_end__1/4;
pF1_2__1  = p1_end__1 - p1_end__1/4;
pF2_1__2  = pF1_1__1;
pF2_2__2  = p1_end__1 - p1_end__1/4;
p_F = [pF1_1__1, pF1_2__1, pF2_1__2, pF2_2__2];

% box vertices in body frame
% f = front, b = back    t = top, b = bottom     l = left, r = right
ftl = [0,-b_b/2, 0]; ftr = [0, b_b/2, 0];
fbr = ftr + [0, 0, -h_b]; fbl = ftl + [0, 0, -h_b];
vert_b = [ftl; ftr; fbr; fbl]; % vertices of front face
vert_b = [vert_b; vert_b - [l_b, 0, 0]]; % append vertices of back face

%% Uncomment this if you want to run the sim without optimisation data

% ngrid = 11;
% N = 100; % number of samples
% Ts = 1/120; % sample time
% simTime = N*Ts;
% time_tau = simTime/(ngrid-1)*(0:ngrid-1)';
% time = simTime/(N-1)*(0:N-1)';
% 
% q0  = zeros(7, 1);
% dq0 = zeros(7, 1);
% 
% v = zeros(N, 3);
% Cd  = 1.2;
% tau_max = 30;
% 
% tau1_th  = zeros(ngrid, 1); 
% tau1_psi = zeros(ngrid, 1);
% tau2_th  = zeros(ngrid, 1); 
% tau2_psi = zeros(ngrid, 1);

%% run simulink model
if length(q0)~=7
    q0 = [zeros(3, 1); q0];
    dq0 = [zeros(3, 1); dq0];
end

vx = v(:, 1);
vy = v(:, 2);
vz = v(:, 3);

A = [eye(3), zeros(3, 4)];

cheetahSim
% if ~exist('Sim')
    Sim = sim('cheetahSim','SimulationMode','Accelerator');
% end
vert_I = Sim.box_verts;
tail_lines = Sim.tail_lines;
Fd_lines = Sim.F_lines;
t = Sim.tout;
N = length(t);

angImp = Sim.impulse(end, :)

%% Box faces
front = [1,2,3,4];  back   = [5,6,7,8];
top   = [1,2,6,5];  bottom = [4,3,7,8];
left  = [1,4,8,5];  right  = [2,3,7,6];
faces = [front; back; top; bottom; left; right];

% box vertex colors
goldenrod =   [218, 165, 32]/255; % RGB normalised
saddlebrown = [139, 69, 19]/255; % RGB normalised
map = goldenrod(ones(1, 8), :);
map([1, 7], :) = saddlebrown(ones(1, 2), :);
c_c = colormap(map);

%% Create figure
figure(1);
title({'Cheetah'}, 'FontWeight', 'bold');
xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)');
axis('equal'); grid on;
time_txt =  text(0, 0, 1, 't = 0 s');

hold on
% torso plot setup
plot3(p_s(1), p_s(2), p_s(3), 'o');
% patch([1 -1 -1 1], [1 1 -1 -1], [0 0 0 0], 'g');
% replace with patch as above
% plot floor
[x, y] = meshgrid(-1:0.025:1); % Generate x and y data
z = zeros(size(x, 1)); % Generate z data
surf(x, y, z, 'FaceColor', '#77AC30', 'LineWidth', 0.05) % Plot the surface

p = patch('Faces', faces, 'FaceVertexCData', c_c, 'FaceColor', 'interp');

% tail lines setup
l1 = line('LineWidth', 8, 'Color', goldenrod, 'Marker', 'o', 'MarkerSize', 3, 'MarkerFaceColor', goldenrod);
l2 = line('LineWidth', 8, 'Color', goldenrod, 'MarkerFaceColor', goldenrod); % 'Marker', 'o', 'MarkerSize', 3
l = [l1, l2];

% drag force setup
f1_1 = quiver3(0, 0, 0, 0, 0, 0, 'LineWidth', 1, 'LineStyle', ':',  'Color', 'r','MaxHeadSize', 0.5);
f1_2 = quiver3(0, 0, 0, 0, 0, 0, 'LineWidth', 1, 'LineStyle', ':',  'Color', 'r','MaxHeadSize', 0.5);
f2_1 = quiver3(0, 0, 0, 0, 0, 0, 'LineWidth', 1, 'LineStyle', ':',  'Color', 'r','MaxHeadSize', 0.5);
f2_2 = quiver3(0, 0, 0, 0, 0, 0, 'LineWidth', 1, 'LineStyle', ':',  'Color', 'r','MaxHeadSize', 0.5);
f = [f1_1, f1_2, f2_1, f2_2];
hold off

% save([data_path, '/WorkspaceVars.mat'])

idx = 0;
while 1
    idx = mod(idx+10, N) + 1;
    
    p.Vertices = vert_I(:,:,idx);
    
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

    axis([-1 0.5, -1 1, -0.2 1.2]);
    time_txt.String = ['t = ', num2str(t(idx), '%2.2f'), ' s'];
    drawnow;
    if t(idx)>0.17
        t=t;
    end
%     pause(2*(t(2)-t(1)));
end