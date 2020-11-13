%% run simulink model
% close all;
% clear;
% clc;

% %% initalise variables 
% % box & tail dimensions
% l_b = 0.4; h_b = 0.24; b_b = 0.2;
% l_t = 0.75; r_t = 0.08/2;
% m_b = 15;
% 
% % Body's COM & normal force positions wrt body frame
% p_s = [0; 0; h_b];
% p_b_com__b = [-l_b/2; 0; -h_b/2]; % position of body's com relative to spine
% 
% % tail base in ?? frame, com and tip in tail's frame
% h_tb = 0; % height of tail base ABOVE BODY COM!! 0 <= h_tb <= h_b/2
% p_tb = p_b_com__b + [-l_b/2; 0; h_tb]; 
% p_tt_t = [-l_t; 0; 0];
% % p_com_t = p_tt_t/2;
% 
% %% box setup
% % initial vertices in body frame
% % f = front, b = back    t = top, b = bottom     l = left, r = right
% ftl = [0,-b_b/2, 0]; ftr = [0, b_b/2, 0];
% fbr = ftr + [0, 0, -h_b]; fbl = ftl + [0, 0, -h_b];
% vert_b = [ftl; ftr; fbr; fbl]; % vertices of front face
% vert_b = [vert_b; vert_b - [l_b, 0, 0]]; % append vertices of back face

% %% Run simulink model
% q0  = [0; 0; 0; 0; 0; 0];
% dq0 = zeros(6, 1);
% A   = [0, 0, 0, 0, 1, 0];
% B   = [0; 0; 0; 0; 0; 1];
% v   = 20;
% Cd  = 0;

% out = sim('cheetahSim', 1);
tail_lines = out.tail_line;
vert_I = out.box_verts;
Fd_lines = out.Fd_line;
t = out.tout;
N = length(t);

%% Box faces
front = [1,2,3,4];  back   = [5,6,7,8];
top   = [1,2,6,5];  bottom = [4,3,7,8];
left  = [1,4,8,5];  right  = [2,3,7,6];
f = [front; back; top; bottom; left; right];

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
time =  text(-0.5, 0, 1, 't = 0 s');

hold on
% torso plot setup
plot3(p_s(1), p_s(2), p_s(3), 'o');
% patch([1 -1 -1 1], [1 1 -1 -1], [0 0 0 0], 'g');
% replace with patch as above
% plot floor
[x, y] = meshgrid(-1.5:0.025:1); % Generate x and y data
z = zeros(size(x, 1)); % Generate z data
surf(x, y, z-0.2, 'FaceColor', '#77AC30', 'LineWidth', 0.05) % Plot the surface

p = patch('Faces', f, 'FaceVertexCData', c_c, 'FaceColor', 'interp');
% tail line setup
l = line('LineWidth', 2, 'Color', goldenrod, 'Marker', 'o', 'MarkerSize', 3, 'MarkerFaceColor', goldenrod);
% drag force setup - this deletes other plots so it must come first
q = quiver3(0, 0, 0, 0, 0, 0, 'LineWidth',0.75, 'LineStyle', ':',  'Color', 'r', 'MaxHeadSize', 0.5);
hold off

% clearvars -except tail_lines vert_I Fd_lines t N l p q
idx = 0;
while 1
    idx = mod(idx+1, N) + 1;
    
    p.Vertices = vert_I(:,:,idx);
    
    tail_line = tail_lines(:,:,idx);
    l.XData = tail_line(1,:);
    l.YData = tail_line(2,:);
    l.ZData = tail_line(3,:);
    
    Fd_line = Fd_lines(:,:,idx);
    q.XData = Fd_line(1, 1);
    q.YData = Fd_line(2, 1);
    q.ZData = Fd_line(3, 1);
    q.UData = Fd_line(1, 2);
    q.VData = Fd_line(2, 2);
    q.WData = Fd_line(3, 2);

    axis([-1.5 0.5, -1 1, -0.5 1.2]);
    time.String = ['t = ', num2str(t(idx), '%2.2f'), ' s'];
    drawnow;
%     pause(2*(t(2)-t(1)));
end