% This file is called for InitFcn in simulink model

%% initalise variables 
% box & tail dimensions
l_b = 0.25; h_b = 0.3; b_b = 0.2;
l_t = 0.75; De = 0.08;
m_b = 15; 

% Body's COM & normal force positions wrt body frame
p_s = [0; 0; h_b/2];
p_b_com__b = [-l_b/2; 0; 0]; % position of body's com relative to spine

% tail base in ?? frame, com and tip in tail's frame
p_tb = p_b_com__b + [-l_b/2; 0; 0]; 
p_tt_t = [-l_t; 0; 0];
% p_com_t = p_tt_t/2;

%% box setup
% initial vertices in body frame
% f = front, b = back    t = top, b = bottom     l = left, r = right
ftl = [0,-b_b/2, h_b/2]; ftr = [0, b_b/2, h_b/2];
fbr = ftr + [0, 0, -h_b]; fbl = ftl + [0, 0, -h_b];
vert_b = [ftl; ftr; fbr; fbl]; % vertices of front face
vert_b = [vert_b; vert_b - [l_b, 0, 0]]; % append vertices of back face

%% Run simulink model
q0  = [0; 0; 0; pi/4; 0];
dq0 = zeros(5, 1);
B   = [0; 0; 0; -5; 0];
v   = 0;
Cd  = 1.2;