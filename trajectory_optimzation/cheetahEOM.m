%% This script calculates the matrices & equations for cheetah dynamics equation.

%   Where:    M*dqq + C + G = Tau

%% Admin
clear;
clc;

%% constants
% g = 9.81;
% back torso
m_b = 15;
l_b = 0.25;
h_b = 0.3;
b_b = 0.2;
% tail
m_t = 0.66;
l_t = 0.75;
I_tb = 0.045;
d_com = 0.2; % distance from com (from thesis) to tail base

% drag force stuff
p = 1.2041; % air density at 20°C
% Cd = 1.2;
De = 0.08; % effective diam

%% Symbolic defs
disp('Defining symbols')

% https://www.mathworks.com/help/symbolic/add-suffixes-to-symbolic-results.html

syms Cd v_x v_y v_z 'real' 'positive' % drag coeff, velocity of spine in inertial
syms F_d_x F_d_y F_d_z 'real'
syms tau_th tau_psi 'real' % tail's pitch & yaw torque
syms th_b phi_b psi_b 'real'     % body angles
syms dth_b dphi_b dpsi_b 'real'   % derivatives of body angles
syms ddth_b ddphi_b ddpsi_b 'real'
syms th_t psi_t 'real'     % tail angles
syms dth_t dpsi_t 'real'   % derivatives of tail angles
syms ddth_t ddpsi_t 'real'

% Generalised coords
q   = [th_b; phi_b; psi_b; th_t; psi_t];
dq  = [dth_b; dphi_b; dpsi_b; dth_t; dpsi_t];
ddq = [ddth_b; ddphi_b; ddpsi_b; ddth_t; ddpsi_t];

% Rotation matrices
R_I__b = RotZ(psi_b) * RotX(phi_b) * RotY(th_b);
R_b__I = R_I__b'; % Rotation of body into inertial
R_I__t = RotZ(psi_t) * RotY(th_t) * R_I__b;
R_t__I = R_I__t'; % Rotation of tail into inertial

%% KINEMATICS - POSITIONS
disp('Kinematics: Positions')

% Body's COM & normal force positions wrt body frame
p_s = [0; 0; 0]; % spine's position in inertial for now
p_b_com__b = [-l_b/2; 0; 0]; % position of body's com relative to spine

p_b_com__I = simplify(p_s + R_b__I*p_b_com__b); % position of body's com in inertial

% tail base relative to com, com and tip in tail's frame
p_tb = p_b_com__b + [-l_b/2; 0; 0]; 
p_tt__t = [-l_t; 0; 0];
p_t_com__t = p_tt__t/2; % or /4

p_tb__I = simplify(p_s + R_b__I*p_tb);
p_t_com__I = simplify(p_tb__I + R_t__I*p_t_com__t);

%% KINEMATICS - VELOCITIES
disp('Kinematics: Velocities')

%%% body
disp(' - Body')
% linear
dp_b_com__I = simplify(jacobian(p_b_com__I, q) * dq);

% angular
% derivative of rotation ratrix
dR_b__I = sym(zeros(size(R_b__I)));
for i = 1:length(R_b__I)
    dR_b__I(:,i) = jacobian(R_b__I(:,i), q) * dq; 
end
dR_b__I = simplify(dR_b__I);

% Skew Symetric Matrix property to express angular rotation in inertial
ss_omega_b = simplify(R_b__I' * dR_b__I);

omega_b = [ss_omega_b(3,2)
           ss_omega_b(1,3)
           ss_omega_b(2,1)]; 
omega_b = simplify(omega_b, 'IgnoreAnalyticConstraints', true);

%%% tail
disp(' - Tail')
% linear
dp_t_com__I = simplify(jacobian(p_t_com__I, q) * dq);

% angular
% derivative of rotation ratrix
dR_t__I = sym(zeros(size(R_t__I)));
for i = 1:length(R_t__I)
    dR_t__I(:,i) = jacobian(R_t__I(:,i), q) * dq; 
end
dR_t__I = simplify(dR_t__I);

% Skew Symetric Matrix property
ss_omega_t = simplify(R_t__I' * dR_t__I);

omega_t = [ss_omega_t(3,2)
           ss_omega_t(1,3)
           ss_omega_t(2,1)]; 
omega_t = simplify(omega_t, 'IgnoreAnalyticConstraints', true);

%% KINETIC ENERGY
disp('Kinetic Energy')
% Moment of inertia for torso (modelled as a box)
I_b_x = m_b * (1/12 * (h_b^2 + b_b^2));
I_b_y = m_b * (1/12 * (l_b^2 + h_b^2));
I_b_z = m_b * (1/12 * (l_b^2 + b_b^2));
I_b = diag([I_b_x I_b_y I_b_z]);

% Moment of inertia for tail
I_t_com = I_tb - m_t*d_com^2; % I around com in the thesis - shift I to be around this model's com?
I_t = diag([0.001 I_t_com I_t_com]);

T_trans = simplify(0.5*m_b*dp_b_com__I'*dp_b_com__I + 0.5*m_t*dp_t_com__I'*dp_t_com__I);
T_rot   = simplify(0.5*omega_b'*I_b*omega_b + 0.5*omega_t'*I_t*omega_t);
T_tot   = simplify(T_trans + T_rot);

%% POTENTIAL ENERGY
% disp('Potential Energy')

% V_b   = m_b*g*p_b_com__I(3) + 0.5*k*(p_b_l__I(3)^2 + p_b_r__I(3)^2);
% V_t   =  m_t * g * p_t_com__I(3);
% V_tot = simplify(V_b + V_t);

%% Mass Matrix
disp('Mass Matrix')

M = hessian(T_tot, dq);
% M = simplify(M) % This takes too long

%% Derivative of Mass Matrix
disp('Derivative of mass matrix')

dM = sym(size(M));
for i = 1:length(M)
    for j = 1:length(M)
        dM(i,j) = jacobian(M(i,j), q) * dq;
    end
end
% dM = simplify(dM); % This takes too long

%% C Matrix -- contains the centrifugal and coriolis accelerations
disp('Coriolis and Centrifugal Matrix')

C = dM * dq - jacobian(T_tot, q)';
% C = simplify(C) % This takes too long

%% G Matrix --> Contains the potential energy
% disp('Gravity Matrix')

% G = simplify(jacobian(V_tot, q)')

%% B input matrix - Torque acts on relative angles! - Velocity formulation!!!
disp('Input Torques')

B = [0; 0; 0; tau_th; 0; tau_psi];

%% Q matrix - generalised forces
disp('Generalised Forces')

v = [v_x; v_y; v_z];

Vtot__I = v + dp_t_com__I; % tail's total velocity in inertial
Vtot__t = R_I__t * Vtot__I; % veloctiy in tail's frame
V_proj = [0; Vtot__t(2); Vtot__t(3)]; % projection onto YZ plane
Fd = simplify(-0.5*p*Cd*l_t*De*norm(V_proj)*V_proj);

F_d = [F_d_x; F_d_y; F_d_z];
Q = jacobian(p_t_com__I, q)' * F_d;
%Q = simplify(Q) % This takes too long

%% Constraint Matrix

% A = jacobian(phi_t, q) % should be [0, 0, 0, 0, 1, 0]
% dA = sym(size(A));
% for i = 1:length(M)
%     dA(i) = jacobian(A(i), q) * dq;
% end
% dA % should be 0

%% Save to functions

matlabFunction(Fd, 'File', 'Fdrag',...
    'Vars', {[th_b; phi_b; psi_b; th_t; psi_t],...
            [dth_b; dphi_b; dpsi_b; dth_t; dpsi_t],...
            Cd, v})

matlabFunction(M, C, Q, 'File', 'MCQ',...
    'Vars', {[th_b; phi_b; psi_b; th_t; psi_t],...
            [dth_b; dphi_b; dpsi_b; dth_t; dpsi_t],...
            F_d})

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