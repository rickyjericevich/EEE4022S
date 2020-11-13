%% This script calculates the matrices & equations for cheetah dynamics equation.

%   Where:    M*dqq + C + G = Tau

%% Admin
clear;
clc;

%% constants
g = 9.81;

% tail
M    = 0.66;
L    = 0.75;
I_tb = 0.045; % radius of tail base
R    = 0.09/2; 

% bottom tail 
m1     = 7*M/8;
I1_com = I_tb - M/160*(107/9*L^2 + 3/4*R^2);
I1_com = diag([93*M/280*R^2, I1_com, I1_com]);

% top tail
m2     = M/8;
r      = R/2;
I2_com = 3*M/2560*(L^2 + 4*R^2);
I2_com = diag([3*m2/10*r^2, I2_com, I2_com]);

% drag force stuff
p  = 1.2041; % air density at 20°C
De = 0.08; % effective diam

%% Symbolic defs
disp('Defining symbols')

% https://www.mathworks.com/help/symbolic/add-suffixes-to-symbolic-results.html

syms Cd v_x v_y v_z 'real' 'positive' % drag coeff, velocity of spine in inertial
syms F1_1_x F1_1_y F1_1_z F1_2_x F1_2_y F1_2_z 'real'
syms F2_1_x F2_1_y F2_1_z F2_2_x F2_2_y F2_2_z 'real'
syms tau1_th tau1_psi tau2_th tau2_psi 'real' % tail's pitch & yaw torque
syms th1 psi1 th2 psi2 'real'     % tail angles
syms dth1 dpsi1 dth2 dpsi2 'real'   % derivatives of tail angles
syms ddth1 ddpsi1 ddth2 ddpsi2 'real'

% Generalised coords
q   = [th1; psi1; th2; psi2];
dq  = [dth1; dpsi1; dth2; dpsi2];
ddq = [ddth1; ddpsi1; ddth2; ddpsi2];

% Rotation matrices
R_0__1 = RotZ(psi1) * RotY(th1);
R_1__0 = R_0__1'; % Rotation of bottom tail into body

R_0__2 = RotZ(psi2) * RotY(th2) * R_0__1;
R_2__0 = R_0__2'; % Rotation of top tail into body

%% KINEMATICS - POSITIONS
disp('Kinematics: Positions')

% p_tb = [0; 0; 0]; 

% bottom tail
p1_end__1 = [-L/2; 0; 0];
p1_com__1 = p1_end__1/3; % com of a frustum where r = R/2
pF1_1__1  = p1_end__1/4;
pF1_2__1  = p1_end__1 - p1_end__1/4;

p1_end__0 = R_1__0 * p1_end__1;
p1_com__0 = R_1__0 * p1_com__1;
pF1_1__0  = R_1__0 * pF1_1__1;
pF1_2__0  = R_1__0 * pF1_2__1;

% top tail
p2_end__2 = [-L/2; 0; 0];
p2_com__2 = p2_end__2/4; % com of a cone
pF2_1__2  = p2_com__2;
pF2_2__2  = p2_end__2 - p2_end__2/4;

p2_end__0 = p1_end__0 + R_2__0 * p2_end__2;
p2_com__0 = p1_end__0 + R_2__0 * p2_com__2;
pF2_1__0  = p1_end__0 + R_2__0 * pF2_1__2;
pF2_2__0  = p1_end__0 + R_2__0 * pF2_2__2;

%% KINEMATICS - VELOCITIES
disp('Kinematics: Velocities')

%%% bottom tail
% linear
dp1_end__0 = simplify(jacobian(p1_end__0, q) * dq);
dp1_com__0 = simplify(jacobian(p1_com__0, q) * dq);
dpF1_1__0  = simplify(jacobian(pF1_1__0, q) * dq);
dpF1_2__0  = simplify(jacobian(pF1_2__0, q) * dq);

% angular
% derivative of rotation ratrix
dR_1__0 = sym(zeros(size(R_1__0)));
for i = 1:length(R_1__0)
    dR_1__0(:,i) = jacobian(R_1__0(:,i), q) * dq; 
end
dR_1__0 = simplify(dR_1__0);

% Skew Symmetric Matrix property
ss_omega1 = simplify(R_0__1 * dR_1__0); % in body frame 

omega1 = [ss_omega1(3,2)
           ss_omega1(1,3)
           ss_omega1(2,1)]; 
omega1 = simplify(omega1);

%%% top tail
% linear
dp2_end__0 = simplify(jacobian(p2_end__0, q) * dq);
dp2_com__0 = simplify(jacobian(p2_com__0, q) * dq);
dpF2_1__0  = simplify(jacobian(pF2_1__0, q) * dq);
dpF2_2__0  = simplify(jacobian(pF2_2__0, q) * dq);

% angular
% derivative of rotation ratrix
dR_2__0 = sym(zeros(size(R_2__0)));
for i = 1:length(R_2__0)
    dR_2__0(:,i) = jacobian(R_2__0(:,i), q) * dq; 
end
dR_2__0 = simplify(dR_2__0);

% Skew Symmetric Matrix property
ss_omega2 = simplify(R_0__2 * dR_2__0);

omega2 = [ss_omega2(3,2)
           ss_omega2(1,3)
           ss_omega2(2,1)]; 
omega2 = simplify(omega2);

%% KINETIC ENERGY
disp('Kinetic Energy')

T_trans = simplify(0.5*m1*dp1_com__0'*dp1_com__0 + 0.5*m2*dp2_com__0'*dp2_com__0);
T_rot   = simplify(0.5*omega1'*I1_com*omega1 + 0.5*omega2'*I2_com*omega2);
T_tot   = simplify(T_trans + T_rot);

%% POTENTIAL ENERGY
% disp('Potential Energy')

V = m1*g*p1_com__0(3) + m2*g*p2_com__0(3);
V = simplify(V);

%% Mass Matrix
disp('Mass Matrix')

M = hessian(T_tot, dq);
M = simplify(M);

%% Derivative of Mass Matrix
disp('Derivative of mass matrix')

dM = sym(size(M));
for i = 1:length(M)
    for j = 1:length(M)
        dM(i,j) = simplify(jacobian(M(i,j), q) * dq);
    end
end
dM = simplify(dM);

%% C Matrix -- contains the centrifugal and coriolis accelerations
disp('Coriolis and Centrifugal Matrix')

C = dM * dq - jacobian(T_tot, q)';
C = simplify(C);

%% G Matrix --> Contains the potential energy
% disp('Gravity Matrix')

G = simplify(jacobian(V, q)');

%% B input matrix - Torque acts on relative angles! - Velocity formulation!!!
disp('Input Torques')

B = [tau1_th; tau1_psi; tau2_th; tau2_psi];

%% Q matrix - generalised forces
disp('Generalised Forces')

v = [v_x; v_y; v_z];

% bottom tail's forces & velocity
vF1_1__0 = v + dpF1_1__0;
vF1_2__0 = v + dpF1_2__0;

vF1_1__1 = R_0__1 * vF1_1__0; % try p1_end_1 x omega1?
vF1_2__1 = R_0__1 * vF1_2__0; % try p1_mid_1 x omega1?

vF1_1_proj = [0; vF1_1__1(2); vF1_1__1(3)];
vF1_2_proj = [0; vF1_2__1(2); vF1_2__1(3)];

F1_1__1 = -0.5*p*Cd*L/4*De*norm(vF1_1_proj)*vF1_1_proj;
F1_2__1 = -0.5*p*Cd*L/4*De*norm(vF1_2_proj)*vF1_2_proj;
F1_1__0 = simplify(R_1__0 * F1_1__1);
F1_2__0 = simplify(R_1__0 * F1_2__1);

% top tail's drag forces & velocity
vF2_1__0 = v + dpF2_1__0;
vF2_2__0 = v + dpF2_2__0;

vF2_1__1 = R_0__2 * vF2_1__0; % try p1_end_1 x omega1?
vF2_2__1 = R_0__2 * vF2_2__0; % try p1_mid_1 x omega1?

vF2_1_proj = [0; vF2_1__1(2); vF2_1__1(3)];
vF2_2_proj = [0; vF2_2__1(2); vF2_2__1(3)];

F2_1__2 = -0.5*p*Cd*L/4*De*norm(vF2_1_proj)*vF2_1_proj;
F2_2__2 = -0.5*p*Cd*L/4*De*norm(vF2_2_proj)*vF2_2_proj;
F2_1__0 = simplify(R_2__0 * F2_1__2);
F2_2__0 = simplify(R_2__0 * F2_2__2);

% symbolic drag forces
Fd1_1 = [F1_1_x; F1_1_y; F1_1_z];
Fd1_2 = [F1_2_x; F1_2_y; F1_2_z];

Fd2_1 = [F2_1_x; F2_1_y; F2_1_z];
Fd2_2 = [F2_2_x; F2_2_y; F2_2_z];

Q = [];
for i = length(q)
    Qi1 = simplify(jacobian(pF1_1__0, q)' * Fd1_1);
    Qi2 = simplify(jacobian(pF1_2__0, q)' * Fd1_2);
    Qi3 = simplify(jacobian(pF2_1__0, q)' * Fd2_1);
    Qi4 = simplify(jacobian(pF2_2__0, q)' * Fd2_2);
    Q   = [Q; Qi1 + Qi2 + Qi3 + Qi4];
end
Q = simplify(Q);

%% Constraint Matrix

% A = jacobian(q(1:3), q) % should be [eye(3), zeros(3, 2)]
% dA = sym(size(A));
% for i = 1:length(A)
%     dA(i) = jacobian(A(i), q) * dq;
% end
% dA % should be 0

%% Save to functions

matlabFunction(F1_1__0, F1_2__0, F2_1__0, F2_2__0, 'File', 'tailDragForces',...
               'Vars', {q, dq, Cd, v})

matlabFunction(M, C, G, Q, 'File', 'tailMCGQ',...
               'Vars', {q, dq, [Fd1_1, Fd1_2, Fd2_1, Fd2_2]})

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