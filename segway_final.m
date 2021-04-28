% Project File

% INPUT OF THE SYSTEM
%Given System parameters
M = 3.5;
m = 2;
R = 0.05;
Iw = 0.004375;
Ir = 0.02667;
l = 2;
g = 9.81;

% Forming the constants in the kinematic equations
k1 = (M + m)*R + Iw/R;
k2 = m*l*R;
k3 = m*l;
k4 = Ir + m*l*l;
k5 = m*g*l;

% USER INPUT
initial_theta = pi/9;

% DESIGNER'S CHOICE
% CALIBERATION -> RELATES THE INPUT THETA TO THE OUTPUT CHANGE IN POSITION
initial_x = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x_final = 30*tanh(1.05*initial_theta);

x_f = [x_final; 0; 0; 0];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Open Loop transfer functions

H1 = tf([k2*k3 - k4*k1 0 k1*k5],[k1+k3]);
H2 = tf([-(k2+k4) 0 k5],[k2*k3 - k4*k1 0 k1*k5 0 0]);

H = zpk(H1*H2);
G = tf([-(k2 + k4) 0 k5],[k1+k3 0 0]);


A = [0      1              0                0;
     0      0     -k2*k5/(k1*k4 - k2*k3)    0;
     0      0              0                1;
     0      0     -k1*k5/(-k1*k4 + k2*k3)   0];
B = [     0;
     k2+k4/(k1*k4 - k2*k3);
          0;
        k1+k3/(-k1*k4 + k2*k3)];
C = [1 0 0 0;
     0 1 0 0;
     0 0 1 0;
     0 0 0 1];
D = [0;0;0;0];

states = {'x' 'x_dot' 'theta' 'theta_dot'};
inputs = {'tau'};
outputs = {'x' 'x_dot' 'theta' 'theta_dot'};

sys_ss = ss(A,B,C,D,'statename',states,'inputname',inputs,'outputname',outputs);
csys = canon(sys_ss,'companion');
syms s
Bb = [0;0;0;1];
% det(s*eye(4)-csys.A)
% s^4 - 6.7425*s^2
% eigs = [-1.1; -3.3;-1.5;-0.7];
% s^4 + 17.6s^3 + 43.6s^2 + 360s = 3.6,4.8  10  0 

Cx = ctrb((csys.A)', Bb);
Cz = ctrb(A, B);

% K_x = [8.7 54.9425 57.9 18];
K_x = [8.7 34.9425 37.9 15];
K = K_x*Cx*inv(Cz);
det(s*eye(4) - (csys.A-csys.B*K));

tspan = 0:0.05:1;
y0 = [initial_x; 0; initial_theta; 0];
[t1,x1] = ode45(@(t,x)segway_state_eqns(x,k1, k2, k3, k4, k5, -K*(x - x_f), pi/12, true),tspan,y0);

tspan = 1:0.05:25;
y0 = x1(length(x1),:);
x_final = x1(length(x1),1) + x_f(1);
x_f(1) = x_final;
velocity = 0;
[t2,x2] = ode45(@(t,x)segway_state_eqns(x,k1, k2, k3, k4, k5, -K*(x - x_f), pi/20, false),tspan,y0);

t = cat(1,t1,t2);
x = cat(1,x1,x2);

for q=1:length(t)
    drawcartpend(x(q,:),m,M,2*l);
end