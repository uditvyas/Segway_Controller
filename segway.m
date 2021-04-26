% Project File

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

initial_theta = pi/12;
initial_x = -5;

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

sys_tf = tf(sys_ss);

eigs = [-1.1; -3.3;-1.5;-0.7];
K = place(A, B, eigs);

tspan = 0:0.05:2;
y0 = [initial_x; 0; initial_theta; 0];
velocity = 0;
[t1,x1] = ode45(@(t,x)segway_state_eqns(x,k1, k2, k3, k4, k5, -K*(x - [initial_x + initial_theta*12 + velocity; 0; 0; 0]), pi/20, true),tspan,y0);
x1(:, 4) = 0;

tspan = 2:0.05:20;
initial_x = x1(length(x1), 1);
y0 = x1(length(x1), :);
y0(4) = 0;
velocity = x1(length(x1),2);
[t,x] = ode45(@(t,x)segway_state_eqns(x,k1, k2, k3, k4, k5, -K*(x - [initial_x + initial_theta*12 + velocity; 0; 0; 0]), pi/100, false),tspan,y0);
final_x = cat(1, x1, x);
time = cat(1,t1, t);
% 
figure
for q=1:length(time)
    drawcartpend(final_x(q,:),m,M,2*l);
end
