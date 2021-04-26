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

tspan = 0:.01:5;
y0 = [2; 0; -pi/6; -.5];
[t,x] = ode45(@(t,x)segway_state_eqns(x,k1, k2, k3, k4, k5, 0),tspan,y0);
% plot(t, x(:, 1))


figure
for q=1:length(t)
    drawcartpend(x(q,:),m,M,2*l);
end


eigs = [-1.1; -1.3;-1.5;-1.7];
K = place(A, B, eigs);  

tspan = 0:0.1:5;
y0 = [2; 0; pi/12; 0];
[t,x] = ode45(@(t,x)segway_state_eqns(x,k1, k2, k3, k4, k5, -K*(x - [-2; 0; 0; 0])),tspan,y0);

