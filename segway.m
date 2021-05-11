clear all
clc

%% INPUT OF THE SYSTEM

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
k6 = M*l*R;
tpc = 6;    % Torque Proportionality Constant (Refer to system modelling for more information)

% USER INPUT
initial_x = 0;
initial_theta = pi/12;
holding_time = 3;

%% STATE SPACE MODEL

A = [0      1              0                0;
     0      0     -k2*k5/(k1*k4 - k2*k3)    0;
     0      0              0                1;
     0      0     -k1*k5/(-k1*k4 + k2*k3)   0];
B = [     0;
     (tpc*k2+k4)/(k1*k4 - k2*k3);
          0;
        (tpc*k1+k3)/(-k1*k4 + k2*k3)];
C = [1 0 0 0;
     0 1 0 0;
     0 0 1 0;
     0 0 0 1];
D = [0;0;0;0];

states = {'x' 'x_dot' 'theta' 'theta_dot'};
inputs = {'tau'};
outputs = {'x' 'x_dot' 'theta' 'theta_dot'};

sys_ss = ss(A,B,C,D,'statename',states,'inputname',inputs,'outputname',outputs);
[b,a] = ss2tf(A,B,C,D);

% State Space to Transfer Functions
H1 = tf(b(1,:),a);
H2 = tf(b(2,:),a);
H3 = tf(b(3,:),a);
H4 = tf(b(4,:),a);


%% OPEN-LOOP ANALYSIS
y0 = [0; 0;initial_theta ; 0];
tspan = 0:0.05:10;
[t0,x0] = ode45(@(t,x)segway_state_eqns(x,k1, k2, k3, k4, k5, k6, tpc, 1, initial_theta, false, false),tspan,y0);

axes( 'Position', [0, 0.95, 1, 0.05] ) ;
set( gca, 'Color', 'None', 'XColor', 'None', 'YColor', 'None' ) ;
text( 0.5, 0, 'Non-Linear System Open-Loop Step Response', 'FontSize', 14', 'FontWeight', 'Bold', ...
      'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom' ) ;

subplot(2,2,1); plot(t0, x0(:,1)); grid on; title("Position"); xlabel("Time (in s)"); ylabel("Position(in m)");
subplot(2,2,2); plot(t0, x0(:,2)); grid on; title("Velocity"); xlabel("Time (in s)"); ylabel("Velocity(in m/s)");
subplot(2,2,3); plot(t0, x0(:,3)); grid on; title("Angle"); xlabel("Time (in s)"); ylabel("Angle of Tilt(in radians)");
subplot(2,2,4); plot(t0, x0(:,4)); grid on; title("Angular Velocity"); xlabel("Time (in s)"); ylabel("Angular Velocity(in radians/s)");
%% Controllability Matrices and gain Values

csys = canon(sys_ss,'companion');
Bb = [0;0;0;1];

% Controllablitity Matrices for Conversion to Controllable Canonical Form
Cx = ctrb((csys.A)', Bb);
Cz = ctrb(A, B);

%% Linearized System

% Holding is True, Controller 1 Design
K_x = [9 15 61 110];
K = K_x*Cx*inv(Cz);
K(1) = 0;

tspan = 0:0.05:holding_time;
v_steady= 10*tanh(1.05*initial_theta);

% Initial values, and required steady state value (x is not affected as the K vector is 0)
x_f = [0; v_steady; initial_theta; 0];
y0 = [0; 0;initial_theta ; 0];

[t1,x1] = ode45(@(t,x)segway_state_eqns(x,k1, k2, k3, k4, k5, k6, tpc, -K*(x_f-x), initial_theta, true, true),tspan,y0);

% Holding is False, Controller 2 Design
K_x = [9 30 38 15];
K = K_x*Cx*inv(Cz);

tspan_1= holding_time:0.05:20;
y0 = x1(length(x1),:);
v_int = x1(length(x1),2);

% Choosing x_final as a linear function of the current position and velocity
x_final = 3.2*v_int + x1(length(x1),1);
x_f = [x_final;0;0;0];
[t2,x2] = ode45(@(t,x)segway_state_eqns(x,k1, k2, k3, k4, k5, k6, tpc, -K*(x -x_f), nan, false, true),tspan_1,y0);

% Complete Graphical Modelling - Linear
t_linear = cat(1,t1,t2);
x_linear = cat(1,x1,x2);

figure
axes( 'Position', [0, 0.95, 1, 0.05] ) ;
set( gca, 'Color', 'None', 'XColor', 'None', 'YColor', 'None' ) ;
text( 0.5, 0, 'Linearised System Closed-Loop Response', 'FontSize', 14', 'FontWeight', 'Bold', ...
      'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom' ) ;

subplot(2,2,1); plot(t_linear, x_linear(:,1));grid on; title("Position"); xlabel("Time (in s)"); ylabel("Position(in m)");
subplot(2,2,2); plot(t_linear, x_linear(:,2));grid on; title("Velocity"); xlabel("Time (in s)"); ylabel("Velocity(in m/s)");
subplot(2,2,3); plot(t_linear, x_linear(:,3));grid on; title("Angle"); xlabel("Time (in s)"); ylabel("Angle of Tilt(in radians)");
subplot(2,2,4); plot(t_linear, x_linear(:,4));grid on; title("Angular Velocity"); xlabel("Time (in s)"); ylabel("Angular Velocity(in radians/s)");

%% Non-Linear Model with Gain Coefficients of Linear System

% Holding is True, Controller 1 Design
K_x = [9 15 61 110];
K = K_x*Cx*inv(Cz);
K(1) = 0;

tspan = 0:0.05:holding_time;
v_steady= 10*tanh(1.05*initial_theta);

% Initial values, and required steady state value (x is not affected as the K vector is 0)
x_f = [0; v_steady; initial_theta; 0];
y0 = [0; 0;initial_theta ; 0];

[t1,x1] = ode45(@(t,x)segway_state_eqns(x,k1, k2, k3, k4, k5, k6, tpc, -K*(x_f-x), initial_theta, true, false),tspan,y0);

% Holding is False, Controller 2 Design
K_x = [9 30 38 15];
K = K_x*Cx*inv(Cz);

tspan_1= holding_time:0.05:20;
y0 = x1(length(x1),:);
v_int = x1(length(x1),2);

% Choosing x_final as a linear function of the current position and velocity
x_final = 3.2*v_int + x1(length(x1),1);
x_f = [x_final;0;0;0];
[t2,x2] = ode45(@(t,x)segway_state_eqns(x,k1, k2, k3, k4, k5, k6, tpc, -K*(x -x_f), nan, false, false),tspan_1,y0);

% Complete Graphical Modelling - Linear
t_nonlinear = cat(1,t1,t2);
x_nonlinear = cat(1,x1,x2);

figure
axes( 'Position', [0, 0.95, 1, 0.05] ) ;
set( gca, 'Color', 'None', 'XColor', 'None', 'YColor', 'None' ) ;
text( 0.5, 0, 'Non-Linear System Closed-Loop Response', 'FontSize', 14', 'FontWeight', 'Bold', ...
      'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom' ) ;

subplot(2,2,1); plot(t_nonlinear, x_nonlinear(:,1));grid on; title("Position"); xlabel("Time (in s)"); ylabel("Position(in m)");
subplot(2,2,2); plot(t_nonlinear, x_nonlinear(:,2));grid on; title("Velocity"); xlabel("Time (in s)"); ylabel("Velocity(in m/s)");
subplot(2,2,3); plot(t_nonlinear, x_nonlinear(:,3));grid on; title("Angle"); xlabel("Time (in s)"); ylabel("Angle of Tilt(in radians)");
subplot(2,2,4); plot(t_nonlinear, x_nonlinear(:,4));grid on; title("Angular Velocity"); xlabel("Time (in s)"); ylabel("Angular Velocity(in radians/s)");

%% Comparison between Linear and Non-linear Plots
figure
axes( 'Position', [0, 0.95, 1, 0.05] ) ;
set( gca, 'Color', 'None', 'XColor', 'None', 'YColor', 'None' ) ;
text( 0.5, 0, 'Comparison between Linear and Non-Linear Systems', 'FontSize', 14', 'FontWeight', 'Bold', ...
      'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom' ) ;

subplot(2,2,1); plot(t_linear, x_linear(:,1));grid on; hold on; title("Position"); xlabel("Time (in s)"); ylabel("Position(in m)");
subplot(2,2,2); plot(t_linear, x_linear(:,2));grid on; hold on; title("Velocity"); xlabel("Time (in s)"); ylabel("Velocity(in m/s)");
subplot(2,2,3); plot(t_linear, x_linear(:,3));grid on; hold on; title("Angle"); xlabel("Time (in s)"); ylabel("Angle of Tilt(in radians)");
subplot(2,2,4); plot(t_linear, x_linear(:,4));grid on; hold on; title("Angular Velocity"); xlabel("Time (in s)"); ylabel("Angular Velocity(in radians/s)");

subplot(2,2,1); plot(t_nonlinear, x_nonlinear(:,1));grid on; title("Position"); xlabel("Time (in s)"); ylabel("Position(in m)");
subplot(2,2,2); plot(t_nonlinear, x_nonlinear(:,2));grid on; title("Velocity"); xlabel("Time (in s)"); ylabel("Velocity(in m/s)");
subplot(2,2,3); plot(t_nonlinear, x_nonlinear(:,3));grid on; title("Angle"); xlabel("Time (in s)"); ylabel("Angle of Tilt(in radians)");
subplot(2,2,4); plot(t_nonlinear, x_nonlinear(:,4));grid on; title("Angular Velocity"); xlabel("Time (in s)"); ylabel("Angular Velocity(in radians/s)");

Lgnd = legend({"Linear System", "Non-Linear System"});
Lgnd.Position(1) = 0;
Lgnd.Position(2) = 0.48;

%% Simulation
for q=1:length(t_nonlinear)
    drawcartpend(x_nonlinear(q,:),m,M,2*l);
end