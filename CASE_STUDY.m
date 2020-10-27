
% Title: Improving the Fidelity of Mixed-Monotone Reachable Set 
%        Approximations via State Transformations
% Submitted to American Controls Conference (ACC), 2021
% Author: Matthew Abate and Samuel Coogan

% Code Author: Matthew Abate
% Date: 9/27/2020
% Description:  This script generates Figure 2a.
%               Parallelogram approximations of forward time reachable sets
%               are computed for a given system by applying Theorem 1 and
%               Proposition 1.

clear all; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% System Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global W
W = [-.5, .5]; % disturbance bound

% Transformations and bounds for safe set
global S T2
T1 = eye(2);                    % Trivial Tansforation (rectangle)
T2 = 1/sqrt(2)*[1, -1; 1, 1];   % 45 degree roation    (Rotated rectangle)
S = [-1, 1; ...
     -1, 1];                    % interval bounds for both systems
 
% Observer parameters: observer is circle
Obs_Rad = .2; % radius

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dt = .005;       % time step
Sim_Time = .75;   % total simulation time
Sim_Time_Vec = 0:dt:Sim_Time;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate random signals using GP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
K_ss=kernel(Sim_Time_Vec, Sim_Time_Vec, 1);
CHOL =chol(K_ss + 1e-12*eye(length(Sim_Time_Vec)),'lower');
% Generate 1 random deisturbance signal from GP
gp_signal = (CHOL*normrnd(0, 1, length(Sim_Time_Vec),1))';
const = (W(1, 2)- W(1, 1))/(max(gp_signal) - min(gp_signal));
w_signal(1, :) = const*(gp_signal - min(gp_signal))  + W(1, 1);
% Generate 2 random signals from GP, for observer
% the observer returns a cicle, of radius Obs_Rad
%   D_signal is the distance from true state to the center of observer set
%   A_Signal is the angle between them
gp_signal = (CHOL*normrnd(0, 1, length(Sim_Time_Vec), 1))';
const = (Obs_Rad)/(max(gp_signal) - min(gp_signal));
D_signal(1, :) = const*(gp_signal - min(gp_signal)); % distance signal
K_ss=kernel(Sim_Time_Vec, Sim_Time_Vec, 1);
gp_signal = (CHOL*normrnd(0, 1, length(Sim_Time_Vec), 1))';
const = (2*pi)/(max(gp_signal) - min(gp_signal));
A_signal(1, :) = const*(gp_signal - min(gp_signal)); % angle signal


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check Uncertain Embedded Invariance of Sets
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dx1 = S(1, 1):(S(1, 2) - S(1, 1))/4:S(1, 2);
dx2 = S(2, 1):(S(2, 2) - S(2, 1))/4:S(2, 2);

[Q1, Q2] = meshgrid(dx1, dx2);
thing = [Q1(:)' ; Q2(:)'];
[Q1, Q2] = meshgrid(dx1, dx2);
thing2 = [Q1(:)' ; Q2(:)'];

holder = [];
for i = 1:size(thing, 2)
    i/size(thing, 2)
    for j = 1:size(thing2, 2)
        if prod(thing(:, i) <= thing2(:, j))
            zu = thing(:, i);
            zo = thing2(:, j);
            if prod(-alpha(zu - S(:, 1)) <= d(zu, W(1), zo, W(2))) && ...
               prod(-alpha(zu - S(:, 2)) >= d(zo, W(2), zu, W(1))) && ...
               prod(-alpha(zu - S(:, 1)) <= dT(zu, W(1), zo, W(2))) && ...
               prod(-alpha(zu - S(:, 2)) >= dT(zo, W(2), zu, W(1)))
            else
                error(); % not embdded invariant
            end
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup Plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%
% Figure 1
%%%%%%%%%%%%%%%%
figure(1); clf;
hold on; grid on;


xlabel('$x_1$','Interpreter','latex')
ylabel('$x_2$','Interpreter','latex')
set(gca,'FontSize',16, 'TickLabelInterpreter','latex')
axis([-1.5, 1.5, -1.5, 1.5]);
%xticks([0 1 2 3])
%yticks([0 1 2 3])

Leg = legend();
set(Leg,'visible','off');

grid on;
ax.Layer = 'top';

% Compute vertexes of octagon safe set
S_points = [];
for t = 0:7
    theta = t*pi/4;
    point = 1/cos(pi/8)*[cos(theta + pi/8); sin(theta + pi/8)];
    S_points = [S_points, point];
end
clear t point

% Plot parallelotope 1
%rect1 = T1*makeRectangle(S);
%patch(rect1(1, :), rect1(2, :), 'b', 'FaceAlpha', .2);

% Plot parallelotope 2
%rect2 = T2*makeRectangle(S);
%patch(rect2(1, :), rect2(2, :), 'b', 'FaceAlpha', .2);

% Plot octagon safe set
patch(S_points(1, :), S_points(2, :), 'w', 'FaceAlpha',  .8);
patch(S_points(1, :), S_points(2, :), 'g', 'FaceAlpha', .2);


% circle of radius Obs_Rad centered at 0
circle = Obs_Rad*[cos(0:.1:2*pi); sin(0:.1:2*pi)];


%%%%%%%%%%%%%%%%
% Figure 2
%%%%%%%%%%%%%%%%
figure(2); clf;
hold on; grid on;

subplot(2,1,1)
hold on; grid on;
xlabel('$t$','Interpreter','latex')
ylabel('$u_1$','Interpreter','latex')
set(gca,'FontSize',16, 'TickLabelInterpreter','latex')
axis([0, Sim_Time, -7, 7]);
Leg = legend();
set(Leg,'visible','off');
grid on;
ax.Layer = 'top';

subplot(2,1,2)
hold on; grid on;
xlabel('$t$','Interpreter','latex')
ylabel('$u_2$','Interpreter','latex')
set(gca,'FontSize',16, 'TickLabelInterpreter','latex')
axis([0, Sim_Time, -7, 7]);
Leg = legend();
set(Leg,'visible','off');
grid on;
ax.Layer = 'top';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x0 = [0; 0]; % initial state
x = x0;
x_nom = x0;
u_holder = [];
time_holder = [];
for i = 1:size(Sim_Time_Vec, 2)
    t_now = Sim_Time_Vec(i) % current time
    x_now = x(:, i);        % current state
    w_now = w_signal(1, i); % current disturbance
    
    
    % plot observer set
    dist  = D_signal(1, i);
    angle = A_signal(1, i);
    Obs_now = x_now + dist*[cos(angle); sin(angle)]; % observer set center
    
    % rectangle approximation
    Obs_rect_1 =  [Obs_now, Obs_now] + Obs_Rad*[-1, 1; -1, 1];
    % filter
    Obs_rect_1(:, 1) = max([Obs_rect_1(:, 1), S(:, 1)]')';
    Obs_rect_1(:, 2) = min([Obs_rect_1(:, 2), S(:, 2)]')';
    
    % parallelotope approximation
    Obs_rect_2 =  inv(T2)*[Obs_now, Obs_now] + Obs_Rad*[-1, 1; -1, 1];
    % filter
    Obs_rect_2(:, 1) = max([Obs_rect_2(:, 1), S(:, 1)]')';
    Obs_rect_2(:, 2) = min([Obs_rect_2(:, 2), S(:, 2)]')';
    
    % PLOT OBSERVER
    % plot rectangular approximation of observer set
    %Obs_rect_1_points = makeRectangle(Obs_rect_1);
    %plot_holder(1) = patch(Obs_rect_1_points(1, :), ...
    %                       Obs_rect_1_points(2, :), ...
    %                       'b', 'FaceAlpha', .3);
    % plot parallelotope approximation of observer set
    %Obs_rect_2_points = T2*makeRectangle(Obs_rect_2);
    %plot_holder(end + 1) = patch(Obs_rect_2_points(1, :), ...
    %                             Obs_rect_2_points(2, :), ...
    %                                'b', 'FaceAlpha', .3);
    
    % plot circle observer set
    figure(1);
    plot_holder(1) = patch(Obs_now(1) + circle(1, :), ...
                                 Obs_now(2) + circle(2, :), ...
                                    'w', 'FaceAlpha', .8);
    plot_holder(end + 1) = patch(Obs_now(1) + circle(1, :), ...
                                 Obs_now(2) + circle(2, :), ...
                                    'r', 'FaceAlpha', .3);
    
    % compute next state in trajectory
    u_des = 6*[cos(pi*t_now); sin(pi*t_now)];
    [u_now, time_taken] = ASIF(u_des, Obs_rect_1, Obs_rect_2);
    x_next = x_now + dt*F(x_now, u_now, w_now);
    x(:, i + 1) = x_next;
    
    x_now_nom = x_nom(:, i);
    x_next = x_now_nom + dt*F(x_now_nom, u_des, w_now);
    x_nom(:, i + 1) = x_next;
    
    % plot
    plot_holder(end + 1) = plot(x_nom(1, :), x_nom(2, :), 'm', 'LineWidth', 2);
    plot_holder(end + 1) = scatter(x_nom(1, end), x_nom(2, end), 80, 'm', ...
                                                   'filled');
    plot_holder(end + 1) = plot(x(1, :), x(2, :), 'b', 'LineWidth', 2);
    plot_holder(end + 1) = scatter(x(1, end), x(2, end), 80, 'b', ...
                                                   'filled');
    drawnow
    
    u_holder = [u_holder, [u_des; u_now]];
    
    figure(2); subplot(2,1,1);
    plot_holder2(1) =  plot(0:dt:t_now, u_holder(1, 1:end), 'm', 'LineWidth', 2);
    plot_holder2(2) =  plot(0:dt:t_now, u_holder(3, 1:end), 'b', 'LineWidth', 2);
    drawnow
    figure(2); subplot(2,1,2);
    plot_holder3(1) =  plot(0:dt:t_now, u_holder(2, 1:end), 'm', 'LineWidth', 2);
    plot_holder3(2) =  plot(0:dt:t_now, u_holder(4, 1:end), 'b', 'LineWidth', 2);
    drawnow
    
    pause(.0001)
    figure(1); delete(plot_holder);
    figure(2); subplot(2,1,1); delete(plot_holder2);
    figure(2); subplot(2,1,2); delete(plot_holder3);
    
    time_holder = [time_holder, time_taken];
end


figure(1)
step = 2;
plot_holder(1) = patch(Obs_now(1) + circle(1, :), ...
                                 Obs_now(2) + circle(2, :), ...
                                    'w', 'FaceAlpha', .8);
plot_holder(end + 1) = patch(Obs_now(1) + circle(1, :), ...
                                 Obs_now(2) + circle(2, :), ...
                                    'r', 'FaceAlpha', .3);
plot_holder(end + 1) = plot(x_nom(1, 1:step:end), x_nom(2, 1:step:end), 'r', 'LineWidth', 2);
plot_holder(end + 1) = scatter(x_nom(1, end), x_nom(2, end), 80, 'r', ...
                                               'filled');
plot_holder(end + 1) = plot(x(1, 1:step:end), x(2, 1:step:end), 'b', 'LineWidth', 2);
plot_holder(end + 1) = scatter(x(1, end), x(2, end), 80, 'b', ...
                                               'filled');
drawnow;
figure(2); subplot(2,1,1);
plot_holder2(1) =  plot(0:step*dt:t_now, u_holder(1, 1:step:end), 'r', 'LineWidth', 2);
plot_holder2(2) =  plot(0:step*dt:t_now, u_holder(3, 1:step:end), 'b', 'LineWidth', 2);
figure(2); subplot(2,1,2);
plot_holder3(1) =  plot(0:step*dt:t_now, u_holder(2, 1:step:end), 'r', 'LineWidth', 2);
plot_holder3(2) =  plot(0:step*dt:t_now, u_holder(4, 1:step:end), 'b', 'LineWidth', 2);

% average time per computation
av_time = mean(time_holder)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out = makeRectangle(X0)
    d = [X0(1, 2) - X0(1, 1); X0(2, 2) - X0(2, 1)];
    [X0_x, X0_y] = meshgrid(X0(1, 1): d(1)/5 :X0(1, 2), ...
                            X0(2, 1): d(2)/5 :X0(2, 2));
    X_int = [X0_x(:), X0_y(:)];
    [k,av] = convhull(X_int);
    out = X_int(k, :)';
end

function [out, time_taken] = ASIF(ud, rect, par)
    % ASIF for MM systems
    %   ud   = candidate input
    %   rect = filtered rectangle
    %   par  = filtered parallelotope
    global S W T2
    
    A = [ -eye(2); eye(2); -inv(T2); inv(T2)];
    
    b = [  alpha(rect(:, 1) - S(:, 1)) + d(rect(:, 1), W(1), rect(:, 2), W(2)); ...
         - alpha(rect(:, 2) - S(:, 2)) - d(rect(:, 2), W(2), rect(:, 1), W(1)); ...
           alpha( par(:, 1) - S(:, 1)) + dT(par(:, 1), W(1), par(:, 2), W(2)); ...
         - alpha( par(:, 2) - S(:, 2)) - dT(par(:, 2), W(2), par(:, 1), W(1))];
    
    options =  optimset('Display','off');
    tic
    out = quadprog(eye(2), -ud, A,b, [], [], [], [], ud,options);
    time_taken = toc;
end

function out = alpha(x)
    out = 5000*x.^3;
end

%%%%%%%%%%%%%%%%%%%%
% Nominal Dynamics
%%%%%%%%%%%%%%%%%%%%

function out = F(x, u, w)
    out(1, 1) = - x(1) - x(1)^3 + x(2)  + w^3 + u(1);
    out(2, 1) = - x(2) - x(2)^3 - x(1)  - w + u(2);
end

function out = d(x, w, xhat, what)
    out(1, 1) = - x(1) - x(1)^3 + x(2)  + w^3;
    out(2, 1) = - x(2) - x(2)^3 - xhat(1)  - what;
end

%%%%%%%%%%%%%%%%%%%%
% Transformed Dynamics
%%%%%%%%%%%%%%%%%%%%

function out = dydt(y, w)
    global T2
    out = inv(T2)*F(T2*y, zeros(2, 1), w);
end

function out = dT(y, w, yhat, what)
    if y(1) <= yhat(1) && y(2) <= yhat(2) && w <= what
        % compute minimum of F1 y2, w 
        
        fun_y1 = @(x) [1, 0] * dydt([y(1); x], 0);
        fun_y2 = @(x) [0, 1] * dydt([x; y(2)], 0);
        fun_w1 = @(v) [1, 0] * dydt([0; 0], v);
        fun_w2 = @(v) [0, 1] * dydt([0; 0], v);
        
        options = optimset('TolX', 1e-10);
        out_y1 = fminbnd(fun_y1, y(2), yhat(2), options);
        out_y2 = fminbnd(fun_y2, y(1), yhat(1), options);
        out_w1 = fminbnd(fun_w1, w, what, options);
        out_w2 = fminbnd(fun_w2, w, what, options);
        
        out(1, 1) = fun_y1(out_y1) + fun_w1(out_w1);
        out(2, 1) = fun_y2(out_y2) + fun_w2(out_w2);
        
        out = out - dydt([0; 0], 0);
        
    elseif yhat(1) <= y(1) && yhat(2) <= y(2) && what <= w
        fun_y1 = @(x) -[1, 0]* dydt([y(1); x], 0);
        fun_y2 = @(x) -[0, 1]* dydt([x; y(2)], 0);
        fun_w1 = @(v) -[1, 0]* dydt([0; 0], v);
        fun_w2 = @(v) -[0, 1]* dydt([0; 0], v);
        
        options = optimset('TolX', 1e-10);
        out_y1 = fminbnd(fun_y1, yhat(2), y(2), options);
        out_y2 = fminbnd(fun_y2, yhat(1), y(1), options);
        out_w1 = fminbnd(fun_w1, what, w, options);
        out_w2 = fminbnd(fun_w2, what, w, options);
        
        out(1, 1) = - fun_y1(out_y1) - fun_w1(out_w1);
        out(2, 1) = - fun_y2(out_y2) - fun_w2(out_w2);
        out = out - dydt([0; 0], 0);
    else
        error
    end
end
