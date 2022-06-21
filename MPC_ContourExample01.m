close all; clear; clc;

% sampling interval
Tz   = 50e-3;

% simulation time
Time = 6;
cnt  = floor(Time/Tz);

% The model of mobile robot.
A  = zeros(2, 2);
B  = eye(2);
Az = eye(size(A)) + A*Tz;
Bz = B*Tz;
Cp = eye(2);

% the state (position) buffer
x  = zeros(2, cnt);

% set the initial position
x(:, 1) = [-1; 0.3];    


% the controller output (velocity) buffer
u  = zeros(2, cnt);

% the obstacles position
po = [[0; -1-1e-3],  [0; 1-1e-3]];

% the radius of obstacles
D  = 1-10e-3;

% figures initialization
H = graphic_init([]);

for k=1:cnt-1
    % Simulation Process

    % simulation time
    t = k*Tz;

    % steps of MPC
    N_MPC = 8;

    % MPC-based safety controller
    [u_mpc, x_mpc, v_max] = MPC_CBF(x(:, k), N_MPC);

    % obtain the output of the safety controller from the solution of MPC
    u(:, k)    = u_mpc(1:size(u, 1));

    % The mobile robot is defined by an integrator
    x(:, k+1)  = model_robot(x(:, k), u(:, k));

    % update figures
    H = graphic_init(H);
end

%% 
%{
    The MPC-CBF function
%}
function [u, x, v_max] = MPC_CBF(x_k, N)
persistent opt opt_nonlinear;
if(isempty(opt))
    opt = optimoptions('quadprog',  'Display','off');
    opt_nonlinear = optimoptions('fmincon', 'Display','off','Algorithm','sqp', 'MaxIterations', 1e5);
end

po    = evalin('base', 'po');
D     = evalin('base', 'D');
Tz    = evalin('base', 'Tz');
Az    = evalin('base', 'Az');
Bz    = evalin('base', 'Bz');
Cp    = evalin('base', 'Cp');

v_max = 1.0;
n_u = size(Bz, 2);
n_x = size(Az, 1);


% Define the solution of MPC-Based controller
u_x = [zeros((n_u)*N, 1); repmat(x_k, N, 1)];

% The primary controller, expecting the mobile robot to move at the speed
% of [1; 0].
u_aim = [repmat([1; 0], [N, 1])];

% cost function
fun = @(u_x) (10*norm(u_x(1:n_u*N)-u_aim)^2 + ... % Tracking the reference control input
    sum(vecnorm(reshape(u_x((n_u+1):(n_u*N)), n_u, N-1) - reshape(u_x(1:n_u*(N-1)), n_u, N-1))) ... % The smoothness cost
    );

% constraints
con = @(u_x) mycon(u_x, n_u, n_x, x_k, po, D, Az, Bz, Cp, N);

% solve the MPC controller
[u_x,fval,exitflag,output,lambda] = fmincon(fun,u_x,[],[],[],[], [], [], con, opt_nonlinear);

% mpc
% if(exitflag == 0)
%     warning('[MPC] Number of iterations exceeded options.MaxIterations or number of function evaluations exceeded options.MaxFunctionEvaluations.');
% end
if(exitflag == -2)
    u_x = [zeros((n_u)*N, 1); repmat(x_k, N, 1)];
    warning('[MPC] infeasible.');
end

u = reshape(u_x(1:(n_u*N), 1), n_u, N);
x = reshape(u_x((n_u*N + 1):end, 1), n_x, N);
end
function [c,ceq] = mycon(u_x, n_u, n_x, x_0, po, D, Az, Bz, Cz, N)
u = reshape(u_x(1:(n_u*N), 1), n_u, N);
x = [x_0, reshape(u_x((n_u*N + 1):end, 1), n_x, N)];

ceq = zeros(n_x*N, 1);
c   = zeros(2*N, 1);
for i = 1:N
    % inequality constrints (safety constraints)
    c(2*(i-1)+1, 1) = -(Cz*x(:, i) - po(:, 1))'/norm(Cz*x(:, i) - po(:, 1))*u(:, i) - (norm(Cz*x(:, i) - po(:, 1))-D);
    c(2*(i-1)+2, 1) = -(Cz*x(:, i) - po(:, 2))'/norm(Cz*x(:, i) - po(:, 2))*u(:, i) - (norm(Cz*x(:, i) - po(:, 2))-D);

    % equality constraints (model constraints)
    ceq((1:n_x)+n_x*(i-1), 1) = -x(:, i+1) + Az*x(:, i) + Bz*u(:, i);
end
end



%{
    The model of mobile robot.
%}
function x = model_robot(x_k, u_k)
Az = evalin('base', 'Az');
Bz = evalin('base', 'Bz');
x = Az*x_k + Bz*u_k;
end

%
function H = graphic_init(H)
if(isempty(H))
    figure(1); hold on; axis equal; grid on;
    set(gcf, 'position', [0,0,600,400], 'color', 'w');

    po   = evalin('base', 'po');
    D    = evalin('base', 'D');
    Time = evalin('base', 'Time');

    H.traj         = animatedline('color', 'k', 'displayname', 'Agent Trajectory');
    H.agent        = plot(0, 0, 'k*', 'linestyle', 'none', 'displayname', 'Mobile Robot', 'handlevisibility', 'off');
    H.position_mpc = plot(0, 0, 'bo', 'linestyle', '-', 'displayname', 'MPC path', 'handlevisibility', 'off');

    H.o_cir       = D*[cos(linspace(0,2*pi,40)); sin(linspace(0,2*pi,40))];
    H.obstacle(1) = patch(po(1,1)+H.o_cir(1, :), po(2,1)+H.o_cir(2, :), 'b', 'facealpha', 0.5,  'displayname', 'Obstacle-1', 'handlevisibility', 'on', 'linestyle', 'none');
    H.obstacle(2) = patch(po(1,2)+H.o_cir(1, :), po(2,2)+H.o_cir(2, :), 'b', 'facealpha', 0.5,  'displayname', 'Obstacle-2', 'handlevisibility', 'on', 'linestyle', 'none');
    H.vel         = quiver(0,0,0,0, 1.0, 'color', 'k', 'LineWidth', 2,  'displayname', 'output of controller', 'handlevisibility', 'on');
    set(gca, 'box', 'on');
    axis([-2, 4, -2,2]);
    H.f1_legend = legend('AutoUpdate','off');


    figure(2);
    set(gcf, 'position', [600,0,600,400], 'color', 'w');
    hold on;
    axis([0, Time, -1.1, 1.1]);
    grid on;
    H.v_traj(1) = animatedline('color', 'r', 'displayname', 'u1', 'marker', '.', 'MarkerSize', 10);
    H.v_traj(2) = animatedline('color', 'b', 'displayname', 'u2');
    H.v_mpc(1)  = plot(0, 0, 'ro', 'linestyle', '-', 'displayname', 'u1 mpc', 'handlevisibility', 'on');
    H.v_mpc(2)  = plot(0, 0, 'bo', 'linestyle', '-', 'displayname', 'u2 mpc', 'handlevisibility', 'on');

    H.f21_legend = legend('AutoUpdate','off','location','northwest');
else
    t     = evalin('base', 't');
    Tz    = evalin('base', 'Tz');
    k     = evalin('base', 'k');
    N_MPC = evalin('base', 'N_MPC');

    x     = evalin('base', 'x');
    x_mpc = evalin('base', 'x_mpc');
    u     = evalin('base', 'u');
    u_mpc = evalin('base', 'u_mpc');

    addpoints(H.traj,  x(1, k), x(2, k));
    H.agent.XData = x(1, k);
    H.agent.YData = x(2, k);
    H.position_mpc.XData = x_mpc(1, :);
    H.position_mpc.YData = x_mpc(2, :);

    addpoints(H.v_traj(1),  t, u(1, k));
    addpoints(H.v_traj(2),  t, u(2, k));
    H.v_mpc(1).XData = t+((1:N_MPC)-1)*Tz;
    H.v_mpc(1).YData = u_mpc(1, :);
    H.v_mpc(2).XData = t+((1:N_MPC)-1)*Tz;
    H.v_mpc(2).YData = u_mpc(2, :);

    H.vel.XData = x_mpc(1, 1);
    H.vel.YData = x_mpc(2, 1);
    H.vel.UData = u_mpc(1, 1);
    H.vel.VData = u_mpc(2, 1);

    drawnow limitrate;
end

end
