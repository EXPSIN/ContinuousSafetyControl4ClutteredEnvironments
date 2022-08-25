clear; clc; close all;

opt=optimset('Display','off');

% The map setting
image      = imread('env_02.bmp');
bwimage    = rgb2gray(image) < 125; % obstacle is 1
p_m        = [-0.6; -0.4];          % The start point
resolution = 500;                   % resolution pixels/meter

% Simulation setting
T          = 0.05;          % step
Time       = 7;             % duration
radarrange = 10;            % radar range, meter

% mobile robot state
p          = [-.3; -0.1];   % position, meter
v          = [ 0;  0];      % velocity, m/s
z          = zeros(8, 1);   % actuation state (if exists)
D          = 0.05;          % radius of mobile robot, meter

% shaping parameters
cntReshape = 7;             % the cardinality of positive basis set.
c_A        = cos(2*pi/cntReshape); 

% initialize the graphics
H = graphic_init();

for i = 1:(Time/T)
    t = i*T;
  
    % velocity command for objective
    v_aim = [0.3; 0.3];
  
    % measurement model
    [angle_radar, radius_radar, obs_flag]   = pic2ar(p, radarrange, bwimage, p_m, resolution);
    
    % choose positive basis
    A = [cos(linspace(0, 2*pi-2*pi/cntReshape, cntReshape)); sin(linspace(0, 2*pi-2*pi/cntReshape, cntReshape))]';
    
    % regularized the measurement model.
    r_c     = measurement_regularization(A, angle_radar, radius_radar, c_A);
    alpha_f = 1.0;
    b       = alpha_f*(r_c-D);
    
    % QP-based safety controller
    v_set = quadprog(eye(2), -v_aim, A, b, [],[],[],[],[],opt);
    
    % update model: no acutation uncertainty
    p(:, 1) = rungekutta(@m_position, p(:, 1), v_set, T);     % 位置模型
    
    % update model: exists acutation uncertainty
%     z(:, 1) = rungekutta(@model_v, z(:, 1), v_set, T);
%     C = [1.668   0.05345         0         0    0.1065  -0.02628         0         0;
%              0         0    -3.971  0.001333         0         0   -0.6368     1.575];
%     v(:, 1) = C*z(:, 1);
%     p(:, 1) = rungekutta(@m_position, p(:, 1), v(:, 1), T);
    
    
    %% Update figures
    % velocity command
    H.v_aim.UData = v_aim(1); H.v_aim.VData = v_aim(2); H.v_aim.XData = p(1); H.v_aim.YData = p(2);
    % velocity
    H.v_set.UData = v_set(1); H.v_set.VData = v_set(2); H.v_set.XData = p(1); H.v_set.YData = p(2);
    
    % positive basis
    H.L_vectors.UData = A(:, 1).*b(:, 1);
    H.L_vectors.VData = A(:, 2).*b(:, 1);
    H.L_vectors.XData = p(1)*ones(cntReshape, 1);
    H.L_vectors.YData = p(2)*ones(cntReshape, 1);
    
    % position trajectory
    addpoints(H.traj, p(1), p(2));
        
    % radar
    p_radar = [cos(angle_radar); sin(angle_radar)].*radius_radar;    % 采样结果
    clearpoints(H.radar(1));  addpoints(H.radar(1), p(1)+p_radar(1, :), p(2)+p_radar(2, :));
    H.radarpatch(1).XData = p(1) + p_radar(1, :);
    H.radarpatch(1).YData = p(2) + p_radar(2, :);
    
    % show the regularized measurement model
    L      = [cos(linspace(0, 2*pi, 360)); sin(linspace(0, 2*pi, 360))]';
    radius = measurement_regularization(L, angle_radar, radius_radar, c_A);
    p_radar_polish = L.*radius;
    clearpoints(H.radar_NR);
    addpoints(H.radar_NR, p(1)+p_radar_polish(:, 1), p(2)+p_radar_polish(:, 2));
    
    % The shape of mobile robot
    theta = atan2(v_set(2), v_set(1));
    car = [cos(theta), -sin(theta); sin(theta), cos(theta)] * H.cat_data;
    H.car_patch.XData = car(1, :)+ p(1);
    H.car_patch.YData = car(2, :)+ p(2);
    
    % show the feasible set.
    H.feasibleRegionOriginal = half_space_2D(H.feasibleRegionOriginal,     A(1:end, :),       b(1:end, :), H.axes_range, 'm', 'domain', 'Original', p);
    
    % the velocity trajectory
    addpoints(H.v_x, t, v_set(1));
    addpoints(H.v_y, t, v_set(2));
    
    drawnow limitrate;
end


%{
    p           -- the position of the agent
    map         -- map
    p_m         -- the origin of the map
    resolution  -- specific pix per meter

    Note that the data of the output data is the data in the inertial
    coordinate system.
%}
function [angle, radius, obs_flag] = pic2ar(p, radar_range, pic, p_m, resolution)
map = fliplr(pic');
xN = size(map, 1);
yN = size(map, 2);


cnt_line = 360;
obs_flag = false;
angle    = linspace(2*pi/cnt_line, 2*pi, cnt_line);
radius   = ones(1, cnt_line)*radar_range;
dir      = [cos(angle); sin(angle)];

% Whether the mobile robot is in the obstacle set.
p_pixel  = int32((p-p_m)*resolution);
p_val    = map(min(max(p_pixel(1), 1), xN), min(max(p_pixel(2), 1), yN));
if(p_val == 1 && ((p_pixel(1) >= 1 && p_pixel(1) <= xN && p_pixel(2) >= 1 && p_pixel(2) <= yN)))
    obs_flag = true;
    return;
end

search_radius = linspace(0, radar_range, resolution*radar_range);

for i = 1:cnt_line
    search_points_m     = dir(:, i)*search_radius + (p-p_m);
    search_points_pixel = int32(search_points_m * resolution);
    val = map(min(max(search_points_pixel(1, :), 1), xN) + xN*(min(max(search_points_pixel(2, :), 1), yN)-1) );
    
    map_flag = (search_points_pixel(1, :) >= 1 & search_points_pixel(1, :) <= xN & search_points_pixel(2, :) >= 1 & search_points_pixel(2, :) <= yN);
    index = find(val==1 & map_flag, 1);
    
    if(isempty(index))
        radius(i) = radar_range;
    else
        radius(i) = norm(double([search_points_pixel(1, index); search_points_pixel(2, index)]) - (p - p_m)*resolution) / resolution;
    end
end
end

%{
    draw map
%}
function h = patch_im(pic, p_m, resolution)
% Z = double(fliplr(pic'));
Z = flipud(double(pic));
Z(Z<=0) = -2;
Z(Z> 0) = -1;

xN = size(Z, 2);
yN = size(Z, 1);

x = (1:xN)/resolution+p_m(1);
y = (1:yN)/resolution+p_m(2);
[X, Y] = meshgrid(x, y);
h = surf('XData',X,'YData',Y,'ZData',Z, 'EdgeColor', 'none', 'displayname', 'obstacle', 'HandleVisibility', 'off');
surf('XData',0,'YData',0,'ZData',0, 'FaceColor', 'black', 'displayname', 'obstacle', 'HandleVisibility', 'on');
hold on;
view(2);
% view(-45, 30);
% colorres = [0.9, 0.9, 0.9;0,0,0];
colorres = [ones(1,3);0,0,0];
colormap(colorres);

% h  = patch([0, 1, -1, 0], [1, 3, 3, 1], 'k', 'FaceColor', 'black', 'displayname', 'obstacle', 'HandleVisibility', 'on');
% circle = 0.5*[cos(linspace(0, 2*pi, 21)); sin(linspace(0, 2*pi, 21))];
% h(1)  = patch( 1.0+circle(1, :), circle(2, :), 'k', 'FaceColor', 'black', 'displayname', 'obstacle', 'HandleVisibility', 'on');
% h(2)  = patch(-1.0+circle(1, :), circle(2, :), 'k', 'FaceColor', 'black', 'displayname', 'obstacle', 'HandleVisibility', 'off');
% colormap(1-gray);
end

function dx = model_v(x, u)
A = [-3.251   -0.1069         0         0         0         0         0         0;
     0.0625         0         0         0         0         0         0         0;
          0         0    -32.63    -13.48         0         0         0         0;
          0         0        16         0         0         0         0         0;
          0         0         0         0  -0.01073    -2.489         0         0;
          0         0         0         0         2         0         0         0;
          0         0         0         0         0         0    -8.898      -6.3;
          0         0         0         0         0         0         8         0];
B = [2      0;
     0      0;
     2      0;
     0      0;
     0  0.125;
     0      0;
     0      4;
     0      0];
% C = [1.668   0.05345         0         0    0.1065  -0.02628         0         0;
%          0         0    -3.971  0.001333         0         0   -0.6368     1.575];
dx = A*x + B*u;
end



function x = rungekutta(fun, x0, u, h)
% FcnHandlesUsed  = isa(fun,'function_handle');
k1 = fun(x0       , u);
k2 = fun(x0+h/2*k1, u);
k3 = fun(x0+h/2*k2, u);
k4 = fun(x0+  h*k3, u);
x = x0 + h/6*(k1 + 2*k2 + 2*k3 + k4);
end

function r_c = measurement_regularization(A, angle, radius, c_A)
N     = size(A, 1);
r_c   = zeros(N, 1);
dir   = [cos(angle);sin(angle)];
for i_A = 1:N
    l_i = A(i_A, :);
    r_c(i_A, 1) = min(radius + limit(c_A - l_i*dir, 0, inf));
end
end

function  res = limit(res, lower, upper)
res = min(max(res, lower), upper);
end

function dx = m_position(x, u)
A  = 0*[-1, 0; 0, -1];
B  = [1, 0; 0 1];
dx = A*x + B*u;
end

function H = graphic_init()
Time    = evalin('base','Time');
D       = evalin('base','D');
p       = evalin('base','p');
p_m     = evalin('base','p_m');
bwimage     = evalin('base','bwimage');
resolution  = evalin('base','resolution');
cntReshape  = evalin('base','cntReshape');

figure(2); hold on; set(gcf, 'position', [800, 410, 800, 400]);
H.minh     = animatedline('marker', '.', 'markersize', 5, 'linestyle', '-', 'color', 'r', 'displayname', 'real min');
plot([0,Time], [0,0], 'k--', 'handlevisibility', 'off');
l2=legend;
l2.AutoUpdate = 'off';

figure(1); hold on; 
% 绘制地图
patch_im(bwimage, p_m, resolution);
H.axes_range = [-0.6, 0.6, -0.4, 0.65]*1.01;
% mobile agent
% 
% H.text_primary_velocity_command    = text(0.29, 0.286, {'primary', 'velocity command'}, 'fontsize', 12);
% H.text_actual_velocity             = text(0.266, 0.156, 'actual velocity'            , 'fontsize', 12);
% H.text_positive_basis              = text(-0.07171, -0.0, {'positive', 'basis'}  , 'fontsize', 12, 'HorizontalAlignment', 'right');
% H.text_measured_obstacle_points    = text(-0.58, -0.364, {'measured', 'obstacle points'}   , 'fontsize', 12);
% H.text_regularized_obstacle_points = text( 0.0, -0.286, {'regularized', 'obstacle points'}, 'fontsize', 12);
% H.text_feasible_region             = text(0.02472, -0.1602, {'feasible', 'set'}            , 'fontsize', 12, 'HorizontalAlignment', 'center');
% 
H.L_vectors = quiver(p(1)*ones(cntReshape, 1), p(2)*ones(cntReshape, 1), 0*ones(cntReshape, 1), 0*ones(cntReshape, 1), 1.0, 'color', 'k', 'maxheadsize', 0.5, 'linewidth', 2, 'autoscale', 'off', 'HandleVisibility', 'on', 'displayname', 'Positive Basis');
H.cat_data      = D*[[-0.4; 0.4], [1.0; 0.0], [-0.4; -0.4]];
H.car_patch     = patch('XData',H.cat_data(1, :),'YData',H.cat_data(2, :), 'EdgeColor', 'm', 'linewidth', 2, 'facecolor', 'm', 'facealpha', 0.5, 'displayname', 'Feedback Radar Area', 'HandleVisibility', 'off');
H.radar         = animatedline('marker', 'none', 'markersize',2, 'linestyle', '-', 'linewidth', 2, 'color', 'r', 'displayname', 'Measured Obstacle Points', 'HandleVisibility', 'on');
H.radarpatch    = patch('XData',0,'YData',0, 'EdgeColor', 'none', 'linewidth', 1, 'facecolor', 'r', 'facealpha', 0.0, 'displayname', 'Feedback Radar Area', 'HandleVisibility', 'off');
H.radar_NR      = animatedline('marker', 'none', 'markersize',  5, 'linestyle', '-', 'linewidth', 2, 'color', 'b', 'displayname', 'NR', 'HandleVisibility', 'on');
H.radarpatch_NR = patch('XData',0,'YData',0, 'facecolor', 'b', 'facealpha', 0.0, 'displayname', 'NR', 'HandleVisibility', 'off');

H.traj     = animatedline('marker', 'none', 'markersize', 5, 'linestyle', '-', 'color', 'r', 'HandleVisibility', 'off');

H.feasibleRegionFull     = half_space_2D([], [0, 0], 0, H.axes_range, 'c', 'domain', '', p);
H.feasibleRegionOriginal = half_space_2D([], [0, 0], 0, H.axes_range, 'c', 'domain', '', p);
H.v_aim = quiver(p(1),p(2), 0, 0, 1.0, 'color', 'm', 'maxheadsize', 0.5, 'linewidth', 2, 'autoscale', 'off', 'linestyle', ':', 'HandleVisibility', 'on', 'displayname', 'Primary Velocity Command');
H.v_set = quiver(p(1),p(2), 0, 0, 1.0, 'color', 'm', 'maxheadsize', 0.5, 'linewidth', 2, 'autoscale', 'off', 'HandleVisibility', 'on', 'displayname', 'Velocity Command');




axis equal; axis(H.axes_range);
set(gca, 'box', 'on', 'color', [1,1,1]*1, 'fontsize', 16);
set(gcf, 'position', [10,10,600,450], 'color', 'w');
xlabel('[p]_1 / m');
ylabel('[p]_2 / m');
l1 = legend('location', 'northeast');
% set(l1, 'position', [0.5-0.315/2,0.5-0.1637/2,0.315,0.1637])
l1.AutoUpdate = 'off';
grid off;
axis off;
legend off;

figure(2); hold on; set(gcf, 'position', [800, 0, 800, 500], 'color', 'w'); grid on;
subplot(2, 1, 1);
H.v_x = animatedline('marker', '.', 'markersize', 5, 'linestyle', '-', 'color', 'b', 'displayname', 'sampled min', 'handlevisibility', 'off');
ylabel('v_x [m/s]');set(gca, 'box', 'on', 'color', [1,1,1]*1, 'fontsize', 16);
axis([0, Time, -0.1, 0.2]);
grid on;

subplot(2, 1, 2);
H.v_y = animatedline('marker', '.', 'markersize', 5, 'linestyle', '-', 'color', 'b', 'displayname', 'sampled min', 'handlevisibility', 'off');
ylabel('v_y [m/s]');
xlabel('t [sec]');
axis([0, Time, -0.1, 0.2]); set(gca, 'box', 'on', 'color', [1,1,1]*1, 'fontsize', 16);
grid on;


% figure(3); hold on;
% H.pic          = plot(0, 0, 'k.', 'linewidth', 2, 'markersize', 5, 'linestyle', 'none', 'displayname', 'Obstacle');
% H.radius       = plot(0, 0, 'r', 'linewidth', 2, 'markersize', 5, 'linestyle', '-', 'displayname', 'Measurement Points');
% H.min_radius   = plot(0, 0, 'b', 'linewidth', 2, 'markersize', 5, 'linestyle', '-', 'displayname', 'NR');
% axis([-pi, pi, 0, 8]);
% l3=legend('location', 'northwest');
% l3.AutoUpdate = 'off';
% set(gca, 'fontsize', 16, 'color', 'w', 'position', [0.06, 0.16, 0.93, 0.82], 'box', 'on');
% set(gcf, 'position', [500,10, 800, 400], 'color', 'w');
% xlabel('Angle  rad'); ylabel('Radius [meter]');
% grid on;
% 
% figure(4); hold on; set(gcf, 'position', [800, 610, 800, 400]); grid on;
% H.v_d_norm = animatedline('marker', '.', 'markersize', 5, 'linestyle', '-', 'color', 'b', 'displayname', 'sampled min');
end