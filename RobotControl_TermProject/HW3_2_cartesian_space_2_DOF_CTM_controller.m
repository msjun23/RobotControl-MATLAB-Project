%%
clc
clear all
close all

%% Set simulation parameters
% Draw flag
flag_sim = 1;
flag_draw = 1;
flag_draw_robot = 1;
flag_draw_graph = 1;

% Global variable
global Iz1 Iz2 L1 L2 g m1 m2 r1 r2 tau1 tau2;

% Simulation parameters
st = 0.000;             % [sec], start time
dt = 0.005;             % [sec], sampling time
ft = 6.000;             % [sec], end time

g = 9.8148;             % [m/s^2], gravitational acceleration

% Robot parameters
m1 = 0.2;       m2 = 0.2;       % [kg], link mass
L1 = 0.5;       L2 = 0.5;       % [m], link length
r1 = 0.1;       r2 = 0.1;       % [m], center of mass
Iz1 = 0.05;     Iz2 = 0.05;     % [kgm^2], link inertia

init_q1 = -pi/2;    init_q2 = pi/2;     % [rad], init joint angle
init_dq1 = 0.00;    init_dq2 = 0.00;    % [rad/s], init angular velocity
q = [init_q1; init_q2];                 % [rad], init joint angle
dq = [init_dq1; init_dq2];              % [rad/s], init angular velocity

init_X = GetKinematics_two_link(q(1), q(2));     % [m], init end-effector position, (2x1)
X = init_X;                             % [m], current position
dX = [0; 0];                            % [m], current Velocity
X_d = init_X;                           % [m], target end-effector position
dX_d = [0; 0];                          % [m/s], target end-effector velocity
ddX_d = [0; 0];                         % [m/s^2], target end-effector acceleration

tau1 = 0.0;             % [Nm], control turque 1
tau2 = 0.0;             % [Nm], control turque 2
tau = [tau1; tau2];     % [Nm], control turque

% Controller gain
Wn = 20;                % [rad/s], natural frequency
Kp = Wn^2;              % propotional gain
Kv = 2*Wn;              % derivative gain
Ki = 0;              % integration gain

%% Simulation
if (flag_sim == 1)
    % Simulation
    n = 1;
    sin_t = 0;
    pre_J = 0;
    for time = st:dt:ft
        % Set target trajectory
        % Print loading message
        cmd = sprintf("loading... %2.2f%%", time/ft*100);
        clc
        disp(cmd);
        
        if (time < 1.0)
            % Waiting for 1.0s
            X_d = init_X;
            dX_d = [0; 0];
            ddX_d = [0; 0];
        elseif (time < 2.0)
            % Move to starting position for 0.5s
            % Stay 0.5s
            X_d(1) = init_X(1);
            if (X_d(2) < init_X(2)+0.1)             % Target pos
                X_d(2) = X_d(2) + (0.1/0.5)*dt;     % Target vel
            else
                X_d(2) = init_X(2) + 0.1;
            end
            dX_d = (X_d - [sim_X_x_d(n-1); sim_X_y_d(n-1)])./dt;
            ddX_d = (dX_d - [sim_dX_x_d(n-1); sim_dX_y_d(n-1)])./dt;
        else
            % Draw a circle every second with the end-effector
            X_d = [0.1*sin((2*pi*sin_t)) + init_X(1);
                   0.1*cos((2*pi*sin_t)) + init_X(2)];
            sin_t = sin_t + dt;
            dX_d = (X_d - [sim_X_x_d(n-1); sim_X_y_d(n-1)])./dt;
            ddX_d = (dX_d - [sim_dX_x_d(n-1); sim_dX_y_d(n-1)])./dt;
        end
        % Get dynamics
        J = GetJacobian_two_link(q(1), q(2));       % (2x2), Get two link Jacobian matrix
        dJ = (J - pre_J)/dt;                        % (2x2), Differential of Jacobian
        pre_J = J;
        
        X = GetKinematics_two_link(q(1), q(2));     % (2x1), Get two link Kinematics
        dX = J*dq;                                  % (2x1)
        D = GetInertia_two_link(q(1), q(2));                % Get Inertia term, D term
        H = GetCoriolis_two_link(q(1), q(2), dq(1), dq(2)); % Get Coriolis term, H term
        G = GetGravity_two_link(q(1), q(2));                % Get Gravity term, G term
        
        % Controller
        u = ddX_d + Kv*(dX_d - dX) + Kp*(X_d - X) + Ki*(X_d - X)*dt;  % (2x1)
        ddq_ref = inv(J)*(u - dJ*dq);               % (2x1)
        tq_ctrl = D*ddq_ref + H + G*0.8;            % (2x1), torque for each link
        
        % Robot model
        % Inverse dynamics
        tau = tq_ctrl;
        tau1 = tau(1);
        tau2 = tau(2);
        
        % Return the dynamics of each link
        [t, y] = ode45('two_link', [0 dt], [q(1); dq(1); q(2); dq(2)]);
        index = length(y);
        
        % Update dynamics of each link
        q = [y(index, 1); y(index, 3)];
        dq = [y(index, 2); y(index, 4)];
        
        % Save data
        sim_time(n) = time;         % [sec]
        sim_q1(n) = q(1);           % [rad]
        sim_q2(n) = q(2);           % [rad]
        sim_dq1(n) = dq(1);         % [rad/s]
        sim_dq2(n) = dq(2);         % [rad/s]
        sim_X_x(n) = X(1);          % [m]
        sim_X_y(n) = X(2);          % [m]
        sim_dX_x(n) = dX(1);        % [m/s]
        sim_dX_y(n) = dX(2);        % [m/s]
        sim_X_x_d(n) = X_d(1);      % [m]
        sim_X_y_d(n) = X_d(2);      % [m]
        sim_dX_x_d(n) = dX_d(1);    % [m/s]
        sim_dX_y_d(n) = dX_d(2);    % [m/s]
        n = n + 1;
    end
end

%% Simulation result graph
if (flag_draw == 1)
    font_size_label = 20;
    font_size_title = 25;
    linewidth_current = 3;
    linewidth_target = 5;
    
    if (flag_draw_robot == 1)
        % Draw robot
        % save as gif
        filename = 'HW3_2_cartesian_sapce_2_DOF_CTM_PID_controller.gif';
        
        x1 = L1*cos(init_q1);               % [m], joint 1 X-axis position
        y1 = L1*sin(init_q1);               % [m], joint 1 Y-axis position
        x2 = L2*cos(init_q1 + init_q2);     % [m], joint 2 X-axis position
        y2 = L2*sin(init_q1 + init_q2);     % [m], joint 2 Y-axis position
        
        FG1 = figure('Color', [1 1 1]);
        AX = axes('parent', FG1);
        hold on;
        
        Px1 = [0 x1];       Py1 = [0 y1];
        Px2 = [x1 x1+x2];   Py2 = [y1 y1+y2];
        
        p1 = plot(Px1, Py1, '-ob', 'Linewidth', linewidth_current);
        p2 = plot(Px2, Py2, '-or', 'Linewidth', linewidth_current);
        
        axis([-0.6 1.0 -1.0 0.4]);
        grid on;
        xlabel('X-axis (m)', 'fontsize', font_size_label);
        ylabel('Y-axis (m)', 'fontsize', font_size_label);
        title('2-DOF Robot', 'fontsize', font_size_title);
        
        n = 1;
        for (time = st:dt:ft)
            % Print run time
            cmd = sprintf("Time: %2.2f", time);
            clc
            disp(cmd);  
            
            q1 = sim_q1(n);                 % [deg], joint 1 angle
            q2 = sim_q2(n);                 % [deg], joint 2 angle
            x1 = L1*cos(q1);                % [m], joint 1 X-axis position
            y1 = L1*sin(q1);                % [m], joint 1 Y-axis position
            x2 = L2*cos(q1+q2);             % [m], joint 2 X-axis position
            y2 = L2*sin(q1+q2);             % [m], joint 2 Y-axis position
            
            Px1 = [0 x1];       Py1 = [0 y1];
            Px2 = [x1 x1+x2];   Py2 = [y1 y1+y2];
            set(p1, 'XData', Px1, 'YData', Py1);
            set(p2, 'XData', Px2, 'YData', Py2);
            drawnow
            n = n + 1;
            
            % save as gif
            frame = getframe(FG1);
            img = frame2im(frame);
            [imind, cm] = rgb2ind(img, 256);
            if time==0
                imwrite(imind, cm, filename, 'gif', 'Loopcount', inf, 'DelayTime', 1/60);
            else
                imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 1/60);
            end
        end
    end
    
    if (flag_draw_graph == 1)
        % Draw position
        FG2 = figure('Color', [1 1 1]);
        plot(sim_time, sim_X_x_d, ':r', 'LineWidth', linewidth_target);
        hold on;
        plot(sim_time, sim_X_y_d, ':b', 'LineWidth', linewidth_target);
        hold on;
        
        plot(sim_time, sim_X_x, 'r', 'LineWidth', linewidth_current);
        hold on;
        plot(sim_time, sim_X_y, 'b', 'LineWidth', linewidth_current);
        hold on;
        
        axis([st ft -1.25 1]);
        grid on;
        
        xlabel('time (s)', 'fontsize', font_size_label);
        ylabel('Position (m)', 'fontsize', font_size_label);
        title('Cartesian Space PID CTM Controller', 'fontsize', font_size_title);
        legend({'tar_x', 'tar_y', 'cur_x', 'cur_y'}, 'location', 'best', 'orientation', 'horizontal', 'fontsize', 15);
        
        % Draw velocity
        FG3 = figure('Color', [1 1 1]);
        plot(sim_time, sim_dX_x_d, ':r', 'LineWidth', linewidth_target);
        hold on;
        plot(sim_time, sim_dX_y_d, ':b', 'LineWidth', linewidth_target);
        hold on;
        
        plot(sim_time, sim_dX_x, 'r', 'LineWidth', linewidth_current);
        hold on;
        plot(sim_time, sim_dX_y, 'b', 'LineWidth', linewidth_current);
        hold on;
        
        axis([st ft -1.25 1.25]);
        grid on;
        
        xlabel('time (s)', 'fontsize', font_size_label);
        ylabel('Velocity (m/s)', 'fontsize', font_size_label);
        title('Cartesian Space PID CTM Controller', 'fontsize', font_size_title);
        legend({'tar_{dx}', 'tar_{dy}', 'cur_{dx}', 'cur_{dy}'}, 'location', 'best', 'orientation', 'horizontal', 'fontsize', 15);
    end
end
