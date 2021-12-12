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
global m I L g tau;

% Simulation parameters
dt = 0.005;             % [sec], sampling time
st = 0.000;             % [sec], start time
ft = 5.000;             % [sec], end time

g = 9.8148;             % [m/s^2], gravitational acceleration

% Robot parameters
m = 1.0000;             % [kg], link mass
L = 1.0000;             % [m], link length
I = (m*L^2)/3;          % [kgm^2], link inertia
tau = 0.0000;           % [Nm], control torque

init_q = pi/4;          % [rad], init joint angle
init_dq = 0.00;         % [rad/s], init angular velocity
q = init_q;             % [rad], current joint angle
dq = init_dq;           % [rad/s], current angular velocity

% Target position parameters
q_d = init_q;           % [rad], target joint angle
dq_d = 0;               % [rad/s], target angular velocity
ddq_d = 0;              % [rad/s^2], target angular acceleration

% Controller gain
Wn = 20;                % [rad/s], natural frequency
Kp = Wn^2;              % [Nm/rad], propotional gain
Kv = 2*Wn;              % [Nm*s/rad], derivative gain

%% Simulation
if (flag_sim == 1)
    % Simulation
    n = 1;
    for time = st:dt:ft
        % Set target trajectory
        cmd = sprintf("loading... %2.2f%%", time/ft*100);
        clc
        disp(cmd);
        if (time < 1)
            q_d = init_q;
            dq_d = 0.0;
            ddq_d = 0.0;
        else
            if (q_d < 90*pi/180)
                q_d = q_d + (45/2*pi/180)*dt;
            else
                q_d = 90*pi/180;
            end
            dq_d = (q_d - sim_q_d(n-1))/dt;
            ddq_d = (dq_d - sim_dq_d(n-1))/dt;
        end
        % Get dynamics
        G = GetGravity(q);
        % Controller
        tau = ddq_d + Kv*(dq_d - dq) + Kp*(q_d - q);
        tq_ctrl = I*tau + G*0.8;
        % Robot model
        % Inverse dynamics
        tau = tq_ctrl;
        [t, y] = ode45('one_link_ex', [0 dt], [q; dq]);
        index = length(y);
        q = y(index, 1);
        dq = y(index, 2);
        
        % Save data
        sim_time(n) = time;     % [sec]
        sim_q(n) = q;           % [rad]
        sim_dq(n) = dq;         % [rad/s]
        sim_q_d(n) = q_d;       % [rad]
        sim_dq_d(n) = dq_d;     % [rad/s]
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
        filename = 'HW3_1_joint_sapce_1_DOF_CTM_PID_controller.gif';
        
        init_x = L*sin(init_q);
        init_y = -L*cos(init_q);
        
        FG1 = figure('Color', [1 1 1]);
        AX = axes('parent', FG1);
        hold on;
        
        p = plot([0 0], [init_x init_y], '-ob', 'Linewidth', linewidth_current);
        
        axis([-1.5 1.5 -1.5 1.5]);
        grid on;
        xlabel('X-axis (m)', 'fontsize', font_size_label);
        ylabel('Y-axis (m)', 'fontsize', font_size_label);
        title('1-DOF Robot', 'fontsize', font_size_title);
        
        n = 1;
        for (time = st:dt:ft)
            cmd = sprintf("Time: %2.2f", time);
            clc
            disp(cmd);  
            q = sim_q(n);
            x = L*sin(q);   y = -L*cos(q);
            Px = [0, x];    Py = [0, y];
            set(p, 'XData', Px, 'YData', Py);
            drawnow
            n = n + 1;
            
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
        % Draw angle
        FG2 = figure('Color', [1 1 1]);
        plot(sim_time, sim_q_d*180/pi, ':k', 'LineWidth', linewidth_target);
        hold on;
        plot(sim_time, sim_q*180/pi, 'r', 'LineWidth', linewidth_current);
        hold on;
        
        axis([st ft 0 120]);
        xticks([st:1:ft]);
        yticks([0:45:90]);
        grid on;
        
        xlabel('time (s)', 'fontsize', font_size_label);
        ylabel('Angle (deg)', 'fontsize', font_size_label);
        title('Joint Space PD CTM Controller', 'fontsize', font_size_title);
        legend('Desired', 'Current');
        
        % Draw angular velocity
        FG3 = figure('Color', [1 1 1]);
        plot(sim_time, sim_dq_d*180/pi, ':k', 'LineWidth', linewidth_target);
        hold on;
        plot(sim_time, sim_dq*180/pi, 'r', 'LineWidth', linewidth_current);
        hold on;
        
        axis([st ft -90 90]);
        xticks([st:1:ft]);
        yticks([-60 0 45/2 60]);
        grid on;
        
        xlabel('time (s)', 'fontsize', font_size_label);
        ylabel('Angle (deg)', 'fontsize', font_size_label);
        title('Joint Space PD CTM Controller', 'fontsize', font_size_title);
        legend('Desired', 'Current');
    end
end
