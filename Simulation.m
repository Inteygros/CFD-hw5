% 参数设置
M = 200;  % x方向网格数
N = 200;  % y方向网格数
dx = 1.0 / M;
dy = 1.0 / N;
dt = 0.001;      % 时间步长
nu = 0.001;      % 运动粘度
max_iter = 50000; % 最大迭代次数
tolerance = 1e-6;% 收敛阈值
w = 1.7;        % 初始松弛因子
w_opt = w;

% 初始化场变量
psi = zeros(M, N);   % 流函数
omega = zeros(M, N); % 涡量
u = zeros(M, N);
v = zeros(M, N);

% 初始边界条件用Woods公式
for i = 2:M-1
    omega(i,N) = -3*(sin(pi*(i-1)*dx))^2/dy+dy*pi^2*cos(2*pi*(i-1)*dx);
end

% 初始速度场
for i = 1:M
    u(i, N) = (sin(pi * (i-1)*dx))^2;
end

% 实时绘图设置
figure;
hold on;
x = linspace(0, 1, M);
y = linspace(0, 1, N);
[X, Y] = meshgrid(x, y);

% 主循环
for t = 1:max_iter
    % 计算新时间步涡量
    omega_new = omega;
    for i = 2:M-1
        for j = 2:N-1
            % 扩散项
            diffusion = (omega(i+1,j) + omega(i-1,j) - 2*omega(i,j))/dx^2 ...
                      + (omega(i,j+1) + omega(i,j-1) - 2*omega(i,j))/dy^2;
            
            % 对流项
            convection = u(i,j)*(omega(i+1,j) - omega(i-1,j))/(2*dx) ...
                       + v(i,j)*(omega(i,j+1) - omega(i,j-1))/(2*dy);
            
            omega_new(i,j) = omega(i,j) + dt*(nu*diffusion - convection);
        end
    end
    omega = omega_new;
    
    % 首次迭代计算最佳松弛因子
    if t == 1
        disp('计算最佳松弛因子...');
        w_opt = SOR_w_opt(psi, omega, dx, dy, w, tolerance, M, N);
        fprintf('最佳松弛因子: %.2f\n', w_opt);
    end
    
    % SOR迭代求解流函数
    [psi, ~] = SOR(psi, omega, dx, dy, w_opt, tolerance, M, N);
    
    % 更新速度场
    for i = 2:M-1
        for j = 2:N-1
            u(i,j) = (psi(i,j+1) - psi(i,j-1))/(2*dy);
            v(i,j) = -(psi(i+1,j) - psi(i-1,j))/(2*dx);
        end
    end
    
    % 边界条件更新（Woods壁涡公式）
    for i = 1:M  % 上下边界
        omega(i,1) = -0.5*omega(i,2)-3*(psi(i,2) - psi(i,1))/dy^2;
        omega(i,N) = -0.5*omega(i,N-1)-3*(psi(i,N-1) - psi(i,N))/dy^2-3*(sin(pi*(i-1)*dx))^2/dy+dy*pi^2*cos(2*pi*(i-1)*dx);
    end
    for j = 1:N  % 左右边界
        omega(1,j) = -0.5*omega(2,j)-3*(psi(2,j) - psi(1,j))/dx^2;
        omega(M,j) = -0.5*omega(M-1,j)-3*(psi(M-1,j) - psi(M,j))/dx^2;
    end
    
    % 流线可视化
    if mod(t,100) == 0
        cla;
        contourf(X, Y, psi', 50, 'LineColor', 'none');
        streamslice(X, Y, u', v', 3);
        axis equal tight;
        title(sprintf('时间步: %d', t));
        drawnow;
    end
end

figure_str = ['figure\m=', num2str(m), '_n=', num2str(n), '.png'];
saveas(gcf, figure_str);

% 定义SOR方法计算的函数
function [psi, iter] = SOR(psi, omega, dx, dy, w, tolerance, M, N)
    beta = dx/dy;
    max_error = 1;
    iter = 0;   %迭代次数
    
    while max_error > tolerance
        iter = iter + 1;
        max_error = 0;
        
        for i = 2:M-1
            for j = 2:N-1
                temp = psi(i,j);
                psi(i,j) = w/(2*(1+beta^2)) * (psi(i+1,j) + psi(i-1,j) + ...
                          beta^2*(psi(i,j+1) + psi(i,j-1)) + dx^2*omega(i,j)) ...
                          + (1-w)*psi(i,j);
                current_error = abs(psi(i,j) - temp);
                if current_error > max_error
                    max_error = current_error;
                end
            end
        end
    end
end

% 定义找最佳松弛因子的函数
function w_opt = SOR_w_opt(p, omega, dx, dy, w, tolerance, M, N)
    [~, min_iter] = SOR(p, omega, dx, dy, w, tolerance, M, N);
    w_opt = w;
    
    for w1 = (w+0.01):0.01:1.99
        [~, iter] = SOR(p, omega, dx, dy, w1, tolerance, M, N);
        if iter < min_iter
            min_iter = iter;
            w_opt = w1;
        end
    end
end