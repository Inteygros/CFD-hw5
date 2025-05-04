% 参数设置
M = 100;  % x方向网格数
N = 100;  % y方向网格数
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
text_handles = gobjects(0); % 重置句柄

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
    
    %%%%%%%%%%%%%%%%%%% 可视化部分 %%%%%%%%%%%%%%%%%%%
    if mod(t,100) == 0
        cla;
        contourf(X, Y, psi', 50, 'LineColor', 'none');
        streamslice(X, Y, u', v', 3);
        
        % 检测局部极值点
        [maxima, minima] = find_vortex_cores(psi); %极值点
        all_cores = [maxima; minima];
    
        % 筛选符合条件的涡心（绝对值阈值过滤）
        if ~isempty(all_cores)
            [selected_cores, magnitudes] = select_vortex_cores(all_cores, dx, dy);
        
            % 可视化标记
            hold on;
            delete(text_handles); % 清除旧文本
            text_handles = gobjects(0); % 重置句柄
            legend_handles = []; % 涡心强度句柄
            if ~isempty(selected_cores)
                h1=plot(selected_cores(1,1), selected_cores(1,2), 'ro-', 'MarkerSize', 3, 'MarkerFaceColor', 'r')  % 主涡心
                text_handles(end+1) = text(selected_cores(1,1)+0.02, selected_cores(1,2)-0.02,...
                    sprintf('(%.2f, %.2f)\n%.3e',...
                    selected_cores(1,1), selected_cores(1,2), magnitudes(1)),...
                    'Color','r', 'FontSize',10, 'VerticalAlignment','bottom');
                if size(selected_cores,1)>=2
                    h2=plot(selected_cores(2,1), selected_cores(2,2), 'go-', 'MarkerSize', 3, 'MarkerFaceColor', 'g')  % 二次涡心
                    text_handles(end+1) = text(selected_cores(2,1)+0.02, selected_cores(2,2)-0.02,...
                        sprintf('(%.2f, %.2f)\n%.3e',...
                        selected_cores(2,1), selected_cores(2,2), magnitudes(2)),...
                        'Color','g', 'FontSize',10, 'VerticalAlignment','top');
                end
                if size(selected_cores,1)>=3
                    h3=plot(selected_cores(3,1), selected_cores(3,2), 'bo-', 'MarkerSize', 3, 'MarkerFaceColor', 'b') % 三次涡心
                    text_handles(end+1) = text(selected_cores(3,1)-0.02, selected_cores(3,2)+0.02,...
                        sprintf('(%.2f, %.2f)\n%.3e',...
                        selected_cores(3,1), selected_cores(3,2), magnitudes(3)),...
                        'Color','b', 'FontSize',10, 'HorizontalAlignment','right');
                end
            end
            hold off
        end
        
        axis equal tight;
        title(sprintf('时间步: %d', t));
        drawnow;
    end
end


figure_str = ['figure\M=', num2str(M), '_N=', num2str(N), '.png'];
saveas(gcf, figure_str);

%%%%%%%%%%%%%%%%%%% 函数部分 %%%%%%%%%%%%%%%%%%%

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

% 检测局部极值的核心函数
function [maxima, minima] = find_vortex_cores(psi)
    [M,N] = size(psi);
    maxima = []; minima = [];
    
    for i = 3:M-2
        for j = 3:N-2
            current = psi(i,j);
            window = psi(i-1:i+1,j-1:j+1); % 在3*3窗口寻找
            window(2,2) = NaN; 
            
            if current > nanmax(window(:))
                maxima = [maxima; [i,j,current]];  % 记录 [极大值点, 流函数值] 
            elseif current < nanmin(window(:))
                minima = [minima; [i,j,current]]; % 记录 [极小值点, 流函数值] 
            end
        end
    end
end

function [selected, magnitudes] = select_vortex_cores(cores, dx, dy)
    % 从极值点中按绝对值大小筛选得到涡心
    
    min_distance = 0.15;  % 物理空间最小间隔
    
    % 将极值点转换为物理坐标
    cores(:,1) = (cores(:,1)-1)*dx;
    cores(:,2) = (cores(:,2)-1)*dy;
    
    % 按强度绝对值排序
    [~, idx] = sort(abs(cores(:,3)), 'descend');
    sorted = cores(idx,:);
    
    % 空间过滤
    selected = [];
    for k = 1:size(sorted,1)
        valid = true;
        for m = 1:size(selected,1)
            if norm(sorted(k,1:2)-selected(m,1:2)) < min_distance % 如果靠得太近，不认为是涡心
                valid = false;
                break;
            end
        end
        if valid
            selected = [selected; sorted(k,:)];
            if size(selected,1)>=3  % 找到前三个即退出
                break;
            end
        end
    end
    
    % 提取强度信息，select每一行为[x,y,极值]
    magnitudes = [0,0,0];
    if size(selected,1)>=1, magnitudes(1) = selected(1,3); end
    if size(selected,1)>=2, magnitudes(2) = selected(2,3); end
    if size(selected,1)>=3, magnitudes(3) = selected(3,3); end
end