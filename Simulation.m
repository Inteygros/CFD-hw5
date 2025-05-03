% ��������
M = 200;  % x����������
N = 200;  % y����������
dx = 1.0 / M;
dy = 1.0 / N;
dt = 0.001;      % ʱ�䲽��
nu = 0.001;      % �˶�ճ��
max_iter = 50000; % ����������
tolerance = 1e-6;% ������ֵ
w = 1.7;        % ��ʼ�ɳ�����
w_opt = w;

% ��ʼ��������
psi = zeros(M, N);   % ������
omega = zeros(M, N); % ����
u = zeros(M, N);
v = zeros(M, N);

% ��ʼ�߽�������Woods��ʽ
for i = 2:M-1
    omega(i,N) = -3*(sin(pi*(i-1)*dx))^2/dy+dy*pi^2*cos(2*pi*(i-1)*dx);
end

% ��ʼ�ٶȳ�
for i = 1:M
    u(i, N) = (sin(pi * (i-1)*dx))^2;
end

% ʵʱ��ͼ����
figure;
hold on;
x = linspace(0, 1, M);
y = linspace(0, 1, N);
[X, Y] = meshgrid(x, y);

% ��ѭ��
for t = 1:max_iter
    % ������ʱ�䲽����
    omega_new = omega;
    for i = 2:M-1
        for j = 2:N-1
            % ��ɢ��
            diffusion = (omega(i+1,j) + omega(i-1,j) - 2*omega(i,j))/dx^2 ...
                      + (omega(i,j+1) + omega(i,j-1) - 2*omega(i,j))/dy^2;
            
            % ������
            convection = u(i,j)*(omega(i+1,j) - omega(i-1,j))/(2*dx) ...
                       + v(i,j)*(omega(i,j+1) - omega(i,j-1))/(2*dy);
            
            omega_new(i,j) = omega(i,j) + dt*(nu*diffusion - convection);
        end
    end
    omega = omega_new;
    
    % �״ε�����������ɳ�����
    if t == 1
        disp('��������ɳ�����...');
        w_opt = SOR_w_opt(psi, omega, dx, dy, w, tolerance, M, N);
        fprintf('����ɳ�����: %.2f\n', w_opt);
    end
    
    % SOR�������������
    [psi, ~] = SOR(psi, omega, dx, dy, w_opt, tolerance, M, N);
    
    % �����ٶȳ�
    for i = 2:M-1
        for j = 2:N-1
            u(i,j) = (psi(i,j+1) - psi(i,j-1))/(2*dy);
            v(i,j) = -(psi(i+1,j) - psi(i-1,j))/(2*dx);
        end
    end
    
    % �߽��������£�Woods���й�ʽ��
    for i = 1:M  % ���±߽�
        omega(i,1) = -0.5*omega(i,2)-3*(psi(i,2) - psi(i,1))/dy^2;
        omega(i,N) = -0.5*omega(i,N-1)-3*(psi(i,N-1) - psi(i,N))/dy^2-3*(sin(pi*(i-1)*dx))^2/dy+dy*pi^2*cos(2*pi*(i-1)*dx);
    end
    for j = 1:N  % ���ұ߽�
        omega(1,j) = -0.5*omega(2,j)-3*(psi(2,j) - psi(1,j))/dx^2;
        omega(M,j) = -0.5*omega(M-1,j)-3*(psi(M-1,j) - psi(M,j))/dx^2;
    end
    
    % ���߿��ӻ�
    if mod(t,100) == 0
        cla;
        contourf(X, Y, psi', 50, 'LineColor', 'none');
        streamslice(X, Y, u', v', 3);
        axis equal tight;
        title(sprintf('ʱ�䲽: %d', t));
        drawnow;
    end
end

figure_str = ['figure\m=', num2str(m), '_n=', num2str(n), '.png'];
saveas(gcf, figure_str);

% ����SOR��������ĺ���
function [psi, iter] = SOR(psi, omega, dx, dy, w, tolerance, M, N)
    beta = dx/dy;
    max_error = 1;
    iter = 0;   %��������
    
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

% ����������ɳ����ӵĺ���
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