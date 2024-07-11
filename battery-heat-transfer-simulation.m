% --- Battery Parameters ---
rho = 2500;  % Battery density (kg/m^3)
cp = 780;    % Battery specific heat (J/kg.K)
k = 1.5;     % Battery thermal conductivity (W/m.K)
alpha = k / (rho * cp);

% --- Battery Dimensions (m) ---
L = 0.254;
H = 0.1524;
D = 0.1016;

% --- Grid Points (Refined for 3D) ---
nx = 100;
ny = 50;
nz = 25;

% --- Non-Uniform Mesh Generation ---
x = linspace(0, L^0.8, nx).^(1/0.8);
y = linspace(0, H^1.2, ny).^(1/1.2);
z = linspace(0, D^1.2, nz).^(1/1.2);
dx = diff(x);
dy = diff(y);
dz = diff(z);
[X, Y, Z] = meshgrid(x, y, z);

% --- Initial and Boundary Conditions (in Kelvin) ---
T_initial = 298.15;
T_inlet = 293.15;
T_outlet = 303.15;
T_amb = 298.15;

% --- Convective Heat Transfer Coefficients (W/m^2.K) ---
h_coolant = 1500;
h_air = 10;

% --- Heat Generation Function (W/m^3) - Simplified Cycle ---
q_dot = @(t) 10000 + 15000 * sin(2 * pi * t / 600).^2;

% --- Time Stepping and Stability ---
dt = 0.005;
max_iter = 50000;

% --- Convergence Criterion ---
convergence_limit = 1e-4;

% --- Storage ---
T = T_initial * ones(ny, nx, nz);
cpu_times = [];

% --- Progress Tracking Variables ---
print_frequency = 100;
last_print_time = tic;

% --- Finite Difference Solution ---
start_time = tic;
for iter = 1:max_iter
    T_old = T;
    
    % Apply Boundary Conditions
    T(:, :, 1) = T_inlet + (T_outlet - T_inlet) * Y(:, :, 1) / H;
    
    % Top, Sides, Front, & Back (Convection to Ambient)
    T(:, :, end) = (h_air * dt * T_amb + rho * cp * dz(end) * T_old(:, :, end)) / (rho * cp * dz(end) + h_air * dt);
    T(1, :, :)  = (h_air * dt * T_amb + rho * cp * dy(1)  * T_old(1, :, :))  / (rho * cp * dy(1)  + h_air * dt);
    T(:, 1, :)  = (h_air * dt * T_amb + rho * cp * dx(1)  * T_old(:, 1, :))  / (rho * cp * dx(1)  + h_air * dt);
    T(:, end, :) = (h_air * dt * T_amb + rho * cp * dx(end) * T_old(:, end, :)) / (rho * cp * dx(end) + h_air * dt);
    
    % Finite Difference Calculation (Interior Nodes)
    for i = 2:nx-1
        for j = 2:ny-1
            for k = 2:nz-1
                T(j, i, k) = T(j, i, k) + alpha * dt * ( ...
                    (T_old(j, i+1, k) - 2*T_old(j, i, k) + T_old(j, i-1, k)) / dx(i)^2 + ...
                    (T_old(j+1, i, k) - 2*T_old(j, i, k) + T_old(j-1, i, k)) / dy(j)^2 + ...
                    (T_old(j, i, k+1) - 2*T_old(j, i, k) + T_old(j, i, k-1)) / dz(k)^2 ...
                ) + q_dot(iter * dt) * dt / (rho * cp);
            end
        end
    end
    
    % Convergence Check & Progress Update
    if max(abs(T(:) - T_old(:))) < convergence_limit
        fprintf('Converged after %d iterations.\n', iter);
        break;
    end
    
    % Enhanced Progress Tracking
    current_time = toc(start_time);
    if current_time - last_print_time >= 1.0 || mod(iter, print_frequency) == 0
        cpu_times = [cpu_times, current_time];
        fprintf('Iteration: %d, Elapsed time: %.2f seconds, Max Temp: %.2f K\n', iter, cpu_times(end), max(T(:)));
        last_print_time = tic;
    end
end

% 3D Plotting (Example at mid-depth)
figure;
mid_depth_idx = floor(nz / 2);
surf(X(:, :, mid_depth_idx), Y(:, :, mid_depth_idx), T(:, :, mid_depth_idx) - 273.15, 'EdgeColor', 'none');
colorbar;
xlabel('Length (m)');
ylabel('Height (m)');
zlabel('Temperature (°C)');
title('Temperature Distribution at Mid-Depth');
