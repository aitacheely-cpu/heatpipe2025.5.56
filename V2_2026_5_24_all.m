% =========================================================================
% TPTL 稳态分布参数模型 (lsqnonlin 求解器 + MultiStart 全局搜索)
% 末尾已添加自动调用 TPTL_postprocess 后处理函数逻辑
% =========================================================================
clear; clc; close all;
%% 1. 几何参数与边界条件
geo.D_o = 0.025; % [m] 外径
geo.D_i = 0.020; % [m] 内径
geo.A_c = pi/4 * geo.D_i^2; % [m^2] 流通截面积
% 管段离散参数
geo.L_evap = 0.5;   geo.N_evap = 50;  geo.dz_evap = geo.L_evap / geo.N_evap;
geo.L_riser_v = 2.0; geo.N_riser_v = 10; geo.dz_riser_v = geo.L_riser_v / geo.N_riser_v; 
geo.L_riser_h = 1.0; geo.N_riser_h = 10; geo.dz_riser_h = geo.L_riser_h / geo.N_riser_h; 
geo.L_cond = 0.5;   geo.N_cond = 50;  geo.dz_cond = geo.L_cond / geo.N_cond;
geo.L_down_v = 2.0;  geo.N_down_v = 10;  geo.dz_down_v = geo.L_down_v / geo.N_down_v;    
geo.L_down_h = 1.0;  geo.N_down_h = 10;  geo.dz_down_h = geo.L_down_h / geo.N_down_h;    
% --- 材质与管件参数 ---
geo.k_steel = 15.0;  % [W/mK] 不锈钢导热系数
geo.K_bend  = 1.2;   % 90度弯头单相局部阻力系数
% 壁温边界条件 (蒸发段 250°C，冷凝段 20°C)
geo.T_wall_evap = 273.15 + 250; % [K]
geo.T_wall_cond = 273.15 + (20); % [K]
M_charge = 0.35;                % [kg] 系统初始充注量

%% 2. 求解器初始猜测值与边界设置
X0 = [2.3, 0.8, 200];        % [P0(MPa), H0(m), G(kg/m²s)]
lb = [0.000611657, 0.02, 5];
% 动态获取蒸发段壁温对应的饱和压力作为压力的绝对上限
P_max = safe_IAPWS('psat_T', geo.T_wall_evap); 
ub = [P_max, geo.L_down_v, 2000];

%% 3. 配置并调用 MultiStart 全局求解器
options = optimoptions('lsqnonlin', 'Display', 'off', 'Algorithm', 'trust-region-reflective', ...
    'FunctionTolerance', 1e-4, 'StepTolerance', 1e-4, 'MaxIterations', 200);
problem = createOptimProblem('lsqnonlin', 'x0', X0, ...
    'objective', @(X) TPTL_residuals(X, geo, M_charge), 'lb', lb, 'ub', ub, 'options', options);
ms = MultiStart('Display', 'iter', 'UseParallel', false);
fprintf('开始调用 MultiStart 全局搜索...\n');
tic;
[X_sol, resnorm, exitflag, output, solutions] = run(ms, problem, 10); 
toc;

fprintf('\n=== 全局求解完成 ===\n');
fprintf('最优解: P0 = %.4f MPa, H0 = %.4f m, G = %.2f kg/m²s\n', X_sol(1), X_sol(2), X_sol(3));
fprintf('最小残差平方和 = %e\n', resnorm);

%% 4. 遍历并输出所有误差较小的局部最优解分布
fprintf('\n=== 收敛的局部最优解分布 (按误差从小到大排序) ===\n');
error_threshold = 1e-3; 
valid_solution_count = 0;
for i = 1:length(solutions)
    x_i = solutions(i).X;
    err_i = solutions(i).Fval;
    
    if err_i < error_threshold
        valid_solution_count = valid_solution_count + 1;
        fprintf('解 %2d: P0 = %.4f MPa, H0 = %.4f m, G = %6.2f kg/m²s | 残差 = %e\n', ...
            i, x_i(1), x_i(2), x_i(3), err_i);
    end
end
if valid_solution_count == 0
    fprintf('除全局最优解外，没有其他残差小于 %e 的解。\n', error_threshold);
else
    fprintf('共找到 %d 个低误差解。\n', valid_solution_count);
end

%% 5. 【新增核心修改】：提取最优解并全自动触发后处理函数
P0_opt = X_sol(1);
H0_opt = X_sol(2);
G_opt  = X_sol(3);

fprintf('\n==================================================\n');
fprintf('>>> 正在全自动调用后处理函数 TPTL_postprocess ...\n');
fprintf('>>> 传入最优参数: P0 = %.4f MPa, H0 = %.4f m, G = %.2f kg/(m²·s)\n', P0_opt, H0_opt, G_opt);
fprintf('==================================================\n');

% 执行后处理函数（自动生成 2x4 分布曲线图，并导出规范命名的 Excel 工作表）
TPTL_postprocess(P0_opt, H0_opt, G_opt);


%% ===================== 局部函数区 =====================
function F = TPTL_residuals(X, geo, M_charge)
    P0 = X(1); H0 = X(2); G = X(3);
    try
        h0 = safe_IAPWS('hL_p', P0); P_curr = P0; h_curr = h0; Total_Mass = 0;
        R_wall = (geo.D_i / (2 * geo.k_steel)) * log(geo.D_o / geo.D_i);
        
        % (1) 下降段 - 竖直
        L_unfilled = max(0, geo.L_down_v - H0); curr_z = 0;
        for i = 1:geo.N_down_v
            z_start = curr_z; curr_z = curr_z + geo.dz_down_v; z_end = curr_z;
            if L_unfilled >= z_end
                rho = 1 / safe_IAPWS('vV_p', P_curr); 
                Total_Mass = Total_Mass + rho * geo.A_c * geo.dz_down_v;
            elseif L_unfilled <= z_start
                [rho, mu, ~, ~, ~, ~, ~] = get_fluid_state(P_curr, h_curr);
                Re = G * geo.D_i / mu; f = calc_fanning_friction(Re);      
                dP_f = f * (geo.dz_down_v/geo.D_i) * (G^2)/(2*rho); dP_g = rho * 9.81 * geo.dz_down_v; 
                P_curr = P_curr + ((dP_g - dP_f) / 1e6); Total_Mass = Total_Mass + rho * geo.A_c * geo.dz_down_v;
            else
                L_gas = L_unfilled - z_start; L_liq = z_end - L_unfilled;
                rho_g = 1 / safe_IAPWS('vV_p', P_curr); mass_g = rho_g * geo.A_c * L_gas;
                [rho_l, mu_l, ~, ~, ~, ~, ~] = get_fluid_state(P_curr, h_curr);
                Re = G * geo.D_i / mu_l; f_l = calc_fanning_friction(Re);
                dP_f_l = f_l * (L_liq/geo.D_i) * (G^2)/(2*rho_l); dP_g_l = rho_l * 9.81 * L_liq;
                Total_Mass = Total_Mass + mass_g + rho_l * geo.A_c * L_liq; P_curr = P_curr + ((dP_g_l - dP_f_l) / 1e6);
            end
        end
        
        % 弯头1
        [~, ~, ~, x_th, ~, ~, ~] = get_fluid_state(P_curr, h_curr);
        P_curr = P_curr - calc_local_dp(P_curr, x_th, G, geo.K_bend) / 1e6;
        
        % (2) 下降段 - 水平
        for i = 1:geo.N_down_h
            [rho, mu, ~, x_th, ~, rho_mix, ~] = get_fluid_state(P_curr, h_curr);
            if (x_th > 0 && x_th < 1)
                rho_l = 1 / safe_IAPWS('vL_p', P_curr); mu_l = safe_IAPWS('mu_ph', P_curr, safe_IAPWS('hL_p', P_curr));
                rho_g = 1 / safe_IAPWS('vV_p', P_curr); mu_g = safe_IAPWS('mu_ph', P_curr, safe_IAPWS('hV_p', P_curr));
                dP_f = calc_two_phase_friction_dp(x_th, G, geo.D_i, rho_l, mu_l, rho_g, mu_g, geo.dz_down_h);
            else
                Re = G * geo.D_i / mu; f = calc_fanning_friction(Re); dP_f = f * (geo.dz_down_h/geo.D_i) * (G^2)/(2*rho);
            end
            P_curr = P_curr - (dP_f / 1e6); Total_Mass = Total_Mass + rho_mix * geo.A_c * geo.dz_down_h;
        end
        
        % 弯头2
        [~, ~, ~, x_th, ~, ~, ~] = get_fluid_state(P_curr, h_curr);
        P_curr = P_curr - calc_local_dp(P_curr, x_th, G, geo.K_bend) / 1e6;
        
        % (3) 蒸发段 - 竖直
        for i = 1:geo.N_evap
            [~, ~, T_fluid, x_th, ~, rho_mix, mu_mix] = get_fluid_state(P_curr, h_curr);
            if x_th <= 0 || x_th >= 1
                Re = G * geo.D_i / mu_mix;
                Cp_fluid = safe_IAPWS('cp_ph', P_curr, h_curr); k_fluid = safe_IAPWS('k_ph', P_curr, h_curr);
                Pr = (Cp_fluid * 1000) * mu_mix / k_fluid;
                h_wall = safe_IAPWS('h_pT', P_curr, geo.T_wall_evap); mu_wall = safe_IAPWS('mu_ph', P_curr, h_wall); 
                z_curr = max(i * geo.dz_evap, 1e-4);
                h_htc = calc_single_phase_htc(Re, Pr, k_fluid, geo.D_i, z_curr, mu_mix, mu_wall, true);
                U_i = 1 / (1 / h_htc + R_wall); q_node = U_i * (geo.T_wall_evap - T_fluid);
            else
                q_guess = 5000; tol_q = 1e-3;
                for iter = 1:500
                    T_wall_in = geo.T_wall_evap - q_guess * R_wall;
                    h_htc = calc_boiling_htc(P_curr, q_guess, x_th); 
                    q_new = h_htc * (T_wall_in - T_fluid); 
                    if abs(q_new - q_guess) / max(q_guess, 1) < tol_q, break; end
                    q_guess = 0.8 * q_guess + 0.2 * q_new; 
                end
                q_node = q_new;
            end
            dQ = q_node * pi * geo.D_i * geo.dz_evap; dh = (dQ / (G * geo.A_c)) / 1000; 
            dP_g = rho_mix * 9.81 * geo.dz_evap;
            if (x_th > 0 && x_th < 1)
                rho_l = 1 / safe_IAPWS('vL_p', P_curr); mu_l = safe_IAPWS('mu_ph', P_curr, safe_IAPWS('hL_p', P_curr));
                rho_g = 1 / safe_IAPWS('vV_p', P_curr); mu_g = safe_IAPWS('mu_ph', P_curr, safe_IAPWS('hV_p', P_curr));
                dP_f = calc_two_phase_friction_dp(x_th, G, geo.D_i, rho_l, mu_l, rho_g, mu_g, geo.dz_evap);
            else
                Re = G * geo.D_i / mu_mix; f = calc_fanning_friction(Re); dP_f = f * (geo.dz_evap/geo.D_i) * (G^2)/(2*rho_mix);
            end
            P_curr = P_curr - ((dP_g + dP_f) / 1e6); h_curr = h_curr + dh; Total_Mass = Total_Mass + rho_mix * geo.A_c * geo.dz_evap;
        end

        % (4) 上升段 - 竖直
        for i = 1:geo.N_riser_v
            [~, ~, ~, x_th, ~, rho_mix, mu_mix] = get_fluid_state(P_curr, h_curr);
            dP_g = rho_mix * 9.81 * geo.dz_riser_v;
            if (x_th > 0 && x_th < 1)
                rho_l = 1 / safe_IAPWS('vL_p', P_curr); mu_l = safe_IAPWS('mu_ph', P_curr, safe_IAPWS('hL_p', P_curr));
                rho_g = 1 / safe_IAPWS('vV_p', P_curr); mu_g = safe_IAPWS('mu_ph', P_curr, safe_IAPWS('hV_p', P_curr));
                dP_f = calc_two_phase_friction_dp(x_th, G, geo.D_i, rho_l, mu_l, rho_g, mu_g, geo.dz_riser_v);
            else
                Re = G * geo.D_i / mu_mix; f = calc_fanning_friction(Re); dP_f = f * (geo.dz_riser_v/geo.D_i) * (G^2)/(2*rho_mix);
            end
            P_curr = P_curr - ((dP_g + dP_f) / 1e6); Total_Mass = Total_Mass + rho_mix * geo.A_c * geo.dz_riser_v;
        end
        
        % 弯头3
        [~, ~, ~, x_th, ~, ~, ~] = get_fluid_state(P_curr, h_curr);
        P_curr = P_curr - calc_local_dp(P_curr, x_th, G, geo.K_bend) / 1e6;
        
        % (5) 上升段 - 水平
        for i = 1:geo.N_riser_h
            [~, ~, ~, x_th, ~, rho_mix, mu_mix] = get_fluid_state(P_curr, h_curr);
            if (x_th > 0 && x_th < 1)
                rho_l = 1 / safe_IAPWS('vL_p', P_curr); mu_l = safe_IAPWS('mu_ph', P_curr, safe_IAPWS('hL_p', P_curr));
                rho_g = 1 / safe_IAPWS('vV_p', P_curr); mu_g = safe_IAPWS('mu_ph', P_curr, safe_IAPWS('hV_p', P_curr));
                dP_f = calc_two_phase_friction_dp(x_th, G, geo.D_i, rho_l, mu_l, rho_g, mu_g, geo.dz_riser_h);
            else
                Re = G * geo.D_i / mu_mix; f = calc_fanning_friction(Re); dP_f = f * (geo.dz_riser_h/geo.D_i) * (G^2)/(2*rho_mix);
            end
            P_curr = P_curr - (dP_f / 1e6); Total_Mass = Total_Mass + rho_mix * geo.A_c * geo.dz_riser_h;
        end
        
        % 弯头4
        [~, ~, ~, x_th, ~, ~, ~] = get_fluid_state(P_curr, h_curr);
        P_curr = P_curr - calc_local_dp(P_curr, x_th, G, geo.K_bend) / 1e6;
        
        % (6) 冷凝段 - 竖直
        for i = 1:geo.N_cond
            [~, ~, T_fluid, x_th, ~, rho_mix, mu_mix] = get_fluid_state(P_curr, h_curr);
            if x_th > 0 && x_th < 1
                k_l = safe_IAPWS('k_ph', P_curr, safe_IAPWS('hL_p', P_curr)); rho_l = 1 / safe_IAPWS('vL_p', P_curr);
                mu_l = safe_IAPWS('mu_ph', P_curr, safe_IAPWS('hL_p', P_curr)); h_fg = (safe_IAPWS('hV_p', P_curr) - safe_IAPWS('hL_p', P_curr)) * 1000; 
                z_curr = max(i * geo.dz_cond, 1e-4); 
                q_guess = 5000; tol_q = 1e-3;
                for iter = 1:500
                    T_wall_in = geo.T_wall_cond + q_guess * R_wall; dT_film = max(T_fluid - T_wall_in, 0.1); 
                    h_htc = 1.13 * ( (k_l^3 * rho_l^2 * h_fg * 9.81) / (mu_l * z_curr * dT_film) )^0.25; q_new = h_htc * dT_film;
                    if abs(q_new - q_guess) / max(q_guess, 1) < tol_q, break; end
                    q_guess = 0.8 * q_guess + 0.2 * q_new; 
                end
                q_node = q_new;
            else
                Re = G * geo.D_i / mu_mix; Cp_fluid = safe_IAPWS('cp_ph', P_curr, h_curr); k_fluid = safe_IAPWS('k_ph', P_curr, h_curr);
                Pr = (Cp_fluid * 1000) * mu_mix / k_fluid; T_wall_in_approx = max(T_fluid - 10, 273.16); 
                h_wall = safe_IAPWS('h_pT', P_curr, T_wall_in_approx); mu_wall = safe_IAPWS('mu_ph', P_curr, h_wall); 
                z_curr = max(i * geo.dz_cond, 1e-4); h_htc = calc_single_phase_htc(Re, Pr, k_fluid, geo.D_i, z_curr, mu_mix, mu_wall, false);
                U_i = 1 / (1 / h_htc + R_wall); q_node = U_i * max(T_fluid - geo.T_wall_cond, 0.1);
            end
            dQ = q_node * pi * geo.D_i * geo.dz_cond; dh = -(dQ / (G * geo.A_c)) / 1000;
            dP_g = rho_mix * 9.81 * geo.dz_cond;
            if (x_th > 0 && x_th < 1)
                rho_l = 1 / safe_IAPWS('vL_p', P_curr); mu_l = safe_IAPWS('mu_ph', P_curr, safe_IAPWS('hL_p', P_curr));
                rho_g = 1 / safe_IAPWS('vV_p', P_curr); mu_g = safe_IAPWS('mu_ph', P_curr, safe_IAPWS('hV_p', P_curr));
                dP_f = calc_two_phase_friction_dp(x_th, G, geo.D_i, rho_l, mu_l, rho_g, mu_g, geo.dz_cond);
            else
                Re = G * geo.D_i / mu_mix; f = calc_fanning_friction(Re); dP_f = f * (geo.dz_cond/geo.D_i) * (G^2)/(2*rho_mix);
            end
            P_curr = P_curr + (dP_g / 1e6) - (dP_f / 1e6); h_curr = h_curr + dh; Total_Mass = Total_Mass + rho_mix * geo.A_c * geo.dz_cond;
        end
        F = [(P_curr - P0) / P0, (h_curr - h0) / h0, (Total_Mass - M_charge) / M_charge];
    catch ME
        fprintf(2, '\n【发现致命报错，求解器强制停止】\n报错内容: %s\n在函数 %s 的第 %d 行\n', ME.message, ME.stack(1).name, ME.stack(1).line);
        F = [1e6, 1e6, 1e6]; 
    end
end

function dP_local = calc_local_dp(P, x, G, K_bend)
    rho_l = 1 / safe_IAPWS('vL_p', P); rho_g = 1 / safe_IAPWS('vV_p', P); dP_l = K_bend * (G^2) / (2 * rho_l);
    if x <= 0, dP_local = dP_l; elseif x >= 1, dP_local = K_bend * (G^2) / (2 * rho_g);
    else
        x_safe = max(min(x, 0.999), 1e-4); Xb = (rho_g/rho_l)^0.5 * (1-x_safe)/x_safe; 
        C = (1 + (3.2-1)*(rho_g)^0.5) * ((rho_g/rho_l)^0.5 + (rho_g/rho_l)^-0.5); phi_c_sq = 1 + C/Xb + 1/Xb^2; dP_local = dP_l * phi_c_sq;
    end
end

function f_fan = calc_fanning_friction(Re)
    Re = max(Re, 1e-4); if Re < 2000, f_fan = 64 / Re; else, f_fan = 0.184 / Re^0.2; end
end

% -------------------------------------------------------------------------
% Lockhart-Martinelli 两相摩擦压降 (Chisholm 表达式 - 严格修正版)
% -------------------------------------------------------------------------
function dP_f = calc_two_phase_friction_dp(x, G, D, rho_l, mu_l, rho_g, mu_g, dz)
    % 1. 实际液相、气相分流质量流速
    G_l = G * (1 - x);
    G_g = G * x;
    
    % 2. 实际分相雷诺数（用于确定各自流态）
    Re_l = max(G_l * D / mu_l, 1e-4);
    Re_g = max(G_g * D / mu_g, 1e-4);
    
    % 3. 计算分相单相摩擦系数
    f_l = calc_fanning_friction(Re_l);
    f_g = calc_fanning_friction(Re_g);
    
    % 4. Martinelli 参数 X
    X = ((1 - x) / x) * sqrt((f_l / f_g) * (rho_g / rho_l));
    X = max(X, 1e-6);  % 避免除零
    
    % 5. 根据实际流态确定 C 值 (与图2表格一致)
    if Re_l < 2000 && Re_g < 2000
        C = 5;
    elseif Re_l < 2000 && Re_g >= 2000
        C = 12;
    elseif Re_l >= 2000 && Re_g < 2000
        C = 10;
    else
        C = 20;
    end
    
    % 6. 【修正点 1】：计算液相单独流动（Liquid-only）的摩擦压降梯度 dpdz_l
    % 公式形式必须与主程序单相压降（f * L/D * G^2 / (2*rho)）保持绝对一致
    dpdz_l = f_l * (G_l^2 / rho_l) / (2 * D);
    
    % 7. 两相摩擦乘子 (Chisholm 表达式)
    phi_l_sq = 1 + C / X + 1 / X^2;
    
    % 8. 【修正点 2】：两相阻力由液相单独流动梯度 dpdz_l 与 phi_l_sq 耦合得到
    dP_f = phi_l_sq * dpdz_l * dz;
end

function [rho_single, mu_single, T_fluid, x_th, alpha, rho_mix, mu_mix] = get_fluid_state(P, h)
    hL = safe_IAPWS('hL_p', P); hV = safe_IAPWS('hV_p', P); Tsat = safe_IAPWS('Tsat_p', P); x_th = (h - hL) / (hV - hL);
    if x_th <= 0 
        x_th = 0; T_fluid = safe_IAPWS('T_ph', P, h); rho_single = 1 / safe_IAPWS('v_ph', P, h); mu_single = safe_IAPWS('mu_ph', P, h); alpha = 0; rho_mix = rho_single; mu_mix = mu_single;
    elseif x_th >= 1
        x_th = 1; T_fluid = safe_IAPWS('T_ph', P, h); rho_single = 1 / safe_IAPWS('v_ph', P, h); mu_single = safe_IAPWS('mu_ph', P, h); alpha = 1; rho_mix = rho_single; mu_mix = mu_single;
    else
        T_fluid = Tsat; vL = safe_IAPWS('vL_p', P); vV = safe_IAPWS('vV_p', P); rho_l = 1 / vL; rho_g = 1 / vV;
        mu_l = safe_IAPWS('mu_ph', P, hL + 0.1); mu_g = safe_IAPWS('mu_ph', P, hV - 0.1);
        alpha = max(0, min(1, 1 / (1 + ((1-x_th)/x_th) * (rho_g/rho_l)^0.89 * (mu_l/mu_g)^0.18)));
        rho_mix = alpha * rho_g + (1 - alpha) * rho_l; mu_mix = 1 / ( (x_th/mu_g) + ((1-x_th)/mu_l) ); rho_single = rho_l; mu_single = mu_l;
    end
end

function val = safe_IAPWS(prop, varargin)
    val = IAPWS_IF97(prop, varargin{:}); if any(isnan(val)), error('TPTL:NaN_Error', 'IAPWS_IF97 out of bounds.'); end
end

function h_htc = calc_boiling_htc(P_MPa, q, x)
    M = 18.015; P_crit = 22.064; Pr = max(P_MPa / P_crit, 1e-4); q = max(q, 10); x = max(0, min(x, 0.99));
    h_htc = 540 * Pr^0.394 *(1-x)^-0.65* M^(-0.5) * q^0.54;
end

function h_htc = calc_single_phase_htc(Re, Pr, k, D, z, mu_fluid, mu_wall, is_heating)
    Re = max(Re, 10); Pr = max(Pr, 0.1); z = max(z, 1e-4);
    if Re < 2200
        h_htc = 1.24 * (k * (1 + 5 * exp(-z / (10 * D))) / D) * (Re * Pr * (D / z))^(1/3) * (mu_fluid / mu_wall);
    else
        n = 0.3; if is_heating, n = 0.4; end
        h_htc = 0.027 * (k * (1 + 7 * (D / z)) / D) * Re^0.8 * Pr^n * (mu_fluid / mu_wall)^0.14;
    end
end