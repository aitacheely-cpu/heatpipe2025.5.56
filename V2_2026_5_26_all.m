% =========================================================================
% TPTL 稳态分布参数模型 (异构变径并联环路专用版)
% 结构：3m垂直下降(单0.05) -> 3m水平下降(单0.05) -> 2.83m竖直蒸发(3并联0.02) 
%      -> 1.2m垂直上升(单0.05) -> 3m水平冷凝(3并联0.02)
% =========================================================================
clear; clc; close all;
%% 1. 几何参数与边界条件
% (1) 垂直下降段 (单根)
geo.L_down_v = 3.0;  geo.N_down_v = 30;  geo.dz_down_v = geo.L_down_v / geo.N_down_v;
geo.D_i_down = 0.05; geo.A_c_down = pi/4 * geo.D_i_down^2;

% (2) 水平下降段 (单根)
geo.L_down_h = 3.0;  geo.N_down_h = 30;  geo.dz_down_h = geo.L_down_h / geo.N_down_h;

% (3) 竖直蒸发段 (3根并联)
geo.N_tubes_evap = 3;
geo.L_evap = 2.83;   geo.N_evap = 50;    geo.dz_evap = geo.L_evap / geo.N_evap;
geo.D_i_evap = 0.02; geo.D_o_evap = 0.025; geo.A_c_evap = pi/4 * geo.D_i_evap^2;

% (4) 垂直上升段 (单根)
geo.L_riser_v = 1.2; geo.N_riser_v = 20; geo.dz_riser_v = geo.L_riser_v / geo.N_riser_v;
geo.D_i_riser = 0.05; geo.A_c_riser = pi/4 * geo.D_i_riser^2;

% (5) 水平冷凝段 (3根并联)
geo.N_tubes_cond = 3;
geo.L_cond = 3.0;    geo.N_cond = 50;    geo.dz_cond = geo.L_cond / geo.N_cond;
geo.D_i_cond = 0.02; geo.D_o_cond = 0.025; geo.A_c_cond = pi/4 * geo.D_i_cond^2;

% 材质与管件参数
geo.k_steel = 22.0;  
geo.K_bend  = 1.2;   

% 壁温边界条件
geo.T_wall_evap = 273.15 + 250; 
geo.T_wall_cond = 273.15 + (20); 
M_charge = 8; % [kg] 系统初始充注量

%% 2. 求解器初始猜测值与边界设置
% 注意：迭代变量X(3)变更为 总质量流量 W (kg/s)
X0 = [2.3, 1.5, 0.1];        % [P0(MPa), H0(m), W(kg/s)]
lb = [0.000611657, 0.02, 1e-4];
P_max = safe_IAPWS('psat_T', geo.T_wall_evap); 
ub = [P_max, geo.L_down_v, 5];

%% 3. 配置并调用 MultiStart 全局求解器
options = optimoptions('lsqnonlin', 'Display', 'off', 'Algorithm', 'trust-region-reflective', ...
    'FunctionTolerance', 1e-4, 'StepTolerance', 1e-4, 'MaxIterations', 200);
problem = createOptimProblem('lsqnonlin', 'x0', X0, ...
    'objective', @(X) TPTL_residuals(X, geo, M_charge), 'lb', lb, 'ub', ub, 'options', options);
ms = MultiStart('Display', 'iter', 'UseParallel', false);
fprintf('开始调用 MultiStart 全局搜索...\n');
tic;
[X_sol, resnorm, ~, ~, solutions] = run(ms, problem, 10); 
toc;

fprintf('\n=== 全局求解完成 ===\n');
fprintf('最优解: P0 = %.4f MPa, H0 = %.4f m, 总流量 W = %.4f kg/s\n', X_sol(1), X_sol(2), X_sol(3));
fprintf('最小残差平方和 = %e\n', resnorm);

%% 4. 后处理调用
fprintf('\n==================================================\n');
fprintf('>>> 正在全自动调用后处理函数 TPTL_postprocess ...\n');
TPTL_postprocess_2026_5_26_all(X_sol(1), X_sol(2), X_sol(3));

%% ===================== 局部函数区 =====================
function F = TPTL_residuals(X, geo, M_charge)
    P0 = X(1); H0 = X(2); W = X(3); % W 为总质量流量 kg/s
    try
        h0 = safe_IAPWS('hL_p', P0);
        P_curr = P0; h_curr = h0; 
        Total_Mass = 0;
        R_wall_evap = (geo.D_i_evap / (2 * geo.k_steel)) * log(geo.D_o_evap / geo.D_i_evap);
        R_wall_cond = (geo.D_i_cond / (2 * geo.k_steel)) * log(geo.D_o_cond / geo.D_i_cond);
        % 计算各段的局部质量流速 G (kg/m²s)
        G_down  = W / geo.A_c_down;
        G_evap  = (W / geo.N_tubes_evap) / geo.A_c_evap;
        G_riser = W / geo.A_c_riser;
        G_cond  = (W / geo.N_tubes_cond) / geo.A_c_cond;
        % (1) 下降段 - 竖直 (单根)
        L_unfilled = max(0, geo.L_down_v - H0); curr_z = 0;
        for i = 1:geo.N_down_v
            z_start = curr_z; 
            curr_z = curr_z + geo.dz_down_v; 
            z_end = curr_z;
            if L_unfilled >= z_end
                rho = 1 / safe_IAPWS('vV_p', P_curr); 
                Total_Mass = Total_Mass + rho * geo.A_c_down * geo.dz_down_v;
            elseif L_unfilled <= z_start
                [rho, mu, ~, ~, ~, ~, ~] = get_fluid_state(P_curr, h_curr);
                Re = G_down * geo.D_i_down / mu; 
                f = calc_fanning_friction(Re);      
                dP_f = f * (geo.dz_down_v/geo.D_i_down) * (G_down^2)/(2*rho);
                dP_g = rho * 9.81 * geo.dz_down_v; 
                P_curr = P_curr + ((dP_g - dP_f) / 1e6); 
                Total_Mass = Total_Mass + rho * geo.A_c_down * geo.dz_down_v;
            else
                L_gas = L_unfilled - z_start; L_liq = z_end - L_unfilled;
                rho_g = 1 / safe_IAPWS('vV_p', P_curr); 
                mass_g = rho_g * geo.A_c_down * L_gas;
                [rho_l, mu_l, ~, ~, ~, ~, ~] = get_fluid_state(P_curr, h_curr);
                Re = G_down * geo.D_i_down / mu_l;
                f_l = calc_fanning_friction(Re);
                dP_f_l = f_l * (L_liq/geo.D_i_down) * (G_down^2)/(2*rho_l); 
                dP_g_l = rho_l * 9.81 * L_liq;
                Total_Mass = Total_Mass + mass_g + rho_l * geo.A_c_down * L_liq; 
                P_curr = P_curr + ((dP_g_l - dP_f_l) / 1e6);
            end
        end
        % 弯头1: 下竖 -> 下水
        [~, ~, ~, x_th, ~, ~, ~] = get_fluid_state(P_curr, h_curr);
        P_curr = P_curr - calc_local_dp(P_curr, x_th, G_down, geo.K_bend) / 1e6;
        
        % (2) 下降段 - 水平 (单根)
        for i = 1:geo.N_down_h
            [rho, mu, ~, x_th, ~, rho_mix, ~] = get_fluid_state(P_curr, h_curr);
            if (x_th > 0 && x_th < 1)
                rho_l = 1 / safe_IAPWS('vL_p', P_curr); 
                mu_l = safe_IAPWS('mu_ph', P_curr, safe_IAPWS('hL_p', P_curr));
                rho_g = 1 / safe_IAPWS('vV_p', P_curr); mu_g = safe_IAPWS('mu_ph', P_curr, safe_IAPWS('hV_p', P_curr));
                dP_f = calc_two_phase_friction_dp(x_th, G_down, geo.D_i_down, rho_l, mu_l, rho_g, mu_g, geo.dz_down_h);
            else
                Re = G_down * geo.D_i_down / mu; f = calc_fanning_friction(Re); dP_f = f * (geo.dz_down_h/geo.D_i_down) * (G_down^2)/(2*rho);
            end
            P_curr = P_curr - (dP_f / 1e6); Total_Mass = Total_Mass + rho_mix * geo.A_c_down * geo.dz_down_h;
        end
        % 分流集管/弯头: 下水 -> 蒸竖
        [~, ~, ~, x_th, ~, ~, ~] = get_fluid_state(P_curr, h_curr);
        P_curr = P_curr - calc_local_dp(P_curr, x_th, G_down, geo.K_bend) / 1e6;
        
        % (3) 蒸发段 - 竖直 (3根并联)
        for i = 1:geo.N_evap
            [~, ~, T_fluid, x_th, ~, rho_mix, mu_mix] = get_fluid_state(P_curr, h_curr);
            if x_th <= 0 || x_th >= 1
                Re = G_evap * geo.D_i_evap / mu_mix;
                Cp_fluid = safe_IAPWS('cp_ph', P_curr, h_curr); k_fluid = safe_IAPWS('k_ph', P_curr, h_curr);
                Pr = (Cp_fluid * 1000) * mu_mix / k_fluid;
                h_wall = safe_IAPWS('h_pT', P_curr, geo.T_wall_evap); mu_wall = safe_IAPWS('mu_ph', P_curr, h_wall); 
                z_curr = max(i * geo.dz_evap, 1e-4);
                h_htc = calc_single_phase_htc(Re, Pr, k_fluid, geo.D_i_evap, z_curr, mu_mix, mu_wall, true);
                U_i = 1 / (1 / h_htc + R_wall_evap); q_node = U_i * (geo.T_wall_evap - T_fluid);
            else
                q_guess = 5000; tol_q = 1e-3;
                for iter = 1:500
                    T_wall_in = geo.T_wall_evap - q_guess * R_wall_evap;
                    h_htc = calc_boiling_htc(P_curr, q_guess, x_th); 
                    q_new = h_htc * (T_wall_in - T_fluid); 
                    if abs(q_new - q_guess) / max(q_guess, 1) < tol_q, break; end
                    q_guess = 0.8 * q_guess + 0.2 * q_new; 
                end
                q_node = q_new;
            end
            dQ_single = q_node * pi * geo.D_i_evap * geo.dz_evap; 
            dQ_total = geo.N_tubes_evap * dQ_single;
            dh = (dQ_total / W) / 1000; 
            dP_g = rho_mix * 9.81 * geo.dz_evap;
            if (x_th > 0 && x_th < 1)
                rho_l = 1 / safe_IAPWS('vL_p', P_curr); mu_l = safe_IAPWS('mu_ph', P_curr, safe_IAPWS('hL_p', P_curr));
                rho_g = 1 / safe_IAPWS('vV_p', P_curr); mu_g = safe_IAPWS('mu_ph', P_curr, safe_IAPWS('hV_p', P_curr));
                dP_f = calc_two_phase_friction_dp(x_th, G_evap, geo.D_i_evap, rho_l, mu_l, rho_g, mu_g, geo.dz_evap);
            else
                Re = G_evap * geo.D_i_evap / mu_mix; f = calc_fanning_friction(Re); dP_f = f * (geo.dz_evap/geo.D_i_evap) * (G_evap^2)/(2*rho_mix);
            end
            P_curr = P_curr - ((dP_g + dP_f) / 1e6); h_curr = h_curr + dh; 
            Total_Mass = Total_Mass + geo.N_tubes_evap * rho_mix * geo.A_c_evap * geo.dz_evap;
        end
        % 汇流集管/弯头: 蒸竖 -> 升竖
        [~, ~, ~, x_th, ~, ~, ~] = get_fluid_state(P_curr, h_curr);
        P_curr = P_curr - calc_local_dp(P_curr, x_th, G_riser, geo.K_bend) / 1e6;

        % (4) 上升段 - 竖直 (单根)
        for i = 1:geo.N_riser_v
            [~, ~, ~, x_th, ~, rho_mix, mu_mix] = get_fluid_state(P_curr, h_curr);
            dP_g = rho_mix * 9.81 * geo.dz_riser_v;
            if (x_th > 0 && x_th < 1)
                rho_l = 1 / safe_IAPWS('vL_p', P_curr); mu_l = safe_IAPWS('mu_ph', P_curr, safe_IAPWS('hL_p', P_curr));
                rho_g = 1 / safe_IAPWS('vV_p', P_curr); mu_g = safe_IAPWS('mu_ph', P_curr, safe_IAPWS('hV_p', P_curr));
                dP_f = calc_two_phase_friction_dp(x_th, G_riser, geo.D_i_riser, rho_l, mu_l, rho_g, mu_g, geo.dz_riser_v);
            else
                Re = G_riser * geo.D_i_riser / mu_mix; f = calc_fanning_friction(Re); dP_f = f * (geo.dz_riser_v/geo.D_i_riser) * (G_riser^2)/(2*rho_mix);
            end
            P_curr = P_curr - ((dP_g + dP_f) / 1e6); Total_Mass = Total_Mass + rho_mix * geo.A_c_riser * geo.dz_riser_v;
        end
        % 分流集管/弯头: 升竖 -> 冷水
        [~, ~, ~, x_th, ~, ~, ~] = get_fluid_state(P_curr, h_curr);
        P_curr = P_curr - calc_local_dp(P_curr, x_th, G_riser, geo.K_bend) / 1e6;
        
        % (5) 冷凝段 - 水平 (3根并联)
        for i = 1:geo.N_cond
            [~, ~, T_fluid, x_th, ~, rho_mix, mu_mix] = get_fluid_state(P_curr, h_curr);
            if x_th > 0 && x_th < 1
                k_l = safe_IAPWS('k_ph', P_curr, safe_IAPWS('hL_p', P_curr)); rho_l = 1 / safe_IAPWS('vL_p', P_curr);
                mu_l = safe_IAPWS('mu_ph', P_curr, safe_IAPWS('hL_p', P_curr)); h_fg = (safe_IAPWS('hV_p', P_curr) - safe_IAPWS('hL_p', P_curr)) * 1000; 
                z_curr = max(i * geo.dz_cond, 1e-4); 
                q_guess = 5000; tol_q = 1e-3;
                for iter = 1:500
                    T_wall_in = geo.T_wall_cond + q_guess * R_wall_cond; dT_film = max(T_fluid - T_wall_in, 0.1); 
                    h_htc = 1.13 * ( (k_l^3 * rho_l^2 * h_fg * 9.81) / (mu_l * z_curr * dT_film) )^0.25; q_new = h_htc * dT_film;
                    if abs(q_new - q_guess) / max(q_guess, 1) < tol_q, break; end
                    q_guess = 0.8 * q_guess + 0.2 * q_new; 
                end
                q_node = q_new;
            else
                Re = G_cond * geo.D_i_cond / mu_mix; Cp_fluid = safe_IAPWS('cp_ph', P_curr, h_curr); k_fluid = safe_IAPWS('k_ph', P_curr, h_curr);
                Pr = (Cp_fluid * 1000) * mu_mix / k_fluid; T_wall_in_approx = max(T_fluid - 10, 273.16); 
                h_wall = safe_IAPWS('h_pT', P_curr, T_wall_in_approx); mu_wall = safe_IAPWS('mu_ph', P_curr, h_wall); 
                z_curr = max(i * geo.dz_cond, 1e-4); h_htc = calc_single_phase_htc(Re, Pr, k_fluid, geo.D_i_cond, z_curr, mu_mix, mu_wall, false);
                U_i = 1 / (1 / h_htc + R_wall_cond); q_node = U_i * max(T_fluid - geo.T_wall_cond, 0.1);
            end
            dQ_single = q_node * pi * geo.D_i_cond * geo.dz_cond; 
            dQ_total = geo.N_tubes_cond * dQ_single;
            dh = -(dQ_total / W) / 1000;
            
            if (x_th > 0 && x_th < 1)
                rho_l = 1 / safe_IAPWS('vL_p', P_curr); mu_l = safe_IAPWS('mu_ph', P_curr, safe_IAPWS('hL_p', P_curr));
                rho_g = 1 / safe_IAPWS('vV_p', P_curr); mu_g = safe_IAPWS('mu_ph', P_curr, safe_IAPWS('hV_p', P_curr));
                dP_f = calc_two_phase_friction_dp(x_th, G_cond, geo.D_i_cond, rho_l, mu_l, rho_g, mu_g, geo.dz_cond);
            else
                Re = G_cond * geo.D_i_cond / mu_mix; f = calc_fanning_friction(Re); dP_f = f * (geo.dz_cond/geo.D_i_cond) * (G_cond^2)/(2*rho_mix);
            end
            P_curr = P_curr - (dP_f / 1e6); h_curr = h_curr + dh; 
            Total_Mass = Total_Mass + geo.N_tubes_cond * rho_mix * geo.A_c_cond * geo.dz_cond;
        end
        % 汇流集管/弯头: 冷水 -> 下竖 (返回闭环起点)
        [~, ~, ~, x_th, ~, ~, ~] = get_fluid_state(P_curr, h_curr);
        P_curr = P_curr - calc_local_dp(P_curr, x_th, G_down, geo.K_bend) / 1e6;
        
        F = [(P_curr - P0) / P0, (h_curr - h0) / h0, (Total_Mass - M_charge) / M_charge];
    catch ME
        % 修复 M-Lint 警告：真实打印抛出错误信息，同时静默优化器中断
        fprintf(2, '求解器内部状态异常抛弃点: %s\n', ME.message);
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
function dP_f = calc_two_phase_friction_dp(x, G, D, rho_l, mu_l, rho_g, mu_g, dz)
    G_l = G * (1 - x); G_g = G * x; Re_l = max(G_l * D / mu_l, 1e-4); Re_g = max(G_g * D / mu_g, 1e-4);
    f_l = calc_fanning_friction(Re_l); f_g = calc_fanning_friction(Re_g); X = max(((1 - x) / x) * sqrt((f_l / f_g) * (rho_g / rho_l)), 1e-6);
    if Re_l < 2000 && Re_g < 2000, C = 5; elseif Re_l < 2000 && Re_g >= 2000, C = 12; elseif Re_l >= 2000 && Re_g < 2000, C = 10; else, C = 20; end
    dpdz_l = f_l * (G_l^2 / rho_l) / (2 * D); dP_f = (1 + C / X + 1 / X^2) * dpdz_l * dz;
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