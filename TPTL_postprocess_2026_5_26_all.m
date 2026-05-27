function TPTL_postprocess_2026_5_26_all(P0, H0, W)
% =========================================================================
% TPTL 稳态后处理函数（变径并联专版）
% 输入参数：P0 (MPa), H0 (m), W (总流量 kg/s)
% =========================================================================
    %% 1. 几何参数与边界条件
    geo.L_down_v = 3.0;  geo.N_down_v = 30;  geo.dz_down_v = geo.L_down_v / geo.N_down_v;
    geo.D_i_down = 0.05; geo.A_c_down = pi/4 * geo.D_i_down^2;

    geo.L_down_h = 3.0;  geo.N_down_h = 30;  geo.dz_down_h = geo.L_down_h / geo.N_down_h;

    geo.N_tubes_evap = 3;
    geo.L_evap = 2.83;   geo.N_evap = 50;    geo.dz_evap = geo.L_evap / geo.N_evap;
    geo.D_i_evap = 0.02; geo.D_o_evap = 0.025; geo.A_c_evap = pi/4 * geo.D_i_evap^2;

    geo.L_riser_v = 1.2; geo.N_riser_v = 20; geo.dz_riser_v = geo.L_riser_v / geo.N_riser_v;
    geo.D_i_riser = 0.05; geo.A_c_riser = pi/4 * geo.D_i_riser^2;

    geo.N_tubes_cond = 3;
    geo.L_cond = 3.0;    geo.N_cond = 50;    geo.dz_cond = geo.L_cond / geo.N_cond;
    geo.D_i_cond = 0.02; geo.D_o_cond = 0.025; geo.A_c_cond = pi/4 * geo.D_i_cond^2;

    geo.k_steel = 22.0; geo.K_bend = 0.5;
    geo.T_wall_evap = 273.15 + 250; geo.T_wall_cond = 273.15 + (20); 
    
    h0 = safe_IAPWS('hL_p', P0); 
    R_wall_evap = (geo.D_i_evap / (2 * geo.k_steel)) * log(geo.D_o_evap / geo.D_i_evap);
    R_wall_cond = (geo.D_i_cond / (2 * geo.k_steel)) * log(geo.D_o_cond / geo.D_i_cond);
    
    N_total = geo.N_down_v + geo.N_down_h + geo.N_evap + geo.N_riser_v + geo.N_cond;
    
    P_arr = zeros(N_total, 1); T_arr = zeros(N_total, 1); x_arr = zeros(N_total, 1);
    alpha_arr = zeros(N_total, 1); h_arr = zeros(N_total, 1); mass_arr = zeros(N_total, 1); dQ_arr = zeros(N_total, 1);
    
    idx = 0; Total_Mass = 0; P_curr = P0; h_curr = h0; Q_evap_total = 0; Q_cond_total = 0;
    
    dP_g_down_v_total  = 0; dP_g_evap_total    = 0; dP_g_riser_v_total = 0; dP_g_cond_total    = 0;
    dP_f_down_v_total  = 0; dP_f_down_h_total  = 0; dP_f_evap_total    = 0; dP_f_riser_v_total = 0; dP_f_cond_total    = 0;

    G_down  = W / geo.A_c_down;
    G_evap  = (W / geo.N_tubes_evap) / geo.A_c_evap;
    G_riser = W / geo.A_c_riser;
    G_cond  = (W / geo.N_tubes_cond) / geo.A_c_cond;

    %% 3. 逐段分布计算
    % ----------------------------- 下降段竖直 -----------------------------
    L_unfilled = max(0, geo.L_down_v - H0); z = 0;
    for i = 1:geo.N_down_v
        z_start = z; z = z + geo.dz_down_v; z_end = z;
        if L_unfilled >= z_end
            dP = 0; dh = 0; P_curr = P_curr + dP/1e6; h_curr = h_curr + dh;
            rho = 1 / safe_IAPWS('vV_p', P_curr); mass_inc = rho * geo.A_c_down * geo.dz_down_v;
            T_rec = safe_IAPWS('Tsat_p', P_curr); x_rec = 1; alpha_rec = 1; h_rec = safe_IAPWS('hV_p', P_curr);
            dP_g_node = rho * 9.81 * geo.dz_down_v; dP_f_node = 0;
        elseif L_unfilled <= z_start
            [rho, mu, ~, ~, ~, ~, ~] = get_fluid_state(P_curr, h_curr);
            Re = G_down * geo.D_i_down / mu; f = calc_fanning_friction(Re);
            dP_f = f * (geo.dz_down_v/geo.D_i_down) * (G_down^2)/(2*rho); dP_g = rho * 9.81 * geo.dz_down_v;
            dP = dP_g - dP_f; dh = 0; mass_inc = rho * geo.A_c_down * geo.dz_down_v;
            P_curr = P_curr + dP/1e6; h_curr = h_curr + dh;
            [~, ~, T_rec, x_rec, alpha_rec, ~, ~] = get_fluid_state(P_curr, h_curr); h_rec = h_curr;
            dP_g_node = dP_g; dP_f_node = dP_f;
        else
            L_gas = L_unfilled - z_start; L_liq = z_end - L_unfilled;
            rho_g = 1 / safe_IAPWS('vV_p', P_curr);
            [rho_l, mu_l, ~, ~, ~, ~, ~] = get_fluid_state(P_curr, h_curr);
            Re = G_down * geo.D_i_down / mu_l; f_l = calc_fanning_friction(Re);
            dP_f_l = f_l * (L_liq/geo.D_i_down) * (G_down^2)/(2*rho_l); dP_g_l = rho_l * 9.81 * L_liq;
            dP = dP_g_l - dP_f_l; dh = 0;
            mass_inc = rho_g * geo.A_c_down * L_gas + rho_l * geo.A_c_down * L_liq;
            P_curr = P_curr + dP/1e6; h_curr = h_curr + dh;
            [~, ~, T_rec, x_rec, ~, ~, ~] = get_fluid_state(P_curr, h_curr); alpha_rec = L_gas / geo.dz_down_v; h_rec = h_curr;
            dP_g_node = rho_g * 9.81 * L_gas + rho_l * 9.81 * L_liq; dP_f_node = dP_f_l;
        end
        dP_g_down_v_total = dP_g_down_v_total + dP_g_node; dP_f_down_v_total = dP_f_down_v_total + dP_f_node;
        Total_Mass = Total_Mass + mass_inc; idx = idx + 1;
        P_arr(idx) = P_curr; T_arr(idx) = T_rec; x_arr(idx) = x_rec; alpha_arr(idx) = alpha_rec; h_arr(idx) = h_rec; mass_arr(idx) = mass_inc; dQ_arr(idx) = 0;
    end
    % 【弯头1】: 垂直下降段 -> 水平下降段
    [~, ~, ~, x_th, ~, ~, ~] = get_fluid_state(P_curr, h_curr); 
    dP_local_bend1 = calc_local_dp(P_curr, x_th, G_down, geo.K_bend); P_curr = P_curr - dP_local_bend1 / 1e6;
    
    % ----------------------------- 下降段水平 -----------------------------
    for i = 1:geo.N_down_h
        [rho_s, mu_s, T_fluid, x_th, alpha, rho_mix, ~] = get_fluid_state(P_curr, h_curr);
        if (x_th > 0 && x_th < 1)
            rho_l = 1 / safe_IAPWS('vL_p', P_curr); mu_l = safe_IAPWS('mu_ph', P_curr, safe_IAPWS('hL_p', P_curr));
            rho_g = 1 / safe_IAPWS('vV_p', P_curr); mu_g = safe_IAPWS('mu_ph', P_curr, safe_IAPWS('hV_p', P_curr));
            dP_f = calc_two_phase_friction_dp(x_th, G_down, geo.D_i_down, rho_l, mu_l, rho_g, mu_g, geo.dz_down_h);
        else
            Re = G_down * geo.D_i_down / mu_s; f = calc_fanning_friction(Re); dP_f = f * (geo.dz_down_h/geo.D_i_down) * (G_down^2)/(2*rho_s);
        end
        dP_f_down_h_total = dP_f_down_h_total + dP_f; P_curr = P_curr - dP_f/1e6;
        mass_inc = rho_mix * geo.A_c_down * geo.dz_down_h; Total_Mass = Total_Mass + mass_inc;
        idx = idx + 1; P_arr(idx) = P_curr; T_arr(idx) = T_fluid; x_arr(idx) = x_th; alpha_arr(idx) = alpha; h_arr(idx) = h_curr; mass_arr(idx) = mass_inc; dQ_arr(idx) = 0;
    end
    % 【弯头2】: 水平下降段 -> 竖直蒸发段
    [~, ~, ~, x_th, ~, ~, ~] = get_fluid_state(P_curr, h_curr); 
    dP_local_bend2 = calc_local_dp(P_curr, x_th, G_down, geo.K_bend); P_curr = P_curr - dP_local_bend2 / 1e6;

    % ----------------------------- 蒸发段竖直 (3并) -----------------------------
    for i = 1:geo.N_evap
        [~, ~, T_fluid, x_th, alpha, rho_mix, mu_mix] = get_fluid_state(P_curr, h_curr);
        if x_th <= 0 || x_th >= 1
            Re = G_evap * geo.D_i_evap / mu_mix; Cp_fluid = safe_IAPWS('cp_ph', P_curr, h_curr); k_fluid = safe_IAPWS('k_ph', P_curr, h_curr);
            Pr = (Cp_fluid * 1000) * mu_mix / k_fluid; h_wall = safe_IAPWS('h_pT', P_curr, geo.T_wall_evap); mu_wall = safe_IAPWS('mu_ph', P_curr, h_wall);
            z_curr = max(i * geo.dz_evap, 1e-4); h_htc = calc_single_phase_htc(Re, Pr, k_fluid, geo.D_i_evap, z_curr, mu_mix, mu_wall, true);
            U_i = 1 / (1 / h_htc + R_wall_evap); q_node = U_i * (geo.T_wall_evap - T_fluid);
        else
            q_guess = 5000; tol_q = 1e-3;
            for iter = 1:500
                T_wall_in = geo.T_wall_evap - q_guess * R_wall_evap; h_htc = calc_boiling_htc(P_curr, q_guess, x_th);
                q_new = h_htc * (T_wall_in - T_fluid); if abs(q_new - q_guess) / max(q_guess, 1) < tol_q, break; end; q_guess = 0.8 * q_guess + 0.2 * q_new;
            end
            q_node = q_new;
        end
        dQ_single = q_node * pi * geo.D_i_evap * geo.dz_evap; dQ_total = geo.N_tubes_evap * dQ_single; 
        Q_evap_total = Q_evap_total + dQ_total; dh = (dQ_total / W) / 1000;
        dP_g = rho_mix * 9.81 * geo.dz_evap; dP_g_evap_total = dP_g_evap_total + dP_g;
        
        if (x_th > 0 && x_th < 1)
            rho_l = 1 / safe_IAPWS('vL_p', P_curr); mu_l = safe_IAPWS('mu_ph', P_curr, safe_IAPWS('hL_p', P_curr)); rho_g = 1 / safe_IAPWS('vV_p', P_curr); mu_g = safe_IAPWS('mu_ph', P_curr, safe_IAPWS('hV_p', P_curr));
            dP_f = calc_two_phase_friction_dp(x_th, G_evap, geo.D_i_evap, rho_l, mu_l, rho_g, mu_g, geo.dz_evap);
        else
            Re = G_evap * geo.D_i_evap / mu_mix; f = calc_fanning_friction(Re); dP_f = f * (geo.dz_evap/geo.D_i_evap) * (G_evap^2)/(2*rho_mix);
        end
        dP_f_evap_total = dP_f_evap_total + dP_f;
        P_curr = P_curr - (dP_g + dP_f)/1e6; h_curr = h_curr + dh; 
        mass_inc = geo.N_tubes_evap * rho_mix * geo.A_c_evap * geo.dz_evap; Total_Mass = Total_Mass + mass_inc;
        [~, ~, T_rec, x_rec, ~, ~, ~] = get_fluid_state(P_curr, h_curr);
        idx = idx + 1; P_arr(idx) = P_curr; T_arr(idx) = T_rec; x_arr(idx) = x_rec; alpha_arr(idx) = alpha; h_arr(idx) = h_curr; mass_arr(idx) = mass_inc; dQ_arr(idx) = dQ_total;
    end
    % (注意：蒸发竖直与上升竖直之间无弯头)

    % ----------------------------- 上升段竖直 (单根) -----------------------------
    for i = 1:geo.N_riser_v
        [~, ~, T_fluid, x_th, alpha, rho_mix, mu_mix] = get_fluid_state(P_curr, h_curr);
        dP_g = rho_mix * 9.81 * geo.dz_riser_v; dP_g_riser_v_total = dP_g_riser_v_total + dP_g;
        if (x_th > 0 && x_th < 1)
            rho_l = 1 / safe_IAPWS('vL_p', P_curr); mu_l = safe_IAPWS('mu_ph', P_curr, safe_IAPWS('hL_p', P_curr)); rho_g = 1 / safe_IAPWS('vV_p', P_curr); mu_g = safe_IAPWS('mu_ph', P_curr, safe_IAPWS('hV_p', P_curr));
            dP_f = calc_two_phase_friction_dp(x_th, G_riser, geo.D_i_riser, rho_l, mu_l, rho_g, mu_g, geo.dz_riser_v);
        else
            Re = G_riser * geo.D_i_riser / mu_mix; f = calc_fanning_friction(Re); dP_f = f * (geo.dz_riser_v/geo.D_i_riser) * (G_riser^2)/(2*rho_mix);
        end
        dP_f_riser_v_total = dP_f_riser_v_total + dP_f;
        P_curr = P_curr - (dP_g + dP_f)/1e6; mass_inc = rho_mix * geo.A_c_riser * geo.dz_riser_v; Total_Mass = Total_Mass + mass_inc;
        idx = idx + 1; P_arr(idx) = P_curr; T_arr(idx) = T_fluid; x_arr(idx) = x_th; alpha_arr(idx) = alpha; h_arr(idx) = h_curr; mass_arr(idx) = mass_inc; dQ_arr(idx) = 0;
    end
    
    % 【弯头3】: 垂直上升段 -> 水平冷凝段
    [~, ~, ~, x_th, ~, ~, ~] = get_fluid_state(P_curr, h_curr); 
    dP_local_bend3 = calc_local_dp(P_curr, x_th, G_riser, geo.K_bend); P_curr = P_curr - dP_local_bend3 / 1e6;

    % ----------------------------- 冷凝段水平 (3并) 无重力 -----------------------------
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
                if abs(q_new - q_guess) / max(q_guess, 1) < tol_q, break; end; q_guess = 0.8 * q_guess + 0.2 * q_new; 
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
        Q_cond_total = Q_cond_total + abs(dQ_total); dh = -(dQ_total / W) / 1000;
        
        dP_g_cond_total = dP_g_cond_total + 0; 
        
        if (x_th > 0 && x_th < 1)
            rho_l = 1 / safe_IAPWS('vL_p', P_curr); mu_l = safe_IAPWS('mu_ph', P_curr, safe_IAPWS('hL_p', P_curr));
            rho_g = 1 / safe_IAPWS('vV_p', P_curr); mu_g = safe_IAPWS('mu_ph', P_curr, safe_IAPWS('hV_p', P_curr));
            dP_f = calc_two_phase_friction_dp(x_th, G_cond, geo.D_i_cond, rho_l, mu_l, rho_g, mu_g, geo.dz_cond);
        else
            Re = G_cond * geo.D_i_cond / mu_mix; f = calc_fanning_friction(Re); dP_f = f * (geo.dz_cond/geo.D_i_cond) * (G_cond^2)/(2*rho_mix);
        end
        dP_f_cond_total = dP_f_cond_total + dP_f;
        P_curr = P_curr - (dP_f / 1e6); h_curr = h_curr + dh; 
        mass_inc = geo.N_tubes_cond * rho_mix * geo.A_c_cond * geo.dz_cond; Total_Mass = Total_Mass + mass_inc;
        [~, ~, T_rec, x_rec, alpha, ~, ~] = get_fluid_state(P_curr, h_curr);
        idx = idx + 1; P_arr(idx) = P_curr; T_arr(idx) = T_rec; x_arr(idx) = x_rec; alpha_arr(idx) = alpha; h_arr(idx) = h_curr; mass_arr(idx) = mass_inc; dQ_arr(idx) = -dQ_total; 
    end
    
    % 【弯头4】: 水平冷凝段 -> 垂直下降段 (返回闭环起点)
    [~, ~, ~, x_th, ~, ~, ~] = get_fluid_state(P_curr, h_curr); 
    dP_local_bend4 = calc_local_dp(P_curr, x_th, G_down, geo.K_bend); P_curr = P_curr - dP_local_bend4 / 1e6;

    %% 4. 终端文字输出
    Q_avg = (Q_evap_total + Q_cond_total) / 2;   
    fprintf('\n====== 函数沿程计算完成 ======\n');
    fprintf('初始 P0 = %.4f MPa, 最终 P = %.4f MPa\n', P0, P_curr);
    fprintf('初始 h0 = %.4f kJ/kg, 最终 h = %.4f kJ/kg\n', h0, h_curr);
    fprintf('压力相对误差 = %.4e\n', (P_curr - P0)/P0);
    fprintf('比焓相对误差 = %.4e\n', (h_curr - h0)/h0);
    fprintf('回路累计总质量 = %.4f kg\n', Total_Mass);
    fprintf('蒸发段总吸热量 = %.2f W\n', Q_evap_total);
    fprintf('冷凝段总放热量 = %.2f W\n', Q_cond_total);
    fprintf('热管平均传热量 = %.2f W\n', Q_avg);

    Net_Gravity_Head = (dP_g_down_v_total + dP_g_cond_total) - (dP_g_evap_total + dP_g_riser_v_total);
    Sum_Friction_DP = dP_f_down_v_total + dP_f_down_h_total + dP_f_evap_total + dP_f_riser_v_total + dP_f_cond_total;
    Sum_Local_DP = dP_local_bend1 + dP_local_bend2 + dP_local_bend3 + dP_local_bend4;

    fprintf('\n====== 1. 循环重力压降 (压头) 统计 ======\n');
    fprintf('垂直下降段重力压降 (驱动力) = %8.2f Pa\n', dP_g_down_v_total);
    fprintf('水平冷凝段重力压降 (无重力) = %8.2f Pa\n', dP_g_cond_total);
    fprintf('竖直蒸发段重力压降 (阻力)   = %8.2f Pa\n', dP_g_evap_total);
    fprintf('竖直上升段重力压降 (阻力)   = %8.2f Pa\n', dP_g_riser_v_total);
    fprintf('----------------------------------------\n');
    fprintf('回路净重力驱动压头 (总和)   = %8.2f Pa\n', Net_Gravity_Head);

    fprintf('\n====== 2. 各管段沿程摩擦压降统计 ======\n');
    fprintf('垂直下降段沿程摩擦压降 = %8.2f Pa\n', dP_f_down_v_total);
    fprintf('水平下降段沿程摩擦压降 = %8.2f Pa\n', dP_f_down_h_total);
    fprintf('竖直蒸发段沿程摩擦压降 = %8.2f Pa\n', dP_f_evap_total);
    fprintf('竖直上升段沿程摩擦压降 = %8.2f Pa\n', dP_f_riser_v_total);
    fprintf('水平冷凝段沿程摩擦压降 = %8.2f Pa\n', dP_f_cond_total);
    fprintf('----------------------------------------\n');
    fprintf('回路摩擦压降总和       = %8.2f Pa\n', Sum_Friction_DP);

    fprintf('\n====== 3. 各弯头局部阻力压降统计 ======\n');
    fprintf('弯头1 (垂直下降 -> 水平下降) = %8.2f Pa\n', dP_local_bend1);
    fprintf('弯头2 (水平下降 -> 竖直蒸发) = %8.2f Pa\n', dP_local_bend2);
    fprintf('弯头3 (垂直上升 -> 水平冷凝) = %8.2f Pa\n', dP_local_bend3);
    fprintf('弯头4 (水平冷凝 -> 垂直下降) = %8.2f Pa\n', dP_local_bend4);
    fprintf('----------------------------------------\n');
    fprintf('回路局部压降总和             = %8.2f Pa\n', Sum_Local_DP);
    fprintf('========================================\n');

    %% 5. 保存数据到 Excel 文件 
    Th_deg = geo.T_wall_evap - 273.15; Tl_deg = geo.T_wall_cond - 273.15;
    Node_Indices = (1:N_total)';
    T_excel = table(Node_Indices, T_arr, P_arr, x_arr, alpha_arr, mass_arr, dQ_arr, h_arr, ...
        'VariableNames', {'节点编号', '温度_K', '压力_MPa', '干度', '空泡率', '节点质量_kg', '节点传热量_W', '比焓_kJ_kg'});
    sheet_name = sprintf('P%.2fH%.2fW%.4fTh%.0fTl%.0f', P0, H0, W, Th_deg, Tl_deg);
    if length(sheet_name) > 31, sheet_name = sheet_name(1:31); end
    excel_filename = 'TPTL_Nodes_Data.xlsx';
    writetable(T_excel, excel_filename, 'Sheet', sheet_name);

    %% 6. 绘图 
    node_sec = [0, geo.N_down_v, geo.N_down_v + geo.N_down_h, ...
                geo.N_down_v + geo.N_down_h + geo.N_evap, ...
                geo.N_down_v + geo.N_down_h + geo.N_evap + geo.N_riser_v, N_total];
    
    figure('Position', [50, 100, 1500, 650]); nodes_axis = 1:N_total;
    
    subplot(2,4,1); plot(nodes_axis, P_arr*1e6, 'b-o', 'MarkerSize', 3, 'LineWidth', 1); hold on;
    for k = 2:length(node_sec)-1, xline(node_sec(k), '--k', 'Alpha', 0.4); end; title('压力沿程分布'); grid on;
    
    subplot(2,4,2); plot(nodes_axis, T_arr, 'r-o', 'MarkerSize', 3, 'LineWidth', 1); hold on;
    for k = 2:length(node_sec)-1, xline(node_sec(k), '--k', 'Alpha', 0.4); end; title('温度沿程分布'); grid on;
    
    subplot(2,4,3); plot(nodes_axis, x_arr, 'g-o', 'MarkerSize', 3, 'LineWidth', 1); hold on;
    for k = 2:length(node_sec)-1, xline(node_sec(k), '--k', 'Alpha', 0.4); end; title('干度沿程分布'); grid on; ylim([-0.1, 1.1]);
    
    subplot(2,4,4); plot(nodes_axis, alpha_arr, 'c-o', 'MarkerSize', 3, 'LineWidth', 1); hold on;
    for k = 2:length(node_sec)-1, xline(node_sec(k), '--k', 'Alpha', 0.4); end; title('空泡率沿程分布'); grid on; ylim([-0.1, 1.1]);
    
    subplot(2,4,5); plot(nodes_axis, h_arr, 'm-o', 'MarkerSize', 3, 'LineWidth', 1); hold on;
    for k = 2:length(node_sec)-1, xline(node_sec(k), '--k', 'Alpha', 0.4); end; title('比焓沿程分布'); grid on;
    
    subplot(2,4,6); plot(nodes_axis, mass_arr, 'k-o', 'MarkerSize', 3, 'LineWidth', 1.5); hold on;
    for k = 2:length(node_sec)-1, xline(node_sec(k), '--k', 'Alpha', 0.4); end; title('节点局部质量分布'); grid on;
    
    subplot(2,4,7); plot(nodes_axis, dQ_arr, 'r-o', 'MarkerSize', 3, 'LineWidth', 1); hold on;
    for k = 2:length(node_sec)-1, xline(node_sec(k), '--k', 'Alpha', 0.4); end; title('节点净传热量分布'); grid on;

    subplot(2,4,8); axis off;
    text(0.05, 0.5, sprintf(['系统运行参数:\n  P0=%.4f MPa, H0=%.3f m\n  总流量W=%.4f kg/s, 充注M=%.4f kg\n\n', ...
        '传热表现:\n  吸热量=%.1f W\n  放热量=%.1f W\n  平均热量=%.1f W\n\n', ...
        '回路压降平衡:\n  净驱动压头 = %.1f Pa\n  沿程总摩擦 = %.1f Pa\n  局部总弯头 = %.1f Pa'], ...
        P0, H0, W, Total_Mass, Q_evap_total, Q_cond_total, Q_avg, Net_Gravity_Head, Sum_Friction_DP, Sum_Local_DP), ...
        'FontSize', 9.5, 'VerticalAlignment','middle', 'FontName', '微軟正黑體');

    sgtitle(sprintf('TPTL 变径异构环路综合性能曲线 (P0=%.3f MPa, H0=%.3f m)', P0, H0));
end

%% ===================== 局部子函数区 =====================
function f = calc_fanning_friction(Re)
    Re = max(Re, 1e-4);
    if Re < 2000, f = 64 / Re; else, f = 0.184 / Re^0.2; end
end
function dP_f = calc_two_phase_friction_dp(x, G, D, rho_l, mu_l, rho_g, mu_g, dz)
    G_l = G * (1 - x); G_g = G * x;
    Re_l = max(G_l * D / mu_l, 1e-4); Re_g = max(G_g * D / mu_g, 1e-4);
    f_l = calc_fanning_friction(Re_l); f_g = calc_fanning_friction(Re_g);
    X = ((1 - x) / x) * sqrt((f_l / f_g) * (rho_g / rho_l)); X = max(X, 1e-6);
    if Re_l < 2000 && Re_g < 2000, C = 5; elseif Re_l < 2000 && Re_g >= 2000, C = 12;
    elseif Re_l >= 2000 && Re_g < 2000, C = 10; else, C = 20; end
    dpdz_l = f_l * (G_l^2 / rho_l) / (2 * D); dP_f = (1 + C / X + 1 / X^2) * dpdz_l * dz;
end
function dP_local = calc_local_dp(P, x, G, K_bend)
    rho_l = 1 / safe_IAPWS('vL_p', P); 
    rho_g = 1 / safe_IAPWS('vV_p', P); 
    dP_l = K_bend * (G^2) / (2 * rho_l);
    
    if x <= 0
        % 单相液态
        dP_local = dP_l; 
    elseif x >= 1
        % 单相气态
        dP_local = K_bend * (G^2) / (2 * rho_g);
    else
        % 两相混合态
        x_safe = max(min(x, 0.999), 1e-4); 
        rho_v = rho_g; % 匹配你公式中的气相密度符号
        
        % ======= 原来的两相计算公式（已注释） =======
        Xb = (rho_g/rho_l)^0.5 * (1-x_safe)/x_safe; 
        C = (1 + (3.2-1)*(rho_g)^0.5) * ((rho_g/rho_l)^0.5 + (rho_g/rho_l)^-0.5); 
        phi_c_sq = 1 + C/Xb + 1/Xb^2; 
        dP_local = dP_l * phi_c_sq;
        % ========================================

        % ======= 新的两相计算公式 =======
        % MF = G^2 / rho_l * (rho_l/rho_v - 1) * x_safe * (1 - x_safe) * (1.1 / (2 + 1.5));
        % dP_local = 1.978 * G^2 / (2 * rho_l) * (1 + x_safe * (rho_l/rho_v - 1)) + MF;
    end
end
function [rho_single, mu_single, T_fluid, x_th, alpha, rho_mix, mu_mix] = get_fluid_state(P, h)
    hL = safe_IAPWS('hL_p', P); hV = safe_IAPWS('hV_p', P); Tsat = safe_IAPWS('Tsat_p', P); x_th = (h - hL) / (hV - hL);
    if x_th <= 0
        x_th = 0; T_fluid = safe_IAPWS('T_ph', P, h); rho_single = 1 / safe_IAPWS('v_ph', P, h); mu_single = safe_IAPWS('mu_ph', P, h); alpha = 0; rho_mix = rho_single; mu_mix = mu_single;
    elseif x_th >= 1
        x_th = 1; T_fluid = safe_IAPWS('T_ph', P, h); rho_single = 1 / safe_IAPWS('v_ph', P, h); mu_single = safe_IAPWS('mu_ph', P, h); alpha = 1; rho_mix = rho_single; mu_mix = mu_single;
    else
        T_fluid = Tsat; vL = safe_IAPWS('vL_p', P); vV = safe_IAPWS('vV_p', P); rho_l = 1/vL; rho_g = 1/vV;
        mu_l = safe_IAPWS('mu_ph', P, hL + 0.1); mu_g = safe_IAPWS('mu_ph', P, hV - 0.1);
        alpha = 1 / (1 + ((1-x_th)/x_th) * (rho_g/rho_l)^0.89 * (mu_l/mu_g)^0.18); alpha = max(0, min(1, alpha));
        rho_mix = alpha * rho_g + (1 - alpha) * rho_l; mu_mix = 1 / ( (x_th/mu_g) + ((1-x_th)/mu_l) ); rho_single = rho_l; mu_single = mu_l;
    end
end
function val = safe_IAPWS(prop, varargin)
    val = IAPWS_IF97(prop, varargin{:}); if any(isnan(val)), error('TPTL:NaN_Error', 'IAPWS_IF97 out of bounds.'); end
end
function h_htc = calc_boiling_htc(P_MPa, q, x)
    M = 18.015; P_crit = 22.064; Pr = max(P_MPa / P_crit, 1e-4); q = max(q, 10); x = max(0, min(x, 0.99));
    h_htc = 540 * Pr^0.394 * (1-x)^-0.65 * M^(-0.5) * q^0.54;
end
function h_htc = calc_single_phase_htc(Re, Pr, k, D, z, mu_fluid, mu_wall, is_heating)
    Re = max(Re, 10); Pr = max(Pr, 0.1); z = max(z, 1e-4);
    if Re < 2200, F1 = 1 + 5 * exp(-z / (10 * D)); h_htc = 1.24 * (k * F1 / D) * (Re * Pr * (D / z))^(1/3) * (mu_fluid / mu_wall);
    else, F2 = 1 + 7 * (D / z); n = 0.4; if ~is_heating, n = 0.3; end; h_htc = 0.027 * (k * F2 / D) * Re^0.8 * Pr^n * (mu_fluid / mu_wall)^0.14; end
end