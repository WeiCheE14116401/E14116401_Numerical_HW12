import numpy as np
import pandas as pd

def solve_heat_equation():

    # 參數設定
    dr = 0.1
    dt = 0.5
    K = 0.1
    N_r = 6         # 範圍 [0.5, 1.0]，共 6 點 (包含端點)
    N_t = 21        # t 從 [0,10]，Δt=0.5 → 21 步
    r_min = 0.5
    r_max = 1.0
    alpha = 20  
    
    print("參數設定:")
    print(f"dr = {dr}, dt = {dt}, K = {K}")
    print(f"N_r = {N_r}, N_t = {N_t}")
    print(f"r範圍: [{r_min}, {r_max}]")
    print(f"alpha = {alpha}")
    
    # 建立 r 座標
    r = np.zeros(N_r)
    for i in range(N_r):
        r[i] = r_min + i * dr
    
    # 建立 T 矩陣 T[時間][空間]
    T = np.zeros((N_t, N_r))
    
    # 初始條件 T(r,0) = 200(r - 0.5)
    for i in range(N_r):
        T[0][i] = 200 * (r[i] - 0.5)
    
    # 時間迴圈
    for n in range(N_t - 1):
        t = n * dt
        
        # 更新內部節點 i = 1 到 N_r - 2
        for i in range(1, N_r - 1):
            r_i = r[i]
            term1 = T[n][i + 1] - 2 * T[n][i] + T[n][i - 1]
            term2 = (T[n][i + 1] - T[n][i]) * (dr / r_i)
            T[n + 1][i] = T[n][i] + alpha * (term1 + term2)
        
        # 邊界條件：右邊 T(1, t) = 100 + 40t
        T[n + 1][N_r - 1] = 100 + 40 * (t + dt)
        
        # 邊界條件：左邊 ∂T/∂r + 3T = 0
        T[n + 1][0] = T[n + 1][1] / (1 + 3 * dr)
    
    return r, T

# 執行計算
r, T = solve_heat_equation()

# 表格顯示結果
print("結果顯示")
print("="*80)
time_indices = [0, 4, 8, 12, 16, 20]
evolution_data = {'r': r}
for t_idx in time_indices:
    t_val = t_idx * 0.5
    evolution_data[f't={t_val:.1f}'] = T[t_idx]

evolution_table = pd.DataFrame(evolution_data)
print(evolution_table.to_string(index=False, float_format='%.4f'))

