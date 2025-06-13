import math
# 常數
Nx = 11        # x: 0.0 to 1.0, step = 0.1
Nt = 21        # 模擬 t=0 ~ 2.0, step = 0.1
dx = 0.1
dt = 0.1
pi = math.pi
lambda2 = (dt * dt) / (dx * dx)  # λ^2 = 1

def initial_u(x):
    return math.cos(2 * pi * x)

def initial_ut(x):
    return 2 * pi * math.sin(2 * pi * x)

def print_horizontal_table(time_values, x_values, u_values_matrix):

    col_width = 12    # 計算每個欄位的寬度

    header = "x/t".ljust(col_width)     # 輸出表頭 - 時間值
    for t in time_values:
        header += f"| t={t:.1f}".ljust(col_width)
    print(header)

    separator = "-" * col_width     # 輸出分隔線
    for _ in range(len(time_values)):
        separator += "+" + "-" * (col_width - 1)
    print(separator)   

    for i, x_val in enumerate(x_values):     # 輸出每一行 (對應每個x值)
        row = f"x={x_val:.1f}".ljust(col_width)
        for t_idx in range(len(time_values)):
            row += f"| {u_values_matrix[t_idx][i]:.5f}".ljust(col_width)
        print(row)
# 初始化陣列
u_prev = [0.0] * Nx  # p^0
u_curr = [0.0] * Nx  # p^1
u_next = [0.0] * Nx  # p^{n+1}
x = [i * dx for i in range(Nx)]
# 初始狀態 t = 0
for i in range(Nx):
    u_prev[i] = initial_u(x[i])
# 邊界條件 t=0
u_prev[0] = 1.0
u_prev[Nx - 1] = 2.0

# 使用 Taylor 展開計算 t = dt 的值（u_curr）
for i in range(1, Nx - 1):
    utt = (u_prev[i+1] - 2*u_prev[i] + u_prev[i-1]) / (dx*dx)
    u_curr[i] = u_prev[i] + dt * initial_ut(x[i]) + 0.5 * dt * dt * utt

# 邊界值
u_curr[0] = 1.0
u_curr[Nx - 1] = 2.0

# 儲存所有時間點的結果
all_time_values = [0.0, 0.1]  # 初始包含 t=0 和 t=0.1
all_u_values = [u_prev.copy(), u_curr.copy()]  # 儲存每個時間點的u值

for n in range(2, 11):
    current_time = n * dt
    all_time_values.append(current_time)
    
    # 計算內部節點
    u_next = [0.0] * Nx
    for i in range(1, Nx - 1):
        u_next[i] = (2 * u_curr[i] - u_prev[i] + 
                    lambda2 * (u_curr[i+1] - 2 * u_curr[i] + u_curr[i-1]))
    
    u_next[0] = 1.0
    u_next[Nx - 1] = 2.0
    
    all_u_values.append(u_next.copy())
    u_prev = u_curr.copy()
    u_curr = u_next.copy()

print_horizontal_table(all_time_values, x, all_u_values)
print("=" * 130)
