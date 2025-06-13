import numpy as np

# 常數定義
dr = 0.1
dt = 0.5
K = 0.1
N_r = 6         # 範圍 [0.5, 1.0]
N_t = 21
r_min = 0.5

alpha = 4 * K * dt / (dr * dr)  # α = 4KΔt / Δr²

def main():
    # 初始化半徑陣列
    r = np.array([r_min + i * dr for i in range(N_r)])
    
    # 初始化溫度矩陣
    T = np.zeros((N_t, N_r))
    
    # 初始條件
    for i in range(N_r):
        T[0, i] = 200 * (r[i] - 0.5)
    
    # Crank-Nicolson
    for n in range(N_t - 1):
        t = n * dt
        
        # 三對角係數
        a = np.zeros(N_r - 2)
        b = np.zeros(N_r - 2)
        c = np.zeros(N_r - 2)
        d = np.zeros(N_r - 2)
        
        for i in range(1, N_r - 1):
            ri = r[i]
            coeff1 = alpha / 2.0
            coeff2 = coeff1 * dr / (2 * ri)
            
            a[i - 1] = -coeff1 + coeff2
            b[i - 1] = 1 + 2 * coeff1
            c[i - 1] = -coeff1 - coeff2
            
            # RHS 計算 (forward half)
            term1 = alpha / 2 * (T[n, i + 1] - 2 * T[n, i] + T[n, i - 1])
            term2 = (alpha * dr / (4 * ri)) * (T[n, i + 1] - T[n, i - 1])
            d[i - 1] = T[n, i] + term1 + term2
        
        # 邊界條件
        T[n + 1, 0] = T[n + 1, 1] / (1 + 3 * dr)
        d[0] -= a[0] * T[n + 1, 0]
        T[n + 1, N_r - 1] = 100 + 40 * (t + dt)
        d[N_r - 3] -= c[N_r - 3] * T[n + 1, N_r - 1]
        

        for i in range(1, N_r - 2):
            m = a[i] / b[i - 1]
            b[i] -= m * c[i - 1]
            d[i] -= m * d[i - 1]
        
        x = np.zeros(N_r - 2)
        x[N_r - 3] = d[N_r - 3] / b[N_r - 3]
        for i in range(N_r - 4, -1, -1):
            x[i] = (d[i] - c[i] * x[i + 1]) / b[i]
        
        for i in range(1, N_r - 1):
            T[n + 1, i] = x[i - 1]
    
    print("T(r, t) 結果表格:")
    print("=" * 100)
    print(f"{'時間 t':>6} | {'r=0.5':>11} | {'r=0.6':>11} | {'r=0.7':>11} | {'r=0.8':>11} | {'r=0.9':>11} | {'r=1.0':>11} |")
    print("-" * 100)

    for j in range(N_t):
        t_val = j * dt
        print(f"{t_val:8.1f} |", end="")
        for i in range(N_r):
            print(f"{T[j, i]:12.4f} |", end="")
        print()
    
if __name__ == "__main__":
    main()
