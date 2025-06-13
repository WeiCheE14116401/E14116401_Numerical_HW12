import numpy as np

dr = 0.1
dt = 0.5
K = 0.1
N_r = 6     # 範圍 [0.5, 1.0]
N_t = 21
r_min = 0.5

alpha = 4 * K * dt / (dr**2)  # α = 4KΔt/Δr²

# 初始化網格
r = np.array([r_min + i*dr for i in range(N_r)])
T = np.zeros((N_t, N_r))

# 初始條件
T[0, :] = 200 * (r - 0.5)

for n in range(N_t-1):
    t = n * dt
    a = np.zeros(N_r-2)
    b = np.zeros(N_r-2)
    c = np.zeros(N_r-2)
    d = np.zeros(N_r-2)
    
    for i in range(1, N_r-1):
        ri = r[i]
        coeff1 = alpha
        coeff2 = alpha * dr / (2 * ri)
        
        a[i-1] = -coeff1 + coeff2
        b[i-1] = 1 + 2*coeff1
        c[i-1] = -coeff1 - coeff2
        d[i-1] = T[n, i]
    
    # 邊界條件處理
    T[n+1, 0] = T[n+1, 1] / (1 + 3*dr)
    d[0] -= a[0] * T[n+1, 0]
    
    T[n+1, -1] = 100 + 40*(t + dt)
    d[-1] -= c[-1] * T[n+1, -1]
    
    for i in range(1, N_r-2):
        m = a[i] / b[i-1]
        b[i] -= m * c[i-1]
        d[i] -= m * d[i-1]
    
    x = np.zeros(N_r-2)
    x[-1] = d[-1] / b[-1]
    for i in range(N_r-3, 0, -1):
        x[i-1] = (d[i-1] - c[i-1]*x[i]) / b[i-1]
    
    T[n+1, 1:-1] = x

header = [' t\\r'] + [f'{val:>8.4f}' for val in r]
print('\t'.join(header))
for j in range(N_t):
    current_t = 0.5 * j
    row = [f'{current_t:>4.1f}'] + [f'{T[j][i]:>8.4f}' for i in range(N_r)]
    print('\t'.join(row))


