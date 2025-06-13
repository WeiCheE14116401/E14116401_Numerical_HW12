import numpy as np
import math

# 常數定義
Nr = 6          # r: 0.5 to 1, step = 0.1
Nth = 6         # theta: 0 to pi/3, step = pi/15
r0 = 0.5
dr = 0.1
dtheta = math.pi / 15
tol = 1e-5
max_iter = 10000

def main():
    # 初始化溫度矩陣 T[i][j] = T(r_i, theta_j)
    T = np.zeros((Nr, Nth))
    
    # 初始化 r, theta 陣列
    r = np.array([r0 + i * dr for i in range(Nr)])
    theta = np.array([j * dtheta for j in range(Nth)])
    
    # 邊界條件設定
    # r = 0.5 和 r = 1.0 的邊界
    for j in range(Nth):
        T[0, j] = 50      # r = 0.5
        T[Nr-1, j] = 100  # r = 1.0
    
    # theta = 0 和 theta = pi/3 的邊界
    for i in range(Nr):
        T[i, 0] = 0       # theta = 0
        T[i, Nth-1] = 0   # theta = pi/3
    
    # Gauss-Seidel 迭代求解
    for iteration in range(max_iter):
        max_err = 0.0
        
        for i in range(1, Nr - 1):
            for j in range(1, Nth - 1):
                ri = r[i]
                
                # 計算各項係數
                term_r = (T[i+1, j] + T[i-1, j]) / (dr * dr)
                term_rr = (T[i+1, j] - T[i-1, j]) / (2 * dr)
                term_th = (T[i, j+1] + T[i, j-1]) / (dtheta * dtheta)
                
                # 新的溫度值計算
                T_new = (
                    term_r +
                    (1.0 / ri) * term_rr +
                    (1.0 / (ri * ri)) * term_th
                ) / (2.0 / (dr * dr) + 2.0 / (ri * ri * dtheta * dtheta))
                
                # 計算誤差並更新
                err = abs(T_new - T[i, j])
                T[i, j] = T_new
                if err > max_err:
                    max_err = err
        
        # 檢查收斂條件
        if max_err < tol:
            print(f"Converged in {iteration + 1} iterations.")
            break
    

    print("r\\theta", end="")
    for j in range(Nth):
        theta_deg = theta[j] * 180 / math.pi  # 轉換為度數顯示
        print(f"{theta_deg:8.1f}°\t", end="")
    print()
    
    print("-" * 50)
    
    # 輸出每一行的結果
    for i in range(Nr):
        print(f"{r[i]:6.1f} |", end="")
        for j in range(Nth):
            print(f"{T[i, j]:8.3f}\t", end="")
        print()
    
    print("=" * 50)

if __name__ == "__main__":
    main()
