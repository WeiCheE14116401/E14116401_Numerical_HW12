import numpy as np

def solve_pde():
    Nx = 11  # x: 0 to π (Δx=0.1π*10)
    Ny = 6   # y: 0 to π/2 (Δy=0.1π*5)
    h = 0.1 * np.pi
    tolerance = 1e-6
    max_iter = 10000

    def f(x, y):
        return x * y

    u = np.zeros((Nx, Ny))
    x = np.array([i * h for i in range(Nx)])
    y = np.array([j * h for j in range(Ny)])

    # BC
    for j in range(Ny):
        u[0, j] = np.cos(y[j])
        u[Nx-1, j] = -np.cos(y[j])
    for i in range(Nx):
        u[i, 0] = np.cos(x[i])
        u[i, Ny-1] = 0.0

    converged = False
    for iter in range(max_iter):
        max_error = 0.0
        for i in range(1, Nx-1):
            for j in range(1, Ny-1):
                u_old = u[i, j]
                u[i, j] = 0.25 * (u[i+1,j] + u[i-1,j] + u[i,j+1] + u[i,j-1] - h**2 * f(x[i], y[j]))
                error = abs(u[i, j] - u_old)
                max_error = max(max_error, error)
        if max_error < tolerance:
            print(f"Converged in {iter+1} iterations.")
            converged = True
            break

    #print u(x, y)
    print("\nu(x, y) values:")
    for j in (range(Ny)):
        print(" ".join([f"{u[i,j]:10.6f}" for i in range(Nx)]))
        
    #print u_x,y
    print("\nCorresponding function indices:")
    for j in (range(Ny)):
        print(" ".join([f"u_{i:02d},{j:02d}" for i in range(Nx)]))

if __name__ == "__main__":
    solve_pde()
