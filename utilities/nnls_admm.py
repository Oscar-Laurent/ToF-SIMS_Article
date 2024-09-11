import numpy as np

def nnls_admm_py(A, b, Y, rho=1.0, max_iter=1000, tol=1e-8):
    """
    Solve `argmin_x ||Ax - b||_2` for `x >= 0` and `Ax <= Y` using ADMM.

    Parameters
    ----------
    A : (m, n) ndarray
        Coefficient array.
    b : (m,) ndarray
        Right-hand side vector.
    Y : (m,) ndarray
        Upper bound constraint for Ax.
    rho : float
        Augmented Lagrangian parameter.
    alpha : float
        Over-relaxation parameter (1.0 means no over-relaxation).
    max_iter : int
        Maximum number of iterations.
    tol : float
        Tolerance for the stopping criterion.

    Returns
    -------
    x : ndarray
        Solution vector.
    rnorm : float
        The 2-norm of the residual, `||Ax - b||_2`.
    """
    
    # Initialize variables
    m, n = A.shape
    x = np.zeros(n)
    z = np.zeros(m)
    u = np.zeros(m)
    
    # Precompute matrix factorizations
    AtA = A.T @ A
    I = np.eye(n)
    
    for k in range(max_iter):
        # x-update (using Cholesky decomposition)
        q = A.T @ (b + z - u) + rho * x
        L = np.linalg.cholesky(AtA + rho * I)
        x = np.linalg.solve(L.T, np.linalg.solve(L, q))
        x = np.maximum(x, 0)
        
        # z-update
        Ax = A @ x
        z_old = z
        z = np.minimum(Y, Ax + u)
        
        # u-update
        u += Ax - z
        
        # Check convergence
        rnorm = np.linalg.norm(Ax - b, 2)
        if np.linalg.norm(z - z_old) < tol:
            break
    
    return x, rnorm

