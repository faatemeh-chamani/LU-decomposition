import numpy as np

# -------------------------------------#
#              economic               #
#            MA=U & PA=LU             #
# -------------------------------------#
A_tmp = np.array([[1, 4, 7], [2, 5, 8], [3, 6, 10]])
n, _ = np.shape(A_tmp)

A = np.zeros(np.shape(A_tmp))
A[:, :] = A_tmp[:, :]

row = np.eye(1, n - 1)

for k in range(0, n - 1):
    r = k
    Max = np.abs(A[k, k])
    for i in range(k + 1, n):
        if np.abs(A[i, k]) > Max:
            Max = np.abs(A[i, k])
            r = i
    row[0][k] = r

    f = A[r, k:].copy()
    A[r, k:] = A[k, k:].copy()
    A[k, k:] = f.copy()

    for i in range(k + 1, n):
        A[i, k] = A[i, k] / A[k, k]
        for j in range(k + 1, n):
            A[i, j] = A[i, j] - A[i, k] * A[k, j]
print('==========================================')
print('Old A:\n', A_tmp)
print('\nnew A:\n', A)
print('\nrow:\n', row)
print('==========================================')
P_stack = np.zeros((n - 1, n, n))
M_stack = np.zeros((n - 1, n, n))

for k in range(n - 1):
    P_stack[k, :, :] = np.eye(n)
    f = P_stack[k, int(row[0][k]), :].copy()
    P_stack[k, int(row[0][k]), 0:] = P_stack[k, k, :].copy()
    P_stack[k, k, :] = f.copy()

    M_stack[k, :, :] = np.eye(n)
    for i in range(k + 1, n):
        M_stack[k, i, k] = -A[i, k]

M = np.eye(n, n)
for k in range(0, n - 1):
    T = P_stack[k, :, :].dot(M)
    M = M_stack[k, :, :].dot(T)

for k in range(1, n):
    print('P' + str(k) + ':\n', P_stack[k - 1, :, :])
    print('\nM' + str(k) + ':\n', M_stack[k - 1, :, :])
    print('==========================================')

U = np.triu(A)

# ==================================================== #
print('M:\n', M)
print('\nU:\n', U)
print('MA-U:\n', M.dot(A_tmp) - U[:, :])
print('==========================================')
# ==================================================== #
P = np.eye(n, n)
for k in range(0, n - 1):
    P = P_stack[k, :, :].dot(P)

L = P.dot(np.linalg.inv(M))

print('P:\n', P)
print('L:\n', L)
print('PA-LU:\n', P.dot(A_tmp) - L.dot(U))
