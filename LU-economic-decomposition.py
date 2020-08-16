import numpy as np

# -------------------------------------#
#            A=LU  economic           #
# -------------------------------------#
A = np.array([[1, 4, 7], [2, 5, 8], [3, 6, 10]])
n, _ = np.shape(A)

A0 = A.copy()

for k in range(0, n - 1):
    for i in range(k + 1, n):
        A[i, k] = A[i, k] / A[k, k]

    for i in range(k + 1, n):
        for j in range(k + 1, n):
            A[i, j] = A[i, j] - A[i, k] * A[k, j]
U = np.triu(A)
L = np.tril(A)
L = L - np.diag(np.diag(L)) + np.eye(n)

print('==========================================')
print('Old A:\n', A0)
print('\nnew A:\n', A)

print('\nU:\n', U)
print('\nL:', L)
print('\nA-LU:')
print(A0 - np.dot(L, U))

print('==========================================')
M_stack = np.zeros((n - 1, n, n))
for k in range(0, n - 1):
    M_stack[k, :, :] = np.eye(n)
    for i in range(k + 1, n):
        M_stack[k, i, k] = -A[i, k]

M = np.eye(n, n)
for k in range(0, n - 1):
    T = M_stack[k, :, :].dot(M[:, :])
    M[:, :] = T[:, :]

for k in range(1, n):
    print('M' + str(k) + ':\n', M_stack[k - 1, :, :])
    print('------------------------------------------')

print('\nM:\n', M)

M_inv = np.linalg.inv(M)
print('\nM_inv:\n', M_inv)
