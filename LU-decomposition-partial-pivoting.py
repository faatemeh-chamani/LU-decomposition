import numpy as np
##########################################################################
#                       PARTIAL PIVOTING                                 #
##########################################################################

n = int(input("Please Enter The Dimension Of The Matrix: "))
A = np.zeros((n, n, n))
for i in range(0, n):
    for j in range(0, n):
        print("A[", i, ",", j, "] :")
        A[0, i, j] = float(input(""))

b = np.zeros(n)
print("Please Enter the Elements of Right Hand Side Vector:")
for i in range(0, n):
    print("b[", i, "]:")
    b[i] = int(input(""))

x = np.zeros(n)
P = np.zeros((n-1, n, n))
M = np.zeros((n-1, n, n))
m = 0

for k in range(0, n - 1):
    Max = abs(A[k, k, k])
    for i in range(k + 1, n):
        if abs(A[k, i, k]) > Max:
            Max = A[k, i, k]
            m = i
    print("Max value in step", k + 1, "is : ", A[k, m, k])
    print("row", k, "should replace by row", m)

    P[k, :, :] = np.eye(n, n)
    temp = P[k][m, :].copy()
    P[k][m, :] = P[k][k, :].copy()
    P[k][k, :] = temp.copy()
    print("The Permutation Matrix P(k) In Step", k + 1, "Is: \n", P[k, :, :])

    A[k+1, :, :] = P[k, :, :].dot(A[k, :, :])
    M[k, :, :] = np.eye(n, n)
    for j in range(k + 1, n):
        M[k, j, k] = -float(A[k+1, j, k] / A[k+1, k, k])
    print("Matrix M include Multipliers in step", k + 1, "is:\n", M[k, :, :])

    A[k + 1, :, :] = M[k, :, :].dot(A[k+1, :, :])                        # A(K) = M(K) * P(K) * A(K-1)
    print("Matrix A in step", k + 1, "is:\n", A[k+1, :, :])

# U = A(n-1)
U = A[n-1, :, :]
print("The Upper Triangular Matrix U is:\n", U)
#######################################################################################

PP = np.eye(n)
for i in range(n-2, -1, -1):
    PP = PP.dot(P[i, :, :])         # P = P(n-1) ... P(1)
print("Row Permutation Matrix Is:\n", PP)

MM = np.eye(n)
for i in range(n-2, -1, -1):
    MM = MM.dot(M[i, :, :]).dot(P[i, :, :])         # M = M(n-1)P(n-1) ... M(1)P(1)
print("Final Preliminary Matrix Is;\n", MM)
L = np.dot(PP, np.linalg.inv(MM))
print("The Unit Lower Triangular Matrix L:\n", L)
print("MA - U:\n", MM.dot(A[0, :, :]) - U)
print("PA - LU:\n", np.dot(PP, A[0, :, :]) - np.dot(L, U))
#####################################################################################
# MA = U
# Ax = b
b1 = np.dot(MM, b)                                                 # b1 = Mb
print("b1 is:\n", b1)
x1 = np.zeros(n)
s1 = 0
x1[n-1] = b1[n-1]/U[n-1, n-1]
for i in range(n-2, -1, -1):                                         # Ux = b1 (backward)
    for j in range(i+1, n):
        s1 += U[i, j] * x1[j]
    if U[i, i] != 0:
        x1[i] = (b1[i] - s1) / U[i, i]
        s1 = 0
print("The System is Solved Using MAQ = U")
print("The Answer is:\n", x1)

# PA = LU
# Ax = b
b2 = np.dot(PP, b)                                                    # b2 = pb
s2 = 0
s3 = 0
x2 = np.zeros(n)

y = np.zeros(n)
y[0] = b2[0]/L[0, 0]
for i in range(1, n):                                                # Ly = b2 (forward)
    for j in range(1, i+1):
        s2 += L[i, j-1]*y[j-1]
    if L[i, i] != 0:
        y[i] = (b2[i] - s2) / L[i, i]
        s2 = 0
print("y", y)
x2[n-1] = y[n-1]/U[n-1, n-1]
for i in range(n-2, -1, -1):                                        # Ux = y (backward)
    for j in range(i+1, n):
        s3 += U[i, j] * x2[j]
    if U[i, i] != 0:
        x2[i] = (y[i] - s3) / U[i, i]
        s3 = 0
print("The System is Solved Using PAQ = LU")
print("The Answer is:\n", x2)
