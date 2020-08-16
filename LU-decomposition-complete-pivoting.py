import numpy as np

#########################################################################
#                   COMPLETE PIVOTING                                   #
#########################################################################

n = int(input("Please Enter The Dimension Of The Matrix: "))                  # ENTER MATRIX A
A = np.zeros((n, n, n))
for i in range(0, n):
    for j in range(0, n):
        print("A[", i, ",", j, "] :")
        A[0, i, j] = float(input(""))

b = np.zeros(n)                                                               # ENTER RHS VECTOR
print(" Please Enter The Elements Of b: ")
for i in range(0, n):
    print("b[", i, "]: ")
    b[i] = float(input(""))
print("--------------------------------------")
P = np.zeros((n-1, n, n))
Q = np.zeros((n-1, n, n))
M = np.zeros((n-1, n, n))
U = np.zeros((n, n))
m1 = 0
n1 = 0

for k in range(0, n-1):
    Max = abs(A[k, k, k])                                                 # FIND THE MAXIMUM VALUE IN MATRIX A[k]
    for i in range(k, n):
        for j in range(k, n):
            if abs(A[k, i, j]) > Max:
                Max = A[k, i, j]
                m1 = i                                                    # THE ROW WHICH INCLUDE MAXIMUM VALUE
                n1 = j                                                    # THE COLUMN WHICH INCLUDE MAXIMUM VALUE
    print("The Maximum Value In Step", k + 1, "Is:", A[k, m1, n1])
    print(" So Row", m1, "And Column ", n1, "Should Replace By Row And Column", k)

    P[k, :, :] = np.eye(n)                            # REPLACE ROWS
    temp1 = P[k, k, :].copy()
    P[k, k, :] = P[k, m1, :].copy()
    P[k, m1, :] = temp1.copy()
    print("The ROW Permutation Matrix(P) In Step", k + 1, "is:\n", P[k, :, :])

    Q[k, :, :] = np.eye(n)                            # REPLACE COLUMNS
    temp2 = Q[k, :, k].copy()
    Q[k, :, k] = Q[k, :, n1].copy()
    Q[k, :, n1] = temp2.copy()
    print("The Column Permutation Matrix(Q) In Step", k + 1, "is:\n", Q[k, :, :])

    A[k+1, :, :] = A[k, :, :]
    A[k+1, :, :] = P[k, :, :].dot(A[k, :, :]).dot(Q[k, :, :])             # A(k) = P(k)A(k-1)Q(k)
    M[k, :, :] = np.eye(n)
    M[k, k + 1:n, k] = -(A[k+1, k + 1:n, k] / A[k+1, k, k])               # COMPUTE PERMUTATION MATRIX M(k)
    print("Preliminary Matrix in step", k + 1, "is:\n", M[k, :, :])
    A[k+1, :, :] = M[k, :, :].dot(A[k+1, :, :])                           # M(k)p(k)A(k-1)Q(k)
    if k == n-2:
        U = A[k+1, :, :]
        print("Matrix A(n-1) = U is:\n", A[k+1, :, :])
    else:
        print("Matrix A in Step ", k + 1, "is:\n", A[k+1, :, :])

######################################################################################################
print("--------------------------------------")
# MAQ = U
# PAQ = LU
MM = np.eye(n)
for i in range(n-2, -1, -1):
    MM = MM.dot(M[i, :, :]).dot(P[i, :, :])                         # MM = M[n-1:0]P[n-1:0]
print("Permutation Matrix M:\n", MM)

QQ = np.eye(n)
for i in range(0, n-1):
    QQ = QQ.dot(Q[i, :, :])                                          # QQ = Q(1) ... Q(n-1)
print("The Final Permutation Matrix Q:\n", QQ)

PP = np.eye(n)                                                       # PP = P(n-1) ... P(1)
for i in range(n-2, -1, -1):
    PP = PP.dot(P[i, :, :])
print("The Final Permutation Matrix P:\n", PP)

L = np.dot(PP, np.linalg.inv(MM))                                         # L = P * M(inv)
print("The Unit Lower Triangular Matrix L Is:\n ", L)
print("The Upper Triangular Matrix U Is:\n ", U)

print("MAQ - U:\n", MM.dot(A[0, :, :].dot(QQ)) - U)                   # MAQ - U ~ 0
print("PAQ - LU:\n", (PP.dot(A[0, :, :])).dot(QQ) - L.dot(U))         # PAQ - LU = 0

#####################################################################################
#                        SOLVE SYSTEM Ax = b                                        #
#####################################################################################
print("--------------------------------------")
# MAQ = U
# SOLVE SYSTEM Ax = b USING MAQ = U
b1 = MM.dot(b)                                                       # b1 = Mb
y1 = np.zeros(n)
x1 = np.zeros(n)
s1 = 0
y1[n-1] = b1[n-1]/U[n-1, n-1]
for i in range(n-2, -1, -1):                                         # Uy = b1 (backward)
    for j in range(i+1, n):
        s1 += U[i, j] * y1[j]
    if U[i, i] != 0:
        y1[i] = (b1[i] - s1) / U[i, i]
        s1 = 0

x1 = QQ.dot(y1)                                                      # x = Qy
print("The System is Solved Using MAQ = U")
print("The Answer is:\n", x1)
print("--------------------------------------")
# PAQ = LU
# Ax = b
# SOLVE SYSTEM Ax = b USING PAQ = LU
b2 = PP.dot(b)                                                       # b2 = P*b
s2 = 0
s3 = 0
x2 = np.zeros(n)
z = np.zeros(n)
z[0] = b2[0]/L[0, 0]
y2 = np.zeros(n)
for i in range(1, n):                                                # Lz = b2 (forward)
    for j in range(1, i+1):
        s2 += L[i, j-1]*z[j-1]
    if L[i, i] != 0:
        z[i] = (b2[i] - s2) / L[i, i]
        s2 = 0

y2[n-1] = z[n-1]/U[n-1, n-1]
for i in range(n-2, -1, -1):                                        # Uy = z (backward)
    for j in range(i+1, n):
        s3 += U[i, j] * y2[j]
    if U[i, i] != 0:
        y2[i] = (z[i] - s3) / U[i, i]
        s3 = 0

x2 = QQ.dot(y2)                                                     # x = Qy
print("The System is Solved Using PAQ = LU")
print("The Answer is:\n", x2)
