import numpy as np

#########################################################################
#              ECONOMICAL PARTIAL PIVOTING                             #
#########################################################################

# MAQ = U
n = int(input("Please Enter The Dimension Of The Matrix: "))
A0 = np.zeros((n, n))
for i in range(0, n):
    for j in range(0, n):
        print("A[", i, ",", j, "] :")
        A0[i, j] = int(input(""))

b = np.zeros(n)
print(" Please Enter The Elements Of b: ")
for i in range(0, n):
    print("b[", i, "]: ")
    b[i] = float(input(""))
print("--------------------------------------")
row = np.zeros(n)
A = A0.copy()
for k in range(0, n-1):
    Max = abs(A[k, k])                                             # FIND THE MAXIMUM VALUE IN MATRIX A[k:n,k]
    for i in range(k, n):
        if abs(A[i, k]) > Max:
            Max = A[i, k]
            row[k] = i                                            # ROW WHICH INCLUDE THE MAXIMUM VALUE
    print("The Maximum Value In Step ", k+1, "Is:", A[int(row[k]), k])
    print("So Row ", k, "Should Replace By Row", row[k])

    temp1 = A[k, k:n].copy()                                       # REPLACE ROWS
    A[k, k:n] = A[int(row[k]), k:n].copy()
    A[int(row[k]), k:n] = temp1.copy()
    print("Matrix A After Pivoting In Step ", k+1, "Is:\n", A)

    for i in range(k+1, n):
        A[i, k] = A[i, k] / A[k, k]

    for i in range(k+1, n):
        for j in range(k+1, n):
            A[i, j] = A[i, j] - A[i, k]*A[k, j]
    if k == n-1:
        print("Matrix A(n-1) = U is:\n", A)
    else:
        print("A After Step ", k+1, "Is:\n", A)


############################################################################################################
M = np.zeros((n-1, n, n))                                            # PRELIMINARY Matrix M(k)
for i in range(0, n-1):
    M[i, :, :] = np.eye(n)
    M[i, i+1:n, i] = -A[i+1:n, i]
    print(" The Preliminary Matrix M In Step ", i+1, "Is:\n", M[i, :, :])

P = np.zeros((n-1, n, n))                                           # PERMUTATION MATRIX P(k)
for r in range(0, n-1):
    P[r, :, :] = np.eye(n)
    t2 = P[r, r, :].copy()
    P[r, r, :] = P[r, int(row[r]), :].copy()
    P[r, int(row[r]), :] = t2.copy()
    print("Row Permutation Matrix P In Step", r+1, "Is:\n", P[r, :, :])
PP = np.eye(n)
for i in range(n-2, -1, -1):                                          # PP = P(n-1) ... P(1)
    PP = PP.dot(P[i, :, :])
print("Final Permutation Matrix P is:\n", PP)

U = np.triu(A)                                                       # UPPER TRIANGULAR MATRIX U
print("The Upper Triangular Matrix U is:\n", U)

MM = np.eye(n)
for i in range(n-2, -1, -1):
    MM = MM.dot(M[i, :, :].dot(P[i, :, :]))                          # MM = M(n-1)P(n-1) ... M(1)P(1)
print("Preliminary Matrix Is;\n", MM)

L = np.dot(PP, np.linalg.inv(MM))                                # UNIT LOWER TRIANGULAR MATRIX L (L = PM(inv))
print("The Unit Triangular Matrix L Is:\n", L)

print("MA - U:\n", np.dot(MM, A0) - U)
print("PA - LU:\n", np.dot(PP, A0) - np.dot(L, U))
#####################################################################################
print("--------------------------------------")
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
print("--------------------------------------")
# PA = LU
# Ax = b
b2 = np.dot(PP, b)                                                    # b2 = pb
print("b2:\n", b2)
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
x2[n-1] = y[n-1]/U[n-1, n-1]
for i in range(n-2, -1, -1):                                        # Ux = y (backward)
    for j in range(i+1, n):
        s3 += U[i, j] * x2[j]
    if U[i, i] != 0:
        x2[i] = (y[i] - s3) / U[i, i]
        s3 = 0
print("The System is Solved Using PAQ = LU")
print("The Answer is:\n", x2)

