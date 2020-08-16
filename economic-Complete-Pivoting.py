import numpy as np

#########################################################################
#              ECONOMICAL COMPLETE PIVOTING                             #
#########################################################################

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

row = np.zeros(n)
column = np.zeros(n)
A = A0.copy()
for k in range(0, n-1):
    Max = abs(A[k, k])                                             # FIND THE MAXIMUM VALUE IN MATRIX A[k]
    for i in range(k, n):
        for j in range(k, n):
            if abs(A[i, j]) > Max:
                Max = A[i, j]
                row[k] = i                                         # ROW WHICH INCLUDE THE MAXIMUM VALUE
                column[k] = j                                      # COLUMN WHICH INCLUDE THE MAXIMUM VALUE
    file = open("file.text", 'a')
    file.write('The Maximum Value In Step' + str(k+1) + ' Is: ' + str(A[int(row[k]), int(column[k])]))
    file.write('\n So Row And Column ' + str(k) + ' Should Replace By Row' + str(row[k]) + 'Column' + str(column[k]))
    file.close()

    temp1 = A[k, k:n].copy()                                        # REPLACE ROWS
    A[k, k:n] = A[int(row[k]), k:n].copy()
    A[int(row[k]), k:n] = temp1.copy()

    temp2 = A[:, k].copy()                                         # REPLACE COLUMNS
    A[:, k] = A[:, int(column[k])].copy()
    A[:, int(column[k])] = temp2.copy()

    file = open("file.text", 'a')
    file.write('\n Matrix A After Pivoting In Step ' + str(k+1) + 'Is: \n' + str(A))
    file.close()
    for i in range(k+1, n):
        A[i, k] = A[i, k] / A[k, k]

    for i in range(k+1, n):
        for j in range(k+1, n):
            A[i, j] = A[i, j] - A[i, k] * A[k, j]
    if k == n-2:
        file = open("file.text", 'a')
        file.write('\n Matrix A(n-1) = U is:\n' + str(A))
        file.close()
    else:
        file = open("file.text", 'a')
        file.write('\n A After Step ' + str(k+1) + ' is:' + '\n' + str(A))
        file.close()

############################################################################################################
M = np.zeros((n-1, n, n))                                                   # PRELIMINARY Matrix M(k)
for i in range(0, n-1):
    M[i, :, :] = np.eye(n)
    M[i, i+1:n, i] = -A[i+1:n, i]
    file = open("file.text", 'a')
    file.write('\n The Preliminary Matrix M In Step ' + str(i+1) + 'Is:\n' + str(M[i, :, :]))
    file.close()

Q = np.zeros((n-1, n, n))                                                   # PERMUTATION MATRIX Q(k)
QQ = np.eye(n)
for r in range(0, n-1):
    Q[r, :, :] = np.eye(n)
    t1 = Q[r, :, r].copy()
    Q[r, :, r] = Q[r, :, int(column[r])].copy()
    Q[r, :, int(column[r])] = t1.copy()
    file = open("file.text", 'a')
    file.write('\n Column Permutation Matrix Q In Step' + str(r+1) + 'Is:\n' + str(Q[r, :, :]))
    file.close()
    QQ = QQ.dot(Q[r, :, :])
file = open("file.text", 'a')
file.write('\n Final Permutation Matrix Q is:\n' + str(QQ))                  # QQ = Q(1) ... Q(n-1)
file.close()

P = np.zeros((n-1, n, n))                                                    # PERMUTATION MATRIX P(k)
for r in range(0, n-1):
    P[r, :, :] = np.eye(n)
    t2 = P[r, r, :].copy()
    P[r, r, :] = P[r, int(row[r]), :].copy()
    P[r, int(row[r]), :] = t2.copy()
    file = open("file.text", 'a')
    file.write('\n Row Permutation Matrix P In Step' + str(r+1) + ' Is:\n' + str(P[r, :, :]))
    file.close()
PP = np.eye(n)
for i in range(n-2, -1, -1):                                                  # PP = P(n-1) ... P(1)
    PP = PP.dot(P[i, :, :])
file = open("file.text", 'a')
file.write('\n Final Permutation Matrix P is:\n' + str(PP))
file.close()

U = np.triu(A)                                                                # UPPER TRIANGULAR MATRIX U
file = open("file.text", 'a')
file.write('\n The Upper Triangular Matrix U is:\n' + str(U))
file.close()

MM = np.eye(n)
for i in range(n-2, -1, -1):                                                 # MM = M(n-1)P(n-1) ... M(1)P(1)
    MM = MM.dot(M[i, :, :]).dot(P[i, :, :])
file = open("file.text", 'a')
file.write('\nFinal Permutation Matrix M Is:\n' + str(MM))
file.close()

L = np.dot(PP, np.linalg.inv(MM))                     # UNIT LOWER TRIANGULAR MATRIX L
file = open("file.text", 'a')
file.write('\n The Unit Lower Triangular Matrix L is:\n' + str(L))
file.close()

file = open("file.text", 'a')
file.write('\n MAQ - U:\n' + str(MM.dot(A0.dot(QQ)) - U))                         # MAQ - U ~ 0
file.write('\n PAQ - LU:\n' + str((PP.dot(A0.dot(QQ)) - L.dot(U))))               # PAQ - LU = 0
file.close()
#####################################################################################
#                        SOLVE SYSTEM Ax = b                                        #
#####################################################################################
# MAQ = U
# Ax = b
b1 = MM.dot(b)                                  # b1 = Mb
y1 = np.zeros(n)
x1 = np.zeros(n)
s1 = 0
y1[n-1] = b1[n-1]/U[n-1, n-1]
for i in range(n-2, -1, -1):                   # Uy = b1 (backward)
    for j in range(i+1, n):
        s1 += U[i, j] * y1[j]
    if U[i, i] != 0:
        y1[i] = (b1[i] - s1) / U[i, i]
        s1 = 0

x1 = QQ.dot(y1)                                   # x = Qy
file = open("file.text", 'a')
file.write('\nThe System is Solved Using MAQ = U\n')
file.write('The Answer is:\n' + str(x1))
file.close()

# PAQ = LU
# Ax = b
# SOLVE SYSTEM Ax = b USING PAQ = LU
b2 = PP.dot(b)                                                # b2 = P*b
s2 = 0
s3 = 0
x2 = np.zeros(n)
z = np.zeros(n)
z[0] = b2[0] / L[0, 0]
y2 = np.zeros(n)
y2[n - 1] = z[n - 1] / U[n - 1, n - 1]
for i in range(1, n):                                          # Lz = b2 (forward)
    for j in range(1, i + 1):
        s2 += L[i, j - 1] * z[j - 1]
    if L[i, i] != 0:
        z[i] = (b2[i] - s2) / L[i, i]
        s2 = 0

y2[n-1] = z[n-1] / U[n-1, n-1]
for i in range(n - 2, -1, -1):                                # Uy = z (backward)
    for j in range(i + 1, n):
        s3 += U[i, j] * y2[j]
    if U[i, i] != 0:
        y2[i] = (z[i] - s3) / U[i, i]
        s3 = 0

x2 = QQ.dot(y2)                                               # x = Qy
file = open("file.text", 'a')
file.write('\nThe System is Solved Using PAQ = LU\n')
file.write('The Answer is:\n' + str(x2))
file.close()

print(" 'All Computations Are Stored In A File'  ")
