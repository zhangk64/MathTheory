#!/usr/bin/python
#-*- coding: utf-8 -*-
#矩阵乘法计算

#欧几里得
def gcd(a, b):
    while b:
        a, b = b, a % b
    return a

#扩展欧几里得
def egcd(a, b):
    # a*x + b*y = r
    if b == 0:
        return (1, 0, a)
    (x, y, r) = egcd(b, a%b)

    tmp = x
    x = y
    y = tmp - (a/b) * y
    return (x, y, r)

# 求二阶阶幺模矩阵
def Unitary_2():
    import random
    a = random.randint(-10,10)
    b = random.randint(-10,10)
    if a == 0:
        while b==0:
            b = random.randint(-10, 10)
    (x,y,r) = egcd(a,b)

    a = a/r
    b = -b/r

    #打印生成的二阶幺模矩阵
    # print (a, b)
    # print (y, x)
    #print r

    return a,b,y,x

# 求n阶幺模矩阵
def Unitary_n(n):
    U = []
    for i in range(n/2):
        a1, b1, a2, b2= Unitary_2()
        list1 =  [0 for x in range(n)]
        list1[2*i],list1[2*i+1] = a1,b1
        list2 =  [0 for x in range(n)]
        list2[2 * i], list2[2 * i + 1] = a2, b2
        U.append(list1)
        U.append(list2)
    return U

# 随机生成矩阵A
def Matrix_A(n):
    import random
    A = []
    for i in range(n):
        list =  map(lambda x: random.randint(-10, 10), [x for x in xrange(n)])
        A.append(list)
    return A

# E乘以q
def Matrix_Eq(n, q):
    Eq = []
    for i in range(n):
        list = [0 for x in range(n)]
        list[i] = q
        Eq.append(list)
    return Eq


#矩阵加法 M + M
def madd(M1, M2):
    if isinstance(M1, (tuple, list)) and isinstance(M2, (tuple, list)):
        return [[m+n for m,n in zip(i,j)] for i, j in zip(M1,M2)]

#矩阵乘法 M x M
def multi(M1, M2):
    if isinstance(M1, (float, int)) and isinstance(M2, (tuple, list)):
        return [[M1*i for i in j] for j in M2]
    if isinstance(M2, (float, int)) and isinstance(M1, (tuple, list)):
        return [[M2*i for i in j] for j in M1]
    if isinstance(M1, (tuple, list)) and isinstance(M2, (tuple, list)):
        n = len(M1)
        Mres = [[0] * n for row in range(n)]  # 初始化c为n行n列的全零矩阵
        for i in range(0, n):
            for j in range(0, n):
                for k in range(0, n):
                    if M1[i][k] and M2[k][j]: # 此处0跳过计算
                        Mres[i][j] = Mres[i][j] + M1[i][k] * M2[k][j]
        return Mres

#矩阵乘以常数
def multi_q(M, q):
    M1 = []
    for i in M:
        list =  map(lambda x: x*q, [x for x in i])
        M1.append(list)
    return M1

#打印行
def print_array(x):
    prn = "\t["
    for j in x:
        if j:
            prn += "%3d, " % j
        else:
            prn += "  0, "

    print prn[:-2] + " ],"

#打印矩阵
def print_matrix(M):
    print "["
    for i in M:
        print_array(i)
    print "]"

#计算表达式
def calculation():
    import random

    n = input("Please input a even number:")
    if (n%2  == 1):
        print "'n' is not even number. Please input a even number!"
        return;
    # 系数q
    q =  random.randint(-100,100)
    r = random.randint(-10, 10)

    U1 = Unitary_n(n)
    U2 = Unitary_n(n)
    U3 = Unitary_n(n)
    A  = Matrix_A(n)
    Eq = Matrix_Eq(n, q)


    #计算A*U1+ q*E*U2=？
    print "1.计算表达式A*U1+ q*E*U2=？"
    print "-----------U1-----------"
    print_matrix(U1)
    print "-----------A-----------"
    print_matrix(A)
    print "-----------U2-----------"
    print_matrix(U2)
    print "-----------Eq-----------"
    print_matrix(Eq)

    M_result1 = madd(multi(A, U1), multi(Eq, U2))
    print "-----A*U1+ q*E*U2计算结果-----"
    print_matrix(M_result1)

    #计算A*U3=？
    print "\n2.计算表达式A*U3=？"
    print "-----------U3-----------"
    print_matrix(U3)
    M_result2 = multi(A, U3)
    print "--------A*U3计算结果--------"
    print_matrix(M_result2)

    #计算A×r=？
    print "\n3.计算表达式A*r=？"
    print "------------r------------ "
    print" r  = ",   r
    print "--------A*r计算结果--------"
    M_result3 = multi_q(A, r)
    print_matrix(M_result3)

if __name__ == "__main__":
    calculation()
exit(0)