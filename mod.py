#!/usr/bin/python
#-*- coding: utf-8 -*-
"""求解多元单模线性方程组
   参考：https://github.com/55-AA/mod_equations
"""
import copy

# 最大公约数（欧几里得）
def gcd(a, b):
    while b:
        a, b = b, a % b
    return a

# 最小公倍数
def lcm(a, b):
    return a * b / (gcd(a, b))

# 扩展欧几里得
def egcd(a, b):
    """
    ax + by = 1
    ax ≡ 1 mod b
    Return a 3-element tuple (g, x, y), the g  = gcd(a, b)
    gmpy2.gcdext(a, b)
    """
    if a == 0:
        return (b, 0, 1)
    else:
        g, y, x = egcd(b % a, a)
        return (g, x - (b // a) * y, y)


def mod_inv(a, m):
    """
    ax ≡ 1 mod m
    gmpy2.invert(a, m)
    """
    g, x, y = egcd(a, m)
    assert g == 1
    return x % m


# class
class Matrix:
    """
    A*X ≡ B (mod p),p为大于0的整数。
    高斯消元求解模线性方程组。先化简为上三角，然后回代求解。
    当r(A) <= n时，一定有多解；
    当r(A) == n时，有多解或唯一解；
    当r(A) != r(A~)时，无解。
    r(A)为系数矩阵的秩，r(A)为增广矩阵的秩，n为未知数的个数。
    http://www.docin.com/p-1063811671.html讨论了gcd(|A|, m) = 1时的LU分解解法，
    本文包括了gcd(|A|, m) > 1时的解法，
    化简原则：
        1、系数与模互质
        2、系数加某一行n次后，对应的系数与模的GCD最小
        3、将1或2得到的系数移到对角线上
    初始化参数：
        matrix：方程组的增广矩阵（最后一列为常数项）。
            matrix = [
                [ 69,  75,  78,  36,  58],
                [ 46,  68,  51,  26,  42],
                [ 76,  40,  42,  49,  11],
                [ 11,  45,   2,  45,   1],
                [ 15,  67,  60,  14,  72],
                [ 76,  67,  73,  56,  58],
                [ 67,  15,  68,  54,  75],
            ]    
        mod：模数
    函数：
        gauss()：求解方程
    输出变量：
        error_str：出错的信息
        count：解的数量
    """

    def __init__(self, matrix, mod):
        self.matrix = copy.deepcopy(matrix)
        self.d = None
        self.row = len(matrix)
        self.col = len(matrix[0]) - 1
        self.N = len(matrix[0])
        self.mod = mod
        self.count = 1
        self.error_str = "unknown error"

    def swap_row(self, ra, rb):
        (self.d[ra], self.d[rb]) = (self.d[rb], self.d[ra])

    def swap_col(self, ca, cb):
        for j in range(self.row):
            (self.d[j][ca], self.d[j][cb]) = (self.d[j][cb], self.d[j][ca])

    def inv_result(self, r, n):
        """
        求解第n个未知数，r已经获得的解。形如：[None,None, ..., n+1, ...]
        a*x ≡ b(mod m)
        x有解的条件：gcd(a,m) | b。也即a,m互质时一定有解，不互质时，b整除gcd(a,m)也有解，否则无解。
        解的格式为：x0+k(m/gcd(a,m))，其中x0为最小整数特解，k为任意整数。
        返回[x0, x1, ...xn]，其中x0 < x1 < xn < m。
        """
        b = self.d[n][self.col]
        a = self.d[n][n]
        m = self.mod
        k = gcd(a, m)
        for j in xrange(n + 1, self.col):
            b = (b - (self.d[n][j] * r[j] % m)) % m

        if 1 == k:
            return [mod_inv(a, m) * b % m]
        else:
            if k == gcd(k, b):
                a /= k
                b /= k
                m /= k

                x0 = mod_inv(a, m) * b % m
                x = []
                for i in xrange(k):
                    x.append(x0 + m * i)
                return x
        return None

    def find_min_gcd_row_col(self, i, j):
        # 查找直接互质的对角线系数
        for k in xrange(i, self.row):
            for l in xrange(j, self.N - 1):
                if (1 == gcd(self.d[k][l], self.mod)):
                    return [k, l]

        def add_min_gcd(a, b, m):
            r = [m, 1]
            g = gcd(a, b)
            if g:
                i = a / g
                for j in xrange(i):
                    g = gcd((a + j * b) % m, m)
                    if g < r[0]:
                        r[0] = g
                        r[1] = j
                    if g == 1:
                        break
            return r

        # 查找加乘后GCD最小的对角线系数
        # [加乘后的最大公约数,加乘的倍数,要化简的行号,加乘的行号,要化简的列号]
        r = [self.mod, 1, i, i + 1, j]
        for k in xrange(i, self.row):
            for kk in xrange(k + 1, self.row):
                for l in range(j, self.N - 1):
                    rr = add_min_gcd(self.d[k][l], self.d[kk][l], self.mod)
                    if rr[0] < r[0]:
                        r[0] = rr[0]
                        r[1] = rr[1]
                        r[2] = k
                        r[3] = kk
                        r[4] = l
                        pass
                    if (1 == rr[0]):
                        break
        g = r[0]
        n = r[1]
        k = r[2]
        kk = r[3]
        l = r[4]

        if n and g < self.mod:
            self.d[k] = map(lambda x, y: (x + n * y) % self.mod, self.d[k], self.d[kk])
        return [k, l]

    def mul_row(self, i, k, j):
        a = self.d[k][j]
        b = self.d[i][j]

        def get_mul(a, b, m):
            k = gcd(a, m)
            if 1 == k:
                return mod_inv(a, m) * b % m
            else:
                if k == gcd(k, b):
                    return mod_inv(a / k, m / k) * (b / k) % (m / k)
            return None

        if b:
            mul = get_mul(a, b, self.mod)
            if None == mul:
                print_matrix(self.d)
                assert (mul != None)
            self.d[i] = map(lambda x, y: (y - x * mul) % self.mod, self.d[k], self.d[i])

    # 求解方程
    def gauss(self):
        """
        返回解向量，唯一解、多解或无解(None)。
        例如：[[61, 25, 116, 164], [61, 60, 116, 94], [61, 95, 116, 24], [61, 130, 116, 129], [61, 165, 116, 59]]
        """

        self.d = copy.deepcopy(self.matrix)
        for i in xrange(self.row):
            for j in xrange(self.N):
                self.d[i][j] = self.matrix[i][j] % self.mod

        if self.row < self.col:
            self.d.extend([[0] * self.N] * (self.col - self.row))

        # 化简上三角
        index = [x for x in xrange(self.col)]
        for i in range(self.col):
            tmp = self.find_min_gcd_row_col(i, i)
            if (tmp):
                self.swap_row(i, tmp[0])
                (index[i], index[tmp[1]]) = (index[tmp[1]], index[i])
                self.swap_col(i, tmp[1])
            else:
                self.error_str = "no min"
                return None

            for k in range(i + 1, self.row):
                self.mul_row(k, i, i)

        # print_matrix(self.d)
        if self.row > self.col:
            for i in xrange(self.col, self.row):
                for j in xrange(self.N):
                    if self.d[i][j]:
                        self.error_str = "r(A) != r(A~)"
                        return None

        # 判断解的数量
        for i in xrange(self.col):
            self.count *= gcd(self.d[i][i], self.mod)

        if self.count > 100:
            self.error_str = "solution too more:%d" % (self.count)
            return None

        # 回代
        result = [[None] * self.col]
        for i in range(self.col - 1, -1, -1):
            new_result = []
            for r in result:
                ret = self.inv_result(r, i)
                if ret:
                    for rr in ret:
                        l = r[:]
                        l[i] = rr
                        new_result.append(l)

                else:
                    self.error_str = "no inv:i=%d" % (i)
                    return None

            result = new_result

        # 调整列变换导致的未知数顺序变化
        for i in xrange(len(result)):
            def xchg(a, b):
                result[i][b] = a

            map(xchg, result[i][:], index)

        return result


# 打印行
def print_array(x):
    prn = "\t["
    for j in x:
        if j:
            prn += "%3d, " % j
        else:
            prn += "  0, "

    print prn[:-2] + " ],"

# 打印矩阵
def print_matrix(x):
    print "["
    for i in x:
        print_array(i)
    print "]"

def run_test(mod, matrix):
    print "row = %d, col = %d" % (len(matrix), len(matrix[0]) - 1)
    print "mod = %d" % (mod)
    print "Augmented matrix ="
    print_matrix(matrix)

    g = Matrix(matrix, mod)

    ret = g.gauss()
    if not ret:
        print "error:"
        print_matrix(g.d)
        print "error_str:", g.error_str
    else:
        print "counts:", g.count
        print "result:"
        print_matrix(ret)

def random_test():
    import random

    #  随机生成初始Test值，随机范围可根据情况修改
    row = random.randint(2, 3)
    col = random.randint(2, 3)
    mod = random.randint(5, 10)
    matrix = [ # 为增广矩阵

    ]
    for y in xrange(row):
        array = map(lambda x: random.randint(-mod, mod), [xc for xc in xrange(col)])
        array.append(random.randint(1-mod, mod - 1))
        matrix.append(array)

    # 手动输入静态Test值
    # row =
    # col =
    # mod =
    # matrix = [ ]

    run_test(mod, matrix)

if __name__ == "__main__":
    random_test()
exit(0)