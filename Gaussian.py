# Gaussian Elimination

class Matrix:
    def __init__(self, shape: tuple, values: list[float]):
        if len(values) != shape[0] * shape[1]:
            raise ValueError("incorrect number of elements")
        self.n = shape[0]
        self.m = shape[1]
        self.values = values

    def getAt(self, i, j):
        return self.values[self.m * i + j]

    def getRow(self, i):
        return Vector('r', self.m, self.values[i * self.m: i * self.m + self.m])

    def setRow(self, i, row):
        self.values[i * self.m: i * self.m + self.m] = row.values

    def getCol(self, i):
        return Vector('c', self.n, self.values[i: self.n * self.m: self.m])

    def setCol(self, i, col):
        A[i: self.n * self.m: self.m] = col.values
    
    def __str__(self):
        s = "["
        s += self.getRow(0).__str__() + '\n'
        for i in range(1, self.n):
            s += " " + self.getRow(i).__str__()
            s += '\n'
        s = s.strip('\n')
        s += ']'
        return s

    def addCol(self, col):
        k = self.m
        self.m += 1
        for i in col.values:
            self.values.insert(k, i)
            k += self.m

    def elimination(self):
        def red(pivotAt):
            if pivotAt != (self.n, self.m - 1):
                row = pivotAt[0]
                col = pivotAt[1]
                piv = self.getAt(row, col)
                if piv == 0:
                    r1 = self.getRow(row)
                    r2 = self.getRow(row + 1)
                    self.setRow(row, r2)
                    self.setRow(row + 1, r1)
                    piv = self.getAt(row, col)
                row += 1
                while row < self.n:
                    k = self.getAt(row, col) / piv
                    self.setRow(row, self.getRow(row) - self.getRow(pivotAt[0]) * k)
                    row += 1
                red((pivotAt[0] + 1, pivotAt[1] + 1))
        red((0,0))

class Vector:
    def __init__(self, type: str, n: int, values: list[float]):
        if len(values) != n:
            raise ValueError("incorrect number of elements")
        self.type = type
        self.n = n
        self.values = values
        if self.type == 'r':
            self.shape = (1, n)
        elif self.type == 'c':
            self.shape = (n, 1)
        else:
            raise ValueError(f"found {type}, expected 'r' or 'c'")

    def __mul__(self, k):
        return Vector(self.type, self.n, [i * k for i in self.values])

    def __sub__(self, other):
        return Vector(self.type, self.n, [self.values[i] - other.values[i] for i in range(self.n)])

    def __str__(self):
        if self.type == 'r':
            s = '[ ' + " ".join(list(map(str, self.values))) + ' ]'
            return s
        else:
            ret = "["
            for i in range(self.n):
                if i == 0:
                    ret += f"{self.values[i]}\n"
                else:
                    ret += f" {self.values[i]}\n"
            ret = ret.strip('\n')
            ret += ']'
            return ret

def findSol(A, b):
    # creates the augmented matrix
    A.addCol(b)
    
    # performs elimination
    A.elimination()
    x = [None] * (A.m - 1)
    for i in range(A.n - 1, -1, -1):
        row = A.getRow(i).values
        k = 0
        j = -1
        while x[j]:
            k += x[j] * row[j - 1]
            j -= 1
        x[j] = (row[-1] - k) / row[j - 1]
    return Vector('c', len(x), x)

# solves Ax = b
A = Matrix((3, 3), [6, -1, 3, 5, 5, -5, 3, -1, 4])
print(A)
b = Vector('c', 3, [-9, 20, -5])
print(b)
x = findSol(A, b)
print(x)