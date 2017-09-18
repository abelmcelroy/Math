import random
import math
import itertools
import numpy
import CryptoSolver as CS
import fractions
class Ten_Twentyfour(object):
        def __init__(self,array):
                self.A = array
        def Print(self):
                print self.A[0]
                print self.A[1]
                print self.A[2]
                print self.A[3]
                print "____________"
        def swipe(self,direction):
                if direction[0] == "L":
                        self.moveL()
                        self.moveL()
                        self.moveL()
                        self.pairL()
                        self.moveL()
                        self.moveL()
                if direction[0] == "R":
                        self.moveR()
                        self.moveR()
                        self.moveR()
                        self.pairR()
                        self.moveR()
                        self.moveR()
                if direction[0] == "U":
                        self.moveU()
                        self.moveU()
                        self.moveU()
                        self.pairU()
                        self.moveU()
                        self.moveU()
                if direction[0] == "D":
                        self.moveD()
                        self.moveD()
                        self.moveD()
                        self.pairD()
                        self.moveD()
                        self.moveD()
                if direction[0] == "P":
                        self.changeVal(int(direction[1]),int(direction[2]),int(direction[3]))
        def moveR(self):
                for x in range(4):
                        for y in range(2,-1,-1):
                                if self.A[x][y+1] == 0:
                                        self.A[x][y+1] = self.A[x][y]
                                        self.A[x][y] = 0
        def moveL(self):
                for x in range(4):
                        for y in range(3):
                                if self.A[x][y] == 0:
                                        self.A[x][y] = self.A[x][y+1]
                                        self.A[x][y+1] = 0
        def moveU(self):
                for y in range(4):
                        for x in range(1,4):
                                if self.A[x-1][y] == 0:
                                        self.A[x-1][y] = self.A[x][y]
                                        self.A[x][y] = 0
        def moveD(self):
                for y in range(4):
                        for x in range(2,-1,-1):
                                if self.A[x+1][y] == 0:
                                        self.A[x+1][y] = self.A[x][y]
                                        self.A[x][y] = 0
        def pairU(self):
                for x in range(3):
                        for y in range(4):
                                if self.A[x][y] == self.A[x+1][y]:
                                        self.A[x][y] = 2*self.A[x][y]
                                        self.A[x+1][y] = 0
        def pairD(self):
                for x in range(2,-1,-1):
                        for y in range(4):
                                if self.A[x][y] == self.A[x+1][y]:
                                        self.A[x+1][y] = 2*self.A[x][y]
                                        self.A[x][y] = 0
        def pairL(self):
                for y in range(3):
                        for x in range(4):
                                if self.A[x][y] == self.A[x][y+1]:
                                        self.A[x][y] = 2*self.A[x][y]
                                        self.A[x][y+1] = 0
        def pairR(self):
                for y in range(2,-1,-1):
                        for x in range(4):
                                if self.A[x][y] == self.A[x][y+1]:
                                        self.A[x][y+1] = 2*self.A[x][y]
                                        self.A[x][y] = 0
        def changeVal(self,a,b,v):
                self.A[a-1][b-1] = v
        def Sum(self):
                return sum(self.A[0])+sum(self.A[1])+sum(self.A[2])+sum(self.A[3])
def isPalendrome(S):
        front = SS[:int(math.floor(len(S)/2.0))]
        back = SS[int(math.ceil(len(S)/2.0)):]
        front = reverse(front)
        return front == back
def reverse(S):
        T = ""
        for i in range(len(S)-1,-1,-1):
                T+=S[i]
        return T
def findLongestPal(S):
        SindexOfLongest = -1
        EindexOfLongest = -1
        Longest = ""
        for i in range(len(S)):
                for j in range(i,len(S)+1):
                        SS = S[i:j]
                        front = SS[:int(math.floor(len(SS)/2.0))]
                        back = SS[int(math.ceil(len(SS)/2.0)):]
                        front = reverse(front)
                        if front == back:
                                SindexOfLongest = i
                                EindexOfLongest = j
                                if abs(i-j) > len(Longest):
                                        Longest = SS
        return Longest
def removePal(S):
        return S.split(findLongestPal(S))
        ### doesnt do it the most efficient way yet
def parsePalendromes(T):
        remaining = [T]
        sol = []
        while len(remaining)>0:
                temp = []
                for string in remaining:
                        temp.append(removePal(string))
                        sol.append(findLongestPal(string))
                remaining = []
                for S in temp:
                        for s in S:
                                if len(s)>0:
                                        remaining.append(s)
        return sol
def MinPathSnails(F,target):
        A = []
        for a in range(len(F)):
                A.append([])
                for b in range(len(F[0])):
                        A[a].append(10**30)             
        A[len(A)-1][0] = F[len(A)-1][0]
        for y in range(0,len(F)):
                for x in range(len(F)-1,-1,-1):
                        if ((x == len(F)-1) and (y == 0)):
                                A[x][y] = int(F[x][y])
                        elif x == len(F)-1:
                                A[x][y] = int(F[x][y])+int(A[x][y-1])
                        elif y == 0:
                                A[x][y] = int(F[x][y])+int(A[x+1][y])
                        else:
                              A[x][y] = min(int(F[x][y])+int(A[x+1][y]),int(F[x][y])+int(A[x][y-1]))
        return A[target[0]][target[1]]
def MinPathFrank(F):
        A = []
        for a in range(len(F)):
                A.append([])
                for b in range(len(F[0])):
                        A[a].append(10**30)
        for x in range(len(F)-1,-1,-1):
                for y in range(0,len(F)):
                        if x == len(F)-1:
                                A[x][y] = int(F[x][y]) + y+1
                        elif y == 0:
                                A[x][y] = int(F[x][y]) + min(int(A[x+1][y])+2,int(A[x+1][y+1])+1)
                        elif y == len(F[x])-1:
                                A[x][y] = int(F[x][y]) + min(int(A[x+1][y])+2,int(A[x+1][y-1])+3)
                        else:
                                A[x][y] = int(F[x][y]) + min(int(A[x+1][y])+2,int(A[x+1][y+1])+1,int(A[x+1][y-1])+3)
        return min(A[0])
def Fnm(n,m):
        if fractions.gcd(n,m) == 1:
                return "1"
        else:
                return "0"
def shortestCycle(code):
        cycle = code[0]+code[1]
        count = 2
        while code != cycle*(len(code)/len(cycle)):
                cycle+=code[count]
                count+=1
        return cycle
def powerSum(L,power):
        sol = 0
        for l in L:
                sol+=l**power
        return sol
def quad_formula(a,b,c):
    if b**2-4*a*c >= 0:
        return [(-b+math.sqrt(b**2-4*a*c))/2*a,(-b-math.sqrt(b**2-4*a*c))/2*a]
def max_regions_of_inscribed_ngon(n):
    return 1+nCk(n,2)+nCk(n,4)
# f(n+1) = f(n-1) + nCk(n,3) + (n-1)
def max_regions_of_ngon(n):
    return int((n*(n-3)*(n*(n-3)+14)/24.0)+1)
def Expected_value_coupon_problem_one_of_each(Ps):
        expectedValue = 0
        probabilities = Ps
        A = 1
        while A <= len(probabilities):
                for combo in itertools.combinations(probabilities, A):
                        denominator = 0
                        for Pi in combo:
                                denominator += Pi
                        expectedValue+= (-1)**(A+1)*1.0/denominator
                A+=1
        return expectedValue
def permutations(code):
    if len(code)==1:
        return [code]
    perms=permutations(code[1:])
    char=code[0]
    result=[]
    for perm in perms:
        for i in range(len(perm)+1):
            result.append(perm[:i] + char + perm[i:])
    return result
def permu(code):
        if len(code) == 2:
                return [code, code[1] + code[0]]
        codes = []
        new_code = code
        while True:
                codes += map(lambda p: new_code[0] + p, permu(new_code[1:]))
                new_code = new_code[1:] + new_code[0]
                if new_code == code:
                        break
        return codes
def Intersections_of_cords_with_N_Points_on_Circumfrance(n):
        count = 0
        if n > 3:
                count+= Intersections_of_cords_with_N_Points_on_Circumfrance(n-1)
        for i in range(1,n-2):
                count+=(n-i-2)*i
        return count
def isPrime(n):
    if n == 1:
        return False
    f = sieve(n/2)
    for p in f:
        if n%p==0:
            return False
    return True
def nCk(n,k):
    if n<k:
        return 0
    else:
        return long(math.factorial(n)/(math.factorial(n-k)*math.factorial(k)))
def sieve(limit):
        nPrime = set()
        primes = []
        for i in range(2,limit+1):
                if i in nPrime:
                        continue
                for j in range(i**2, limit+1, i):
                        nPrime.add(j)
                primes.append(i)
        return primes
def good_numbers_less_than_N(n):
        primes = sieve(n+100)
        goods = []
        for x in xrange(10,n+1):
                count = 0
                for y in xrange(int(x/math.log(x))+5):
                        for z in xrange(y,int(x/math.log(x))+5):
                                if primes[y] + primes[z] == x:
                                        count+=1
                                if count>2:
                                        break
                        if count>2:
                                break
                if count == 2:
                        goods.append(x)
        return goods
def distance(pointOne,pointTwo):
        return math.sqrt((pointTwo[0]-pointOne[0])**2+(pointTwo[1]-pointOne[1])**2+(pointTwo[2]-pointOne[2])**2)
def fibonacci(index):
        five = math.sqrt(5)
        return (((1+five)/2)**index-((1-five)/2)**index)/five
def Stirling_first(n,k):
        A = [[0 for y in range(101)] for x in range(101)]
        A[0][0]=1
        for i in xrange(1,101):
                A[i][0] = 0
                A[i][i] = 1
                for j in xrange(1,101):
                    A[i][j] = (((i-1)*(A[i-1][j])+A[i-1][j-1]))
        return A[n][k]
def Stirling_second(n,k):
    num = 0
    for j in range(k+1):
        num+=(-1)**(k-j)*(nCk(k,j))*j**n
    return num/math.factorial(k)
class TiledRectangle(object):
        def __init__(self, size1, size2):
                self.width = size1
                self.height = size2
                self.dimentions = [self.width, self.height]
        def WaysToTile(self):
                total = 1
                for i in xrange(1,self.height/2+1):
                        for j in xrange(1,(self.width/2)+1):
                                total = total*4*(math.cos((math.pi*i)/(self.height+1))**2+math.cos((math.pi*j)/(self.width+1))**2)
                return total
def count_domino_tilings(m,n):
    result = 1
    for k in range(1,m+1):
        for j in range(1,n+1):
            result = result*math.sqrt(4*math.cos((k*math.pi/(m+1)))**2+4*math.cos((j*math.pi/(n+1)))**2)
    return math.sqrt(result)
class RationalSpiral(object):
        def __init__(self, index):
                self.index = index
                self.orientation = "R"
                self.coordinates = self.findCoordinates()
        def findCoordinates(self):
                coordinates = [0,0]
                movesLeft = self.index
                while movesLeft > 0:
                        if abs(coordinates[0]) == abs(coordinates[1]) and coordinates[0] != 0:
                                if coordinates[0] < 0 or coordinates[1] > 0:
                                        self.turn()
                        if coordinates[0] > coordinates[1] and coordinates[0] == -coordinates[1]+1:
                                        self.turn()
                        if self.orientation == "R":
                                coordinates[0] += 1
                        if self.orientation == "U":
                                coordinates[1] += 1
                        if self.orientation == "L":
                                coordinates[0] -= 1
                        if self.orientation == "D":
                                coordinates[1] -= 1
                        movesLeft-=1
                return coordinates
        def turn(self):
                newOrientation = ""
                if self.orientation == "R":
                        newOrientation = "U"
                if self.orientation == "U":
                        newOrientation = "L"
                if self.orientation == "L":
                        newOrientation = "D"
                if self.orientation == "D":
                        newOrientation = "R"
                self.orientation = newOrientation
        def specs(self):
                return [self.index, self.coordinates, self.coordinates[0]+self.coordinates[1], self.orientation]
def Find_arithmatic_progressions_of_N_in_array(array,N):
        progressions = []
        for num1 in array:
                i = 1
                progression = [num1]
                while num1+(i*N) in array:
                        progression.append(num1+(i*N))
                        i+=1
                if len(progression)>1:
                        progressions.append(len(progression))
        L = progressions.sort()
        if len(progressions)>0:
                return progressions[len(progressions)-1]
        else:
                return 0
def tetrateAndMod(number):
        result = 2
        for x in xrange(3,number+1):
                result = result**x
        result = result%number+1
        return result
def countTrailing0s(num):
        count = 0
        for x in xrange(1,num+1):
                y = x
                subcount = 0
                while y%5==0:
                        subcount+=1
                        y=y/5
                count += subcount
        return count
def countTrailing0sFactorial(num):
        count = 0
        for x in xrange(1,num+1):
                count += countTrailing0s(x)  
        return count
def convertBase(number,base):
    highestExp = 0
    newRepresentation = ""
    while base**highestExp<=number:
        highestExp+=1
    highestExp-=1
    while highestExp>=0:
        nextDigit = "0"
        for x in range(base-1,0,-1):
            if x*(base**highestExp)<=number:
                if x < 10:
                    nextDigit = str(x)
                elif x == 10:
                    nextDigit = "A"
                elif x == 11:
                    nextDigit = "B"
                elif x == 12:
                    nextDigit = "C"
                elif x == 13:
                    nextDigit = "D"
                elif x == 14:
                    nextDigit = "E"
                elif x == 15:
                    nextDigit = "F"
                elif x == 16:
                    nextDigit = "G"
                elif x == 17:
                    nextDigit = "H"
                elif x == 18:
                    nextDigit = "I"
                elif x == 19:
                    nextDigit = "J"
                elif x == 20:
                    nextDigit = "K"
                elif x == 21:
                    nextDigit = "L"
                elif x == 22:
                    nextDigit = "M"
                elif x == 23:
                    nextDigit = "N"
                elif x == 24:
                    nextDigit = "O"
                elif x == 25:
                    nextDigit = "P"
                elif x == 26:
                    nextDigit = "Q"
                elif x == 27:
                    nextDigit = "R"
                elif x == 28:
                    nextDigit = "S"
                elif x == 29:
                    nextDigit = "T"
                elif x == 30:
                    nextDigit = "U"
                elif x == 31:
                    nextDigit = "V"
                elif x == 32:
                    nextDigit = "W"
                elif x == 33:
                    nextDigit = "X"
                elif x == 34:
                    nextDigit = "Y"
                elif x == 35:
                    nextDigit = "Z"
                number-=x*(base**highestExp)
        newRepresentation+=nextDigit
        highestExp-=1
    return newRepresentation
def convertBinary(number):
    highestExp = 0
    newRepresentation = ""
    while 2**highestExp<=number:
        highestExp+=1
    highestExp -=1
    while highestExp>=0:
        nextDigit = "0"
        if (2**highestExp)<=number:
            nextDigit = "1"
            number-=(2**highestExp)
        newRepresentation+=nextDigit
        highestExp-=1
    return newRepresentation
def XOR_binary(O,T):
    One = convertBinary(O)
    Two = convertBinary(T)
    result = ""
    if len(One) > len(Two):
        anex = "0"*(len(One)-len(Two))
        Two = anex+Two
    else:
        anex = "0"*(len(Two)-len(One))
        One = anex+One
    for i in range(len(One)):
        if One[i] == Two[i]:
            result+="0"
        else:
            result+="1"
    return result
def count_triangles_given_N_segments(n):
#   might need work
    i = 3
    while nCk(i,2) <= n:
        i+=1
    i-=1
    return nCk(i,3) 
def coupon_collection_probability_ON(n,d):
    return ((1.0*math.factorial(d))/d**n)*Stirling_second(n-1,d-1)
def coupon_collection_probability_loosely_before(n,d):
    total = 0.0
    for x in range(d,n+1):
        total+=coupon_collection_probability_ON(x,d)
    return total
def harmonic_number(i):
    total = 0.0
    for x in range(1,i+1):
        total+=1.0/x
    return total
def uniform_coupon_EV(n):
    return harmonic_number(n)*n
def coupon_expected_value_single_set(couponPs):
    total = 0.0
    parody = 1
    for x in range(1,len(couponPs)+1):
        for combination in itertools.combinations(couponPs,x):
            subtotal = 0.0
            for p in combination:
                subtotal+=p
            total+=1.0/(subtotal*parody)
        parody=parody*(-1)
    return total
def Euler_phi(n):
    phi = n
    primeContenders = sieve(n)
    primeFactors =[]
    for p in primeContenders:
        if n%p==0:
            primeFactors.append(p)
    for P in primeFactors:
        phi = phi*(1-(1.0/P))
    return int(phi)
def points_per_side_one_divided_diamond(n):
    return 1 + 2**(n-1)
def triangles_in_divided_diamond(n):
    triangles = 4
    points = pointsPerSide(n)
    triangles += 4*PPL.nCk(points,2)
    return triangles