#file: find_points.py
#montgomery curve : By^2 = x^3 + Ax^2 + x

from sympy import mod_inverse, nextprime
import random
# import graph_points as graph

#hasse's theorem
def hassesTheorem(prime):
    upperBound = int(prime + 1 + 2*(prime ** 0.5))
    lowerBound = int(prime + 1 - 2*(prime ** 0.5))
    print("According to hasse's theorem the total number of points should be in the range of {} {}".format(lowerBound, upperBound))

#gets us the next prime if number isn't prime
def getPrime(number):
    new_prime = nextprime(number - 1)
    if number == new_prime:
        return number
    else:
        print("{} isn't a prime so we consider the next prime {} for calculations".format(number, new_prime))
        return new_prime

#fast modular exponentiation
def modular_pow(base, exponent, modulus):
    result = 1
    while exponent > 0:
        if exponent & 1:
            result = (result * base) % modulus
        exponent = exponent >> 1
        base = (base * base) % modulus
    return result

#euler formula to calculate y
def euler(quadraticResidue, prime):
    return modular_pow(quadraticResidue, (prime+1) // 4, prime)

#y^2 = m (mod p); here m is the function of x
def findM(a, b, x, p):
    return ((x*x*x + a*x*x + x) * (mod_inverse(b, p))) % p

def findQRForSW(a, b, x, p):
    return (x*x*x + a*x + b) % p

def legendre(a, p):
    return modular_pow(a, (p - 1) // 2, p)


def power_modulo(a: int, b: int, n: int) -> int:
    """ Computes a ^ b mod n """
    result = 1
    while b != 0:
        if b % 2 == 1:
            # b odd
            result = (result * a) % n
        a = (a * a) % n
        b >>= 1
    return result


def extended_gcd(a: int, b: int) -> (int, int, int):
    # optional check
    if a == 0:
        return b, 0, 1

    # without this check the first iteration will divide by zero
    if b == 0:
        return a, 1, 0

    un_prev = 1
    vn_prev = 0
    un_cur = 0
    vn_cur = 1

    while True:
        qn = a // b
        new_r = a % b
        a = b
        b = new_r

        if b == 0:
            return a, un_cur, vn_cur

        # Update coefficients
        un_new = un_prev - qn * un_cur
        vn_new = vn_prev - qn * vn_cur

        # Shift coefficients
        un_prev = un_cur
        vn_prev = vn_cur
        un_cur = un_new
        vn_cur = vn_new


def inverse_modulo(a: int, n: int) -> int:
    _, b, _ = extended_gcd(a, n)
    return b % n


def legendre_symbol(a: int, p: int, /) -> int:
    return power_modulo(a, (p - 1) >> 1, p)


def _choose_b(p: int, /) -> int:
    b = 2
    while legendre_symbol(b, p) == 1:
        b = random.randrange(2, p)
    return b


def _tonelli_shanks_recursive(a: int, k: int, p: int, b: int, b_inverse: int, /):
    """
    Computes a square root of a modulo prime p
    :param a: the number to take the square root of
    :param k: positive integer, such that a^m = 1 (mod p) where m = (p-1)/(2^k)
    :param p: odd prime p modulo which we are working
    :param b: an arbitrary non-square modulo p
    :param b_inverse: the inverse of b modulo p, i.e., b * b_inverse = 1 (mod p)
    :return: one of the square roots of a modulo p (the other can be obtained via negation modulo p)
    """

    m = (p - 1) >> k
    a_m = 1

    while m % 2 == 0 and a_m == 1:
        m >>= 1
        k += 1
        a_m = power_modulo(a, m, p)

    if a_m == p - 1:
        # a^m = -1 (mod p)
        b_power = 1 << (k - 1)
        b_power_half = 1 << (k - 2)
        a_next = (a * power_modulo(b, b_power, p)) % p
        a_next_root = _tonelli_shanks_recursive(a_next, k, p, b, b_inverse)
        a_root = a_next_root * power_modulo(b_inverse, b_power_half, p)
        return a_root % p

    # we now handle the case when m is odd
    # this case is easy, a^((m+1)/2) is a square root of a
    return power_modulo(a, (m + 1) >> 1, p)


def tonelli_shanks(a: int, p: int, /):
    """
    Computes a square root of a modulo prime p
    :param a: the number to take the square root of
    :param p: odd prime p modulo which we are working
    :return: one of the square roots of a modulo p (the other can be obtained via negation modulo p)
    """

    if legendre_symbol(a, p) != 1:
      # a is not not a square modulo p
      return None

    b = _choose_b(p)
    b_inverse = inverse_modulo(b, p)
    return _tonelli_shanks_recursive(a, 1, p, b, b_inverse)

def generatePoints(a, b, p, start = 0):
    x_coordinates = []
    y_coordinates = []

    if start > p:
      start = 0

    totalPoints = 0
    p = getPrime(p)
    print("Prime field is:", p)
    print("")
    for i in range(start, min(start+1000, p)):
        m = findM(a, b, i, p)
        
        quadraticResidue = legendre(m, p)
        if quadraticResidue == 1:
            if p % 4 == 3:			
                y = euler(m, p)
            else:			
                y = tonelli_shanks(m, p)

            print("({}, {})  ({}, {})".format(i, y, i, p-y))
            x_coordinates.extend([int(i), int(i)])
            y_coordinates.extend([int(y), int(p-y)])

            totalPoints += 2
        elif quadraticResidue == 0:
            print("({}, {})".format(i, 0))
            x_coordinates.append(i)
            y_coordinates.append(0)
            totalPoints += 1

    print("\nNumber of total points (including 0 points) is:", totalPoints)
    print("-----------------------------------------------------------------------------------------")
    return (x_coordinates, y_coordinates)


#add_points.py
def addpoints(a, b, p, p1, p2):
    x1, y1 = p1
    x2, y2 = p2
    # x1 = p1[0]
    # y1 = p1[1]
    # x2 = p2[0]
    # y2 = p2[1]

    if x1 == x2 and y1 == y2 :
        print('Both the points are same! Perform point doubling operation instead addition.')
        return doublepoint(a, b, p, p1)
    elif x1 == x2 :
        print("Addition of two points ({}, {}) and ({}, {}) is = ({}, {})".format(x1, y1, x2, y2, 0, 0))
        return (0, 0)
    else :
        try :
            k = (((y2 - y1) % p) * mod_inverse((x2 - x1) % p, p)) % p
            x3 = (((((((b * ((k ** 2) % p)) % p) - a) % p) - x1) % p) - x2) % p
            y3 = ((((((((((2 * x1) % p) + x2) % p) + a) % p) * k) % p) - (b * ((k ** 3) % p)) % p) - y1) % p
            
            if x3 < 0:
                x3 = x3 + p
            if y3 < 0:
                y3 = y3 + p
            
            print("Addition of two points ({}, {}) and ({}, {}) is = ({}, {})".format(x1, y1, x2, y2, x3, y3))
            return (x3, y3)
        except Exception as e: 
            print(e)
            return (0, -1)

#subtraction
def substractpoints(a, b, p, p1, p2):
    x1, y1 = p1
    x2, y2 = p2

    print("Subtracting two points ({}, {}) and ({}, {})".format(x1, y1, x2, y2))
    y2 = -y2
    if y2 < 0 :
        y2 = y2 + p

    p2 = x2, y2
    return addpoints(a, b, p, p1, p2)

#doubling
def doublepoint(a,b,p,p1):
    x, y = p1
    if y == 0:
      return (0, 0)
    try :
        k = (((((((((x ** 2) % p) * 3) % p) + (((2 * a) % p) * x) % p) % p) + 1) % p) * mod_inverse((((2 * b) % p) * y) % p, p)) % p
        x3 = (((((((b * ((k ** 2) % p)) % p) - a) % p) - x) % p) - x) % p
        y3 = ((((((((((2 * x) % p) + x) % p) + a) % p) * k) % p) - (b * ((k ** 3) % p)) % p) - y) % p

        if x3 < 0:
            x3 = x3 + p
        if y3 < 0:
            y3 = y3 + p
        
        print('Performing Point Doubling on the point ({}, {})'.format(x, y))
        print('Resultant Point : ({}, {})'.format(x3, y3))
        return (x3, y3)
    except Exception as e: 
        print(e)
        return (0, -1)

#scalar multiplication
def multiplypoint(a,b,p,p1,k):
    print("k = ", k)
    # x0 = 0
    # y0 = 0
    x, y = p1
    p0 = (0, 0)
    idx = k.bit_length()
    while idx >= 0:
        print("idx = ", idx)
        if k & (1 << idx):
            # x0, y0 = pointAddition(x0, y0, x1, y1, a, b, p)
            p0 = addpoints(a, b, p, p0, p1)
            # x1, y1 = pointDoubling(x1, y1, a, b, p)
            p1 = doublepoint(a, b, p, p1)
            
        else:
            # x1, y1 = pointAddition(x0, y0, x1, y1, a, b, p)
            p1 = addpoints(a, b, p, p0, p1)
            # x0, y0 = pointDoubling(x0, y0, a, b, p)
            p0 = doublepoint(a, b, p, p0)
        idx = idx - 1
    
    print("({}, {}) multiplied by {} is = ({}, {})".format(x, y, k, p0[0], p0[1]))
    return p0