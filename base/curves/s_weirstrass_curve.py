from gmpy2 import mpz, legendre, powmod, add, invert
from math import sqrt
from sympy import mod_inverse
import numpy as np
from datetime import datetime

# 
# tonelli_shanks() :- implementation of Tonelli-Shanks algorithm
# 
# @n : the quadratic residue for which the solution 
# is to be found
# 
# @p : the prime number with respect to which the
# solution has to be found
# 
# finds x : (x*x) ≡ n mod p 
# where p = 4k + 3
# 

def tonelli_shanks(n, p):
  x = p-1
  s = 0
  while(x%2!=1):
    x=x//2
    s+=1
  q = (p-1)//(2**s)
  z = 2
  while(isResidue(z,p)):
    z+=1
  M = s
  c = pow(z,q,p)
  t = pow(n,q,p)
  R = pow(n,(q+1)/2,p)
  while(True):
    if(t==0):
      return 0
    if(t==1):
      return R
    i = 0
    while(pow(t,2**i,p)!=1):
      i+=1
    b = pow(c,2**(M-i-1),p)
    M = i
    c = pow(b,2,p)
    t = (t*c)%p
    R = (R*b)%p

# 
# generatePoints() :- generates points according to the 
# given curve parameters and stores them in file
# 
# @a, @d : parameters of the Short Weiestrass Curve
# 
# @p : the prime number with respect to which the
# points have to be found
# 
# @fname : takes the name of the file to be created
# to store the points (default : points.txt)
# 
# The general equation for Short Weiestrass Curve is :
# 
# y^2 = x^3+ax+d
# 
# 
# 
# opens a new file with given name
# 
# for every x from 0 to p, checks if the corresponding
# value of y^2 is a quadratic residue
# 
# if it is, uses appropriate function to find y
# depending on nature of p
# 
# writes the points into the opened file
# 

def generatePoints(a, d, p, start=0):
  #Creating the output file
  x_array = []
  y_array = []

  # $$$
  if start > p:
    start = 0

  # take at max 1000 points
  for x in range(start, min(start+1000, p)):
      fx = ((x*x*x)+(a*x)+d)%p   #finding values of y^2 mod p for every integer value of x in range
      if(isResidue(fx, p) or fx == 0):
        # 4k+3 form
        if((p-3)%4 == 0):
          y = pow(fx, (p+1)/4, p)   #euler's method
        # 4k+1 form
        else:
          y = tonelli_shanks(fx,p)   #Tonneli-Shank's method

        x_array.append(x)
        y_array.append(y)
        if fx != 0:
          x_array.append(x)
          y_array.append(p-y)
  return (x_array,y_array)
#


def find_points(a, b, p):
    no_of_points = 0
    for x in range(p):
        temp1 = x**3
        temp2 = x * a
        temp3 = temp2 + b
        xx = temp3 + temp1
        xx %= p

        if xx == 0:
            y = xx
            no_of_points += 1

        ret = legendre(xx, p)
        if ret == 1:
            temp1 = (p + 1) // 4
            y = pow(xx, temp1, p)
            no_of_points += 1

            y = p - y
            no_of_points += 1

    no_of_points += 1  # Including the point at infinity.
    return no_of_points
#
# pow() :- helper function to execute
# fast exponentiation with modular operation
# 
# @a, @b : base and index values respectively for
# exponentiation operation
# 
# @m : the value to be used for modular part of
# the operation
# 
# returns a^b mod m
# time complexity : O(logn)
#
 
def pow(a, b, m):   
  if(b == 0):
    return 1%m
  ans = pow(a, b//2, m) 
  if(b%2 == 0):
    return (ans*ans)%m
  else:
    return ((ans*ans)%m*a)%m

# 
# isResidue() :- helper function to check whether a 
# given number is a quadratic residue with respect to 
# a given prime number
# 
# @x : the number to check for quadratic residue
# 
# @p : the prime number with respect to which the
# condition has to be checked
# 
# Euler's criterion for QR:
# 
# if x^((p-1)/2) ≡ 1 mod p, then x is QR
# if x^((p-1)/2) ≡ -1 mod p, then x is QNR
# 

def isResidue(x, p):
  return pow(x,(p-1)/2,p) == 1

# 
# addpoints() :- function to perform addition operation 
# on two points in the curve using the affine addition formula
# 
# @p1, @p2 : two input points to perform addition operation
# 
# formula for affine addition is as follows:-
# 
# gradient = (y2 - y1) / (x2 - x1)
# x3 = (gradient^2) - x2 - x1
# y3 = (gradient*(x2 - x3)) - y2
# 
# this formula is altered for finite field and used here
# final answer is returned as a tuple of (x3,y3)
# 

def addpoints(a,d,p,p1,p2):
  # print(p1)
  # print(p2)   

  try:
    # pm = p - p2[1]
    # if p1[0] == 0 and p1[1] == 0:
    #   return p2
    # if p2[0] == 0 and p2[1] == 0:
    #   return p1
    # if p1[1] == pm:
    #   return (0,0)
    gradient = (p2[1]-p1[1])*mod_inverse((p2[0]-p1[0]),p)
    x = (gradient**2-p2[0]-p1[0])%p
    y = (gradient*(p1[0]-x)-p1[1])%p
  except:
    print("-*-*-*-*-*-*-*Inverse does not exist, please change the point!!*-*-*-*-*-*-*-")
    return (0,-1)
  return (x,y)

# 
# subtractpoints() :- function to subtract one given point from another
# 
# @p1, @p2 : two inputs to perform subtraction operation (p1-p2)
# 
# it can be deduced from original curve equation that the 
# additive inverse of (x,y) is (x,-y)
# 
# ∴ (x1,y1) - (x2,y2) = (x1,y1) + (-x2,y2)
# 
# here we are directly using this fromula and a previously defined 
# function addpoints to perfrom the subtraction operation
# 

def substractpoints(a,d,p,p1,p2):
  p3 = (p2[0],-1*p2[1])
  return addpoints(a,d,p,p1,p3)

# 
# doublepoint() :- function to perform doubling operation 
# on a given point in the curve using the affine addition formula
# 
# @p1 : two input points to perform addition operation
# 
# substituting the value of x2 and y2 with x1 and y1
# respectively and substituting values from original 
# equation of the curve gives us the following formula
# 
# lamda = ((3*x1*x1) + a) / (2*y1)
# x = (lamda^2) - (2*x1)
# y = (lamda*(x1-x)) - y1
# 
# this formula is altered for finite field and used here
# final answer is returned as a tuple of (x,y)
# 

def doublepoint(a,d,p,p1):
  try:                      
    lam = (3*p1[0]*p1[0]+a)*mod_inverse((2*p1[1]),p)                      
    x = (lam**2-2*p1[0])%p
    y = ((lam*(p1[0]-x))-p1[1])%p
  except:
    print("-*-*-*-*-*-*-*Inverse does not exist, please change the point!!*-*-*-*-*-*-*-")
    return (0,-1)
  return (x,y)

# 
# multiplypoint() :- function to perform scalar multiplication
# with a given a point and scalar value
# 
# @p1, @scalar : the input point and scalar value to perform
# saclar multiplication of point
# 
# multiplication is perfromed with repeated addition operations 
# with some optimization using doubling and addition method
# time complexity : O(logn)
# 

def multiplypoint(a,d,p,p1, scalar):
  pt = (0,1)

  if scalar == 1:
    pt = p1
  elif scalar%2 == 1:
    pt = addpoints(a,d,p,p1,multiplypoint(a,d,p,p1,scalar-1))
  else:
    pt = multiplypoint(a,d,p,doublepoint(a,d,p,p1),scalar//2)
  return pt


def bsgs(a,d,p,p1,p2):
  p1,p2 = p2, p1
  x_arr, y_arr = generatePoints(a,d,p)
  # n = len(x_arr) + 1
  n = find_points(a,d,p)
  m = sqrt(n)+1
  m = int(m)
  count = m
  mp = multiplypoint(a, d, p, p1, m)
  # ip = [[mpz(0) for _ in range(3)] for _ in range(m+1)]
  # jmp = [[mpz(0) for _ in range(3)] for _ in range(m+1)]
  ip = [[0 for _ in range(3)] for _ in range(count+1)]
  jmp = [[0 for _ in range(3)] for _ in range(count+1)]

  for i in range(1,count+1):
        # t3 = mpz(i)
        pd = multiplypoint(a, d, p, p1, i)
        ip[i][0], ip[i][1], ip[i][2] = i, pd[0], pd[1]

  for i in range(1,count+1):
        # t5 = mpz(i)
        pd = multiplypoint(a, d, p, mp, i)
        jmp[i][0], jmp[i][1], jmp[i][2] = i, pd[0], pd[1]

  for i in range(1,count+1):
        if jmp[i][2] != 0: 
            jmp[i][2] = p - jmp[i][2]

  for i in range(1,count+1):
        # t5 = jmp[i][1]
        # t6 = jmp[i][2]
        ptmp = (jmp[i][1], jmp[i][2])
        ps = addpoints(a, d, p, p2, ptmp)
        jmp[i][1] = ps[0]
        jmp[i][2] = ps[1]

  match, I, J = 0, 0, 0
  for i in range(1,count+1):
            # J = 0
            # I += 1
            if match == 1:
                break
            for j in range(1,count+1):
                J += 1
                if ip[i][1] == jmp[j][1]:
                    if ip[i][2] == jmp[j][2]:
                        J = j
                        match = 1
                        break
            I += 1
  
  if match == 0:
     return -1
  else:
     k = int(((m*J) + I) % n)
    #  k = int(n)
    #  return (k, 0)
     return k


