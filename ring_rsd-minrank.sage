  #####################################################################
  #  Solving Systems of Algebraic Equations Over
  #   Finite Commutative Rings and Applications  
  #####################################################################

           # H. Tchatchiem Kamche (hermann.tchatchiem@gmail.com) and
           # H. Tale Kalachi  (herve.tale@univ-yaounde1.cm)


####################################################
# Example on computing Grobner basis using SageMath#
####################################################

Z8xy.<x,y>=PolynomialRing(Integers(8),2,order='lex')
L1=[4*x^2*y +y^3 + 2*y+4,  4*x*y^2]   
L2=ideal(L1)
L3=L2.groebner_basis()
"""
print(L3)
[4*x^2*y + y^3 + 2*y + 4, 4*x*y^2, y^4 + 2*y^2 + 4*y, 2*y^3 + 4*y]
"""

P0=[(x^2-x)^2-2*(x^2-x), (y^2-y)^2-2*(y^2-y)]
P1=[4*x^2*y +y^3 + 2*y+4,  4*x*y^2]+P0  
P2=ideal(P1)
P3=P2.groebner_basis()
"""
print(P3)
x^4 - 2*x^3 - x^2 + 2*x, y^2 + 4, 2*y + 4]
"""

#######################################
# First Example on the MinRank Problem#
#######################################
#
# In this example we will use the Kipmis-Shamir Modeling and
# the Support-Minors Modeling to solve a MinRank Problem.
#
m=4
r=1
n=4
k=3
p=2
nu=3
s=binomial( n , r)
R =Integers(p^nu) 
RX.<z1,z2,z3,z4,x1,x2,x3>= PolynomialRing(R,k+s, order='lex')
X1=[x1,x2,x3]
X2=[z1,z2,z3,z4]

M1=matrix(R,
[[0,0,0,7],
[1,0,0,5],
[0,1,0,2],
[0,0,1,4]])

M2=matrix(R,
[[0,0,7,4],
[0,0,5,3],
[1,0,2,5],
[0,1,4,2]])

M3=matrix(R,
[[2,2,0,4],
[4,2,0,6],
[0,4,2,4],
[0,6,6,0]])

Mx=x1*M1+x2*M2+x3*M3
###
# Solving using the Kipmis-Shamir Modeling
##

Z=block_matrix([[identity_matrix(RX,n-r)],[matrix(RX,list(X2[0:n-r]))]])

"""
print(Z)
[ 1  0  0]
[ 0  1  0]
[ 0  0  1]
[--------]
[z1 z2 z3]
"""
K1=Mx*Z
K2=[RX(K1[i,j]) for i in range(m) for j in range(n-r)]
K3=RX.ideal(K2)
K4=K3.groebner_basis()

"""
print(K4)
[2*z1*x3 - 2*x3, 2*z2*x3 - 2*x3, 2*z3*x3 - 2*x3, x1 + 2*x3,
 x2 + 2*x3, 4*x3]
"""


###
# Solving using the Support-Minors Modeling
##

L1=[]
for i in range(4):
    L1=L1+[Mx[i,0]*z2-Mx[i,1]*z1,
           Mx[i,0]*z3-Mx[i,2]*z1,
           Mx[i,0]*z4-Mx[i,3]*z1,
           Mx[i,1]*z3-Mx[i,2]*z2,
           Mx[i,2]*z4-Mx[i,3]*z2,
           Mx[i,2]*z4-Mx[i,2]*z4]

Mon=[X1[i]*X2[j]  for j in range(len(X2)) for i in range(len(X1))]
"""
print(Mon)
[z1*x1, z1*x2, z1*x3, z2*x1, z2*x2, z2*x3, z3*x1, z3*x2, z3*x3,
z4*x1, z4*x2, z4*x3]
"""

L2= [[RX(poly).monomial_coefficient(c) for c in Mon] for poly in L1]
L3=matrix(L2)
L4=matrix(ZZ,L3).echelon_form()
L5=matrix(R,L4)

"""
print(L5[:len(Mon),:])
[1 0 0 0 0 0 0 0 0 0 0 2]
[0 1 0 0 0 0 0 0 0 0 0 2]
[0 0 2 0 0 0 0 0 0 0 0 2]
[0 0 0 1 0 0 0 0 0 0 0 2]
[0 0 0 0 1 0 0 0 0 0 0 2]
[0 0 0 0 0 2 0 0 0 0 0 2]
[0 0 0 0 0 0 1 0 0 0 0 2]
[0 0 0 0 0 0 0 1 0 0 0 2]
[0 0 0 0 0 0 0 0 2 0 0 2]
[0 0 0 0 0 0 0 0 0 1 0 2]
[0 0 0 0 0 0 0 0 0 0 1 2]
[0 0 0 0 0 0 0 0 0 0 0 4]
"""
#

########################################
# Second Example on the MinRank Problem#
########################################
#
# In this example we will see that adding some rings equations simplifies the resolution
#
p=2
m=3
nu=3
k=3
n=3
r=1
R=Integers(p^nu)
RX.<z1,z2,x1,x2,x3>= PolynomialRing(R,k+(n-r), order='lex')
X1=[x1,x2,x3]
X2=[z1,z2]

M0=matrix(R,
[[5,2,3],
[5,1,4],
[4,3,6]])

M1=matrix(R,
[[1,2,0],
[0,1,3],
[0,2,1]])

M2=matrix(R,
[[0,2,1],
[1,0,3],
[0,5,5]])

M3=matrix(R,
[[0,5,5],
[0,1,0],
[1,2,5]])

Mx=M0+x1*M1+x2*M2+x3*M3

###
# First case
###
# In this case, We use 'Z' in the classical form but
# we don't get a solution because '2' is in a Grobner basis.
#
Z=block_matrix([[identity_matrix(RX,n-r)],[matrix(RX,list(X2[0:n-r]))]])
"""
print(Z)
[ 1  0]
[ 0  1]
[-----]
[z1 z2]
"""

K1=Mx*Z
K2=[RX(K1[i,j]) for i in range(m) for j in range(n-r)]
K3=RX.ideal(K2)
K4=K3.groebner_basis()

"""
print(K4)
[z1*x1 - z1 + 3*x1 + x3 + 3, z1*x2 - 3*z1 - 3*x1 - x2 + 3*x3,
 z1*x3 - 3*x2 + x3 + 1, z2*x1 - z2 + 3*x2 + x3 + 1,
 z2*x2 - 3*z2 - x1 + x2, z2*x3 - x1 - 3*x2 + 3*x3,
 x1^2 - x2^2 - 3*x2*x3 + 3*x3^2 - x3,
 x1*x2 + 3*x1 + x2 - 3*x3^2 + 1,
 x1*x3 + 3*x2^2 - x3^2 + 3*x3 - 3,
 x2^3 - x2^2 - 3*x2*x3^2 + 3*x2 + x3^3 + x3^2 - 3,
 2]
"""

#
###
# Second case
###
# In this case, we choose a suitable permutation but
# we obtain a system which is not easy to solve.
#
Z=block_matrix([[matrix(RX,list(X2[0:n-r]))],[identity_matrix(RX,n-r)]])

"""
print(Z)
[z1 z2]
[-----]
[ 1  0]
[ 0  1]
"""
K1=Mx*Z
K2=[RX(K1[i,j]) for i in range(m) for j in range(n-r)]
K3=RX.ideal(K2)
K4=K3.groebner_basis()

"""
print(K4)
[z1*x1 - z1 - 3*x3 + 2,
 z1*x2 - 3*z1 + x1 + x3 + 1,
 z1*x3 + x2 - 3, 2*z1,
 z2*x1 - 3*z2 + x2 - 3*x3 + 3,
 z2*x2 - 3*z2 + 3*x1 + 3*x2 + 4,
 z2*x3 - 3*x1 + x2 + x3 - 2,
 2*z2 + 4,
 x1^2 - x2^2 - 3*x2*x3 + 3*x3^2 + x3 + 4,
 x1*x2 - 3*x1 + 3*x2 + x3^2 + 3,
 x1*x3 - x2^2 + x3^2 + 3*x3 - 3,
 2*x1 - 2,
 x2^3 - x2^2 + x2*x3^2 + 3*x2 + x3^3 - x3^2 - 3,
 2*x2 + 2,
 2*x3 + 4]
"""

#
###
# Third Case
###
# In this case, we add the rings equations '(x^2-x)^2-2*(x^2-x)' and
# we directly obtain the solution.

K1=Mx*Z
K2_1=[RX(K1[i,j]) for i in range(m) for j in range(n-r)]
K2_2=[(x^2-x)^2-2*(x^2-x) for x in X1+X2]
K2=K2_1+K2_2
K3=RX.ideal(K2)
K4=K3.groebner_basis()

"""
print(K4)
[z1^4 - z1^2, 2*z1, z2^4 + 3*z2^2 + 4, 2*z2 + 4, x1 - 1, x2 - 3, x3 + 2]
"""
#
########################################
# Example for the Rank Decoding Problem#
########################################

def HenselLiftOfPrimitivePolynomial(p,nu,m):
    """
    Input: 'p' the characteristic of the residue field,
    'm' the dimension of the Galois extension,
    'nu' the nilpotency index.
    Output: a monic polynomial 'h' in 'Z_ {p ^ nu}[z]' of degree 'm'
    such that 'h'  divises 'z^(p^m-1)-1' and
    the projection of 'h' in 'GF(p)[z]' is a primitive polynomial.
    """
    Zpz.<z>=QQ[]
    Hensel=Zpz(z^(p^m-1)-1).hensel_lift(p, nu)
    Conway=conway_polynomial(p,m)
    Fpz.<z>=GF(p)[]
    i=0
    while  Fpz(Conway)!=Fpz(Hensel[i]) :
        i=i+1
    return Hensel[i]

def MatrixRepresentationOf(v):
    """
    Input:'v' a list with coefficient in 'S=R[a]'
    Output: the matrix representation of 'v' in the
    ring 'R' relative to the basis '(1,a,...,a^(m-1))'
    """
    return matrix(R,len(v),m,[v[j].list()
                              for j in [0..len(v)-1]]).transpose()

def RankOf(A):
    # This function computes the rank of a matrix
    R=A.base_ring()
    L=matrix(ZZ,A)
    D=matrix(R,L.smith_form()[0])
    ar=0
    r=min(D.nrows(),D.ncols())
    while ar<r and R(D[ar,ar])!=R(0) :
        ar=ar+1
    return ar

def BasisRingRepresentationOf(P):
    # This function expands polynomial 'P' from 'S' to 'R'
    DicP=P.dict()
    ValP=list(DicP.values())
    KeyP=list(DicP.keys())
    l=len(ValP)
    ValPnew=[]
    for i in range(l):
        ValPnew=ValPnew+[ValP[i].list()]
    return [RX(dict((KeyP[k],ValPnew[k][i])
                    for k in range(l))) for i in range(m)]

def sigma(i):
    # This function is the power 'i'
    # of a generator of the Galois group
    return S.hom([a^(p^i)])


############
# Example I#
############

# The following example is given in our manuscript

p=2
m=3
nu=3
k=1
n=3
r=1

R=Integers(ZZ(p**nu))
Rz.<z>=R[]
Bpoly=HenselLiftOfPrimitivePolynomial(p,nu,m)
S.<a>=Rz.quotient(Bpoly)

"""
print(S)
Univariate Quotient Polynomial Ring in a over
Ring of integers modulo 8 with modulus z^3 + 6*z^2 + 5*z + 7
"""

g=matrix(S,[ 1,2*a^2 + a + 2,a^2 + 3*a])
y=matrix(S,[4*a^2 + 3*a + 3, 5*a^2 + 7*a + 6, 2*a^2 + 4*a + 5])

M0=MatrixRepresentationOf([-y[0,i] for i in range(3)])
M1=MatrixRepresentationOf([g[0,j] for j in range(3) ])
M2=MatrixRepresentationOf([a*g[0,j] for j in range(3) ])
M3=MatrixRepresentationOf([a^2*g[0,j] for j in range(3) ])

"""
print(M0)
[5 2 3]
[5 1 4]
[4 3 6]

print(M1)
[1 2 0]
[0 1 3]
[0 2 1]

print(M2)
[0 2 1]
[1 0 3]
[0 5 5]

print(M3)
[0 5 5]
[0 1 0]
[1 2 5]

"""

###
# Solving Using Support-Minors Modeling
###
SX=PolynomialRing(S,6,'z1,z2,z3,x0,x1,x2',order='lex')
[z1,z2,z3,x0,x1,x2]=SX.gens()
x=x0+x1*a+x2*a^2
[g1,g2,g3]=[g[0,j] for j in range(3)]
[y1,y2,y3]=[y[0,j] for j in range(3)]
T=[(x*g1-y1)*z2-(x*g2-y2)*z1,(x*g1-y1)*z3-(x*g3-y3)*z1,(x*g2-y2)*z3-(x*g3-y3)*z2]

T2=[(x*g1-y1)*z2-(x*g2-y2), (x*g1-y1)*z3-(x*g3-y3),(x*g2-y2)*z3-(x*g3-y3)*z2]
RX=SX.change_ring(R)
[z1,z2,z3,x0,x1,x2]=RX.gens()

T3=[]
for j in range(len(T2)):
    T3=T3+BasisRingRepresentationOf(T2[j])

X=[z2,z3,x0,x1,x2]
T4=T3+[(X[i]^2-X[i])^2-2*(X[i]^2-X[i]) for i in range(5)]
T5=RX.ideal(T4)
T6=T5.groebner_basis()
"""
print(T3)
[z2*x0 - 3*z2 - 2*x0 - 2*x1 + 3*x2 - 2,
 z2*x1 - 3*z2 - x0 - x2 - 1,
 z2*x2 + 4*z2 - 2*x0 + 3*x1 - 2*x2 - 3,
 z3*x0 - 3*z3 - x1 + 3*x2 - 3,
 z3*x1 - 3*z3 - 3*x0 - 3*x1 + 4,
 z3*x2 + 4*z3 - x0 + 3*x1 + 3*x2 + 2,
 -z2*x1 + 3*z2*x2 - 3*z2 + 2*z3*x0 + 2*z3*x1 - 3*z3*x2 + 2*z3,
 -3*z2*x0 - 3*z2*x1 + 4*z2 + z3*x0 + z3*x2 + z3,
 -z2*x0 + 3*z2*x1 + 3*z2*x2 + 2*z2 + 2*z3*x0 - 3*z3*x1 + 2*z3*x2 + 3*z3]

print(T6)
[z2^4 - z2^2, 2*z2, z3^4 + 3*z3^2 + 4, 2*z3 + 4, x0 - 1, x1 - 3, x2 + 2]
"""

###
# Solving by Linearization
###

A=matrix(S,[[-y[0,j],g[0,j],sigma(1)(g[0,j]),-sigma(1)(y[0,j])]
            for j in range(3)])

"""
print(A)
[4*a^2 + 5*a + 5               1               1   a^2 + 4*a + 5]
[  3*a^2 + a + 2   2*a^2 + a + 2 7*a^2 + 6*a + 6     6*a^2 + 5*a]
[6*a^2 + 4*a + 3       a^2 + 3*a 2*a^2 + 7*a + 2 6*a^2 + 2*a + 7]
"""
## We calculate the row echelon form using Magma
#
#

"""
p:=2;
nu:=3;
m:=3;
R:=IntegerRing(p^nu); //GaloisRing(p,nu,1);
Zz<z>:=PolynomialRing(IntegerRing(),1);
S<a> := GaloisRing(2,3, PolynomialRing(IntegerRing())!(z^3 + 6*z^2 + 5*z + 7));
A:=Matrix(S,[[4*a^2 + 5*a + 5,1,1,a^2 + 4*a + 5],
[3*a^2 + a + 2,2*a^2 + a + 2,7*a^2 + 6*a + 6, 6*a^2 + 5*a],
[6*a^2 + 4*a + 3,a^2 + 3*a,2*a^2 + 7*a + 2 ,6*a^2 + 2*a + 7]]);
EchelonForm(A);
[              1         a^2 + a               0       2*a^2 + 4]
[              0               2               0     6*a^2 + 4*a]
[              0               0               1 3*a^2 + 6*a + 3]
"""
#
#
"""
print(-sigma(2)(3*a^2 + 6*a + 3))
6*a^2 + 3*a + 1
"""

###
# Solving With Grobner Bases
###

SX=PolynomialRing(S,k*m+r*m,'x0,x1,x2,z0,z1,z2')
[x0,x1,x2,z0,z1,z2]=SX.gens()
X=[x0,x1,x2,z0,z1,z2]
W=[[X[i+m*j] for i in range(m)] for j in range(k)]
P=[[X[i+m*j] for i in range(m)] for j in [k..r+k-1]
  ]+[[1]+[0 for i in range(m-1)]]
L1=[sum([P[u][v]*W[i][l]*a^v*sigma(u)(a^l*g[i,j])
         for u in range(r+1) for v in range(m)
         for i in range(k) for l in range(m)])
    -sum([P[u][v]*a^v*sigma(u)(y[0,j])
          for u in range(r+1) for v in range(m)])
    for j in range(n)]
RX=SX.change_ring(R)
L2=[]
for j in range(n):
    L2=L2+BasisRingRepresentationOf(L1[j])
L3= RX.ideal(L2)
L4=L3.groebner_basis()

"""
print(L4)
[x0 - 1, x1 - 3, x2 + 2, 2*z0 + 2, 2*z1, 2*z2 + 2]
"""

###################################
# Example II
####

# In this example we assume that 'G' and 'e' are random.

p=2    # the characteristic of finite field
nu=3   # the nilpotency index
m=8    # the degree of Galois extension
k=2    # dimension of the linear code
n=8    # the length of the linear code
r=2    # the rank of error
R=Integers(p^nu)    # base ring
Rz.<z>=R[]
Conway=HenselLiftOfPrimitivePolynomial(p,nu,m) #Galois extension of 'R'
S.<a>=Rz.quotient(Conway)
SX=PolynomialRing(S,k*m+r*m,'x')
X=SX.gens()
W=[[X[i+m*j] for i in range(m)] for j in range(k)]
P=[[X[i+m*j] for i in range(m)] for j in [k..r+k-1]
  ]+[[1]+[0 for i in range(m-1)]]

# a generator matrix
G=random_matrix(S,k,n)

# a transmitted message
w=random_matrix(S,1,k)

# a random error of rank 'r'
D=diagonal_matrix(R,[R.random_element() for _ in range(r)])
E=random_matrix(R,m,r)*D*random_matrix(R,r,n)
while  RankOf(E)!=r:
    D=diagonal_matrix(R,[R.random_element() for _ in range(r)])
    E=random_matrix(R,m,r)*D*random_matrix(R,r,n)
e=matrix(S,[a^i for i in range(m)])*matrix(S,E)

# the received word
y=w*G+e

# Decoding using Gröbner bases

# 1) Algebraic modelling
L1=[sum([P[u][v]*W[i][l]*a^v*sigma(u)(a^l*G[i,j])
         for u in range(r+1) for v in range(m)
         for i in range(k) for l in range(m)])
    -sum([P[u][v]*a^v*sigma(u)(y[0,j])
          for u in range(r+1) for v in range(m)])
    for j in range(n)]

# 2) Expanding Equations from 'S' to 'R'
RX=SX.change_ring(R)
L2=[]
for j in range(n):
    L2=L2+BasisRingRepresentationOf(L1[j])
L3= RX.ideal(L2)

# 3) Computing a Gröbner basis
L4=L3.groebner_basis()

# 4) Computing the transmitted message
v=[S([-L4[i+j*m].constant_coefficient()
      for i in range(m)]) for j in range(k)]

# Check if 'v' is equal to 'w'
"""
print(v==list(w[0]))
True
"""

