###---Imports---###

import pylab as p

####-----------Functions-----------#######

####--------Basic_Math_Functions--------####

pi = 3.141592653589793

def Even(n):
    N = [x for x in xrange(1,n+1) if x%2==0]
    return N 

def Odd(n):
    N = [x for x in xrange(1,n+1) if x%2!=0]
    return N

def nrt(x,n=2.):
    A = (x)**(1/float(n))
    if A == int(A):
        A = int(A)
    return A
    
def SigSum(a,b,f):
    Sigma = 0
    for i in xrange(a,b+1):
        Sigma += f(i)
    return Sigma
    
def Facto(n):
    Ftl = 1
    for i in xrange(n,0,-1):
          Ftl *= i
    return Ftl

def Bpi(f,b,a=0):
    Product = 1
    for i in xrange(a,b+1):
        Product *= f(i)
    return Product

def Binomi(n,k):
    Binomial = (Facto(n)/1.)/((Facto(k))*(Facto(n-k)))
    return Binomial
    
def e(x=1,p=15):
    Sum = 1
    for i in xrange(1,p+1):
        Sum += 1./Facto(i)
    Sum = Sum**x
    return Sum 

def pif(k):
    Fsum = ((Facto(4*k))/((Facto(k))**(4)))*(((26390.*k)+(1103.))/((396.)**(4*k)))
    return Fsum

def pii(a=0,b=21):
    A = ((2.)/((99.)**(2)))*(nrt(2))
    B = SigSum(a,b,pif)
    ProductSum = ((A*B)**(-1))
    return ProductSum
 
def Sign(n):
    if n >= 0:
        return True
    else:
        return False
        
        
###___________Trigonometry__________###


def SinE(x,p=20):
    Sine = 0 
    for i in xrange(p+1):
        Sine += (((-1)**i)*((x/1.)**(1+2*i)))/(Facto(1+2*i))
    return Sine

def SinEA(K,p=20):
    Sine = []
    for k in K:
        S = SinE(k,p)
        Sine.append(S)
    return Sine

    
def CoSinE(x,p=20):
    CoSine = 0 
    for i in xrange(p+1):
        CoSine += (((-1)**i)*((x/1.)**(2*i)))/(Facto(2*i))
    return CoSine 

def CoSinEA(K,p=20):
    Cosine = []
    for k in K:
        C = CoSinE(k,p)
        Cosine.append(C)
    return Cosine


def SiN(x,p=20):
    if (type(x)==int) or (type(x)==float):
        Sine = SinE(x,p)
        return Sine
    else:
        Sine = SinEA(x,p)
        return Sine

def CoS(x,p=20):
    if (type(x)==int) or (type(x)==float):
        Cosine = CoSinE(x,p)
        return Cosine
    else:
        Cosine = CoSinEA(x,p)
        return Cosine    

def posCscE(x,p=20):
    cosecant = 1./SinE(x,p) 
    return cosecant

def negCscE(x,p=20):
    x = x*(-1)
    cosecant = (1./SinE(x,p))*(-1)
    return cosecant

def CscE(x,p=20):
    if x > 0:
        cosecant = posCscE(x,p)
        return cosecant
    else:
        cosecant = negCscE(x,p)
        return cosecant
        
def CscEA(X,p=20):
    cosecant = []
    for x in X:
        Csc = CscE(x,p)
        cosecant.append(Csc)
    return cosecant

def CsC(x,p=20):
    if (type(x)==int) or (type(x)==float):
        cosecant = CscE(x,p)
        return cosecant
    else:
        cosecant = CscEA(x,p)
        return cosecant
 
def SeC(x,p=20):
    Sum = 1./CoS(x,p)
    return Sum

def TaN(x,p=20):
    tan = (SiN(x,p))/(CoS(x,p))
    return tan
    
def CoT(x,p=20):
    Cotan = (CoS(x,p))/(SiN(x,p))
    return Cotan

##_____________Inverse_Trig____________###

def ASinF(x,i=0):
    F = (Binomi(2*i,i))*(((x/1.)**(2.*i+1))/(((4.)**i)*(2.*i+1)))
    return F
    
def arcSiN(x,p=85):
    if (x > 1) or (x <-1):
        Sum = 0
    elif (x == 1) or (x ==-1):
        if (x == 1):
            Sum = 0.060867358358260404
            for i in xrange(0,p+1):
                Sum += ASinF(x,i)
        else:
            Sum = -0.060867358358260404
            for i in xrange(0,p+1):
                Sum += ASinF(x,i)
    elif ((x < 1) and (x > .9)) or ((x > -1) and (x < -.9)):
        Sum = (x*(.005))*(0.060867358358260404)
        for i in xrange(0,p+1):
                Sum += ASinF(x,i)
    else:
        Sum = 0.0
        for i in xrange(0,p+1):
            Sum += ASinF(x,i)
    return Sum

def AsinL(K,p=85):
    Asin = [] 
    for k in K:
        A = arcSiN(k,p)
        Asin.append(A)
    return Asin

def AcosL(K,p=85):
    Acos = []
    for k in K:
        A = (pi/2.)-(arcSiN(k,p))
        Acos.append(A)
    return Acos
    
def posArct(x,p=80):
    if x > 1:
        Sum = 0.0
        for i in xrange(0,p+1):
            Sum += ((((-1.)**i)*(1./x)**(2.*i+1))/(2.*i+1))
        Sum = (pi/2.)-Sum
    if x == 1:
        Sum = -8.789025894295754e-05
        for i in xrange(0,p+1):
            Sum += ((((-1.)**i)*(x)**(2.*i+1))/(2.*i+1))
    elif x < 1:
        Sum = 0.0
        for i in xrange(0,p+1):
            Sum += ((((-1.)**i)*(x)**(2.*i+1))/(2.*i+1))
    return Sum

def negArct(x,p=80):
    X = (x)*(-1)
    Sum = posArct(X,p)*(-1)
    return Sum
    
def arcTaN(x,p=80):
    if x >= 0:
        arctan = posArct(x,p)
        return arctan
    else:
        arctan = negArct(x,p)
        return arctan

def AtanL(X,p=80):
    Atan = []
    for x in X:
        A = arcTaN(x,p)
        Atan.append(A)
    return Atan

def arcoTaN(x,p=80):
    arcotan = (pi/2.)-arcTaN(x,p)
    return arcotan
    
def ActanL(X,p=80):
    Actan = []
    for x in X:
        A = arcoTaN(x,p)
        Actan.append(A)
    return Actan

def ArcSin(X,p=85):
    if (type(X)==int) or (type(X)==float):
        arcsin = arcSiN(X,p)
        return arcsin
    else:
        arcsin = AsinL(X,p)
        return arcsin

def ArcCos(X,p=85):
    if (type(X)==int) or (type(X)==float):
        arccos = (pi/2.)-arcSiN(X,p)
        return arccos
    else:
        arccos = AcosL(X,p)
        return arccos

def ArcTan(X,p=80):
    if (type(X)==int) or (type(X)==float):
        arctan = arcTaN(X,p)
        return arctan
    else:
        arctan = AtanL(X,p)
        return arctan

def ArcoTan(x,p=80):
    if (type(x)==int) or (type(x)==float):
        arcotan = arcoTaN(x,p)
        return arcotan
    else:
        arcotan = ActanL(x,p) 
        return arcotan
    
def ArcSec(x,p=85):
    Sum = ArcCos(1./x,p)
    return Sum

def ArcCsc(x,p=85):
    Sum = ArcSin(1./x,p)
    return Sum
    
##---------------Hyperbolics--------------##  
   
def SinH(x):
    Val = ((1-e(-2*x))/(2*e(-x)))
    return Val
    
def CosH(x):
    Val = ((1+e(-2*x))/(2*e(-x)))
    return Val 

def TanH(x):
    Val = SinH(x)/CosH(x)
    return Val

def CotH(x):
    Val = CosH(x)/SinH(x)
    return Val

def SecH(x):
    Val = 1./CosH(x)
    return Val
    
def CscH(x):
    Val = 1./SinH(x)
    return Val 

def D2R2D(O,C=0):
    Deg = O*(180./pi)
    Rad = O*(pi/180.)
    if C==0:
        if Deg == int(Deg):
            Deg = int(Deg)
            return Deg
        else:
            return Deg
    else:
        if Rad == int(Rad):
            Rad = int(Rad)
            return Rad
        else:
            return Rad

def dx(f,a):
    h = 1e-5
    A = (f(a+h)-f(a))/(h)
    if A == (int(A) + (.0001)) or (int(A) - (.0001)):
        A = int(A)
    return A

def VAddSub(A,B,p=0):
    C = range(len(A))
    if len(A) == len(B):
        if p==0:
            for i in xrange(len(A)):
                C[i] = A[i]+B[i]
        else:
            for i in xrange(len(A)):
                C[i] = A[i]-B[i]
    return C

def VectLen(A):
    L = 0
    for a in A:
        L += a**2
    L = nrt(L)
    return L  

def VectDot(A,B):
    C = 0
    if len(A) == len(B):
        for i in xrange(len(A)):
            C += A[i]*B[i]
    return C

def CoSDot(A,B,O,p=1):
    AL = VectLen(A)
    BL = VectLen(B)
    O  = D2R2D(O,p)
    DotA = AL*BL*CoSinE(O)
    return DotA

def VectUnt(A):
    C = range(len(A))
    AL= VectLen(A)
    for i in xrange(len(A)):
        C[i] = A[i]/AL
    return C
    
def CrossPr(A,B):
    C = range(3)
    if len(A)==3 and len(B)==3:
        C[0] = (A[1]*B[2])-(A[2]*B[1])
        C[1] = (A[2]*B[0])-(A[0]*B[2])
        C[2] = (A[0]*B[1])-(A[1]*B[0])
    elif len(A)!=len(B):
        print "Only input '3-D' and '2-D' column vectors. Ex: For '2-D' vector, input [x,y,0],[x,0,z], or [0,y,z]."
    return C

##-----Useful-In-for-loop---functions----####

###---Error-function---###

def Eerror(A,B):
    '''Plot-point Error value comparson
    A=(Exact Value)
    B=(Value for comparison)
    Calculates with this form--> F = (A-B)/A'''
    F = (A-B)/A
    return F
    
###--Step-function--###

def Stp(A,B,C):
    '''Simple exp with scalar
    useful with for loop for stepping up incremental changes in delt dT
    A=(scaler)
    B=(base)
    C=(exponent)
    Calculates with this form--> F=A*B**(C)'''
    F=A*B**(C)
    return F
    
####________Simple_Functions_______####


def Fibonacci(N):
    A = [1,2]
    for i in xrange(2,N):
        A.append(A[i-1]+A[i-2])
        if A[i] > N:
            break
    return A
    
def Summer(A):
    total = 0
    for i in xrange(len(A)):
        total += A[i]
    return total
    
def FactorFlip(N,a):
    ''' This function is merely an embeded 
    list of length 2 that easily integrates 
    a math exspression into the Prime numbers 
    list builder.
    A = [6*N+(-1),6*N+(1)]'''
    A = [6*N+(-1),6*N+(1)]
    return A[a] 

def PrimeChk(N):
    ''' This function checks whether 
    a number is prime or composite by
    verifing the "len" of a list.
    Prime = True
    Composite = False
    N = Input number'''
    A = [1]
    for i in xrange(2,N+1):
        if (N%i == 0) and len(A) < 3:
            A.append(i)
    if len(A) == 3:
        A = False
    else:
        A = True
    return A

def Divsors(N):
    '''This function outputs an array 
    of factors up to the input value.'''
    A = [x for x in xrange(1,N+1) if N%x==0]
    return A

def PrimeFactors(N):
    '''This function returns a list of
    prime factors up to the input value.
    The output will be all prime numbers!!
    This function uses the "Divsors" and
    "PrimeChk" functions.  
    '''
    A = [x for x in xrange(1,N+1) if N%x==0 and PrimeChk(x)==True]
    return A

def Primality(N):
    ''' This fuction uses FactorFlip along with 
    PrimeChk to build a list of prime numbers.
    N = The "Aim value" number.
    Ex: if N = 709 then this fuction will 
    build an array that includes prime numbers less or equal to 709.
    list includes [1,2,3] to start.'''
    A = [1,2,3]
    num = 1
    while num < N:
        for j in xrange(0,2):
            P = FactorFlip(num,j)
            if PrimeChk(P)==True:
                if P > N:
                    break
                else:
                    A.append(P)
        num += 1
    return A

def PrimalistS(N):
    '''This fuction uses FactorFlip along with 
    PrimeChk to build a list of prime numbers.
    N = The "Nth" prime number.
    Ex: if N = 1000 then this function will 
    build an array up to the 1000th prime number.
    list includes [1,2,3] to start'''
    A = [1,2,3]
    num = 1
    K = 3
    while K < (N+1):
        for j in xrange(0,2):
            P = FactorFlip(num,j)
            if PrimeChk(P)==True:
                if K > N:
                    break
                else:
                    A.append(P)
        num += 1
        K    = len(A)
    return A
      
##-------Triangular-root-test------###

def Tri(x):
    '''This function REQUIRES numpy import.
    Triangular root test. 
    If N is an interger 
    then the input value is triangular
    x = (is the triangular value to be tested)
    Calculates with this form--> N = n.sqrt((8*x)+1) - 1/(2)'''
    N = (nrt((8.*x)+1) - 1.)/(2.)
    return N

##-------Polynomial-function------###

def f_x(A,B,C,D,x):
    '''This is a simple, single variable 
    polynomial function for 3rd degree polynomials.
    A=(coefficient of the 3rd power variable)
    B=(coefficient of the 2nd power variable)
    C=(coefficient of the 1st power variable)
    D=(Constant)
    x=(Independent Variable of the function)
    Calculates with this form--> f = (A*x**3)+(B*x**2)+(C*x)+(D)'''
    f = (A*x**3)+(B*x**2)+(C*x)+(D)
    return f
    
##-----Reimann-Sum-Integration-functions-----###

def LineCharge(L = 1.,H = 1.,lam0 = 1.):
    """Line charge function. numpy needed!!
       lambda0 = 1 by default
       H = distance of point away from line of charge
       L = lenght of line charge.
     Calculates in this form--->(1/k*lam0)((H*n.cos((n.pi*H)/(2.*L)))/(((2.*L)**2+(L)**2)**(3./2.)))"""
    eps0 = 8.854187817*10**(-12)
    k = 1./4.*pi*eps0
    E = (H*CoS((pi*H)/(2.*L)))/(((2.*L)**2+(L)**2)**(3./2.))
    return E
          
def ReimannR(function,a,b,N): 
    """ Right-hand Reimann Sum Function. 
        Imbeded with no function!! Requires a defined function as Input!!!
        function = function to be evaluated.Exclude args. when entering function. Ex: If h(x) = (e^x3)/2x + x^7;then input should be 'h'
        a = starting point
        b = ending point
        N = partition number **!!! Must be an int value
        This function uses a loop to calculate the area under a curve
        Calculates in this form----> F = function(i*h)*h """
    h = (b-a)/(N)
    AreaR = 0
    for i in xrange(0,N):
        AreaR += function(i*h)*h
    return AreaR
    
def ReimannL(function,a,b,N): 
    """ Left-hand Reimann Sum Function. 
        Imbeded with no function!! Requires a defined function as Input!!!
        function = function to be evaluated.Exclude args. when entering function. Ex: If g(x) = e^x3 + sin(x) + x^3;then input should be 'g'
        a = starting point
        b = ending point
        N = partition number **!!! Must be an int value
        This function uses a loop to calculate the area under a curve
        Calculates in this form----> F = function(i*h)*h """
    h = (b-a)/(N)
    AreaL = 0
    for i in xrange(0,N+1):
        AreaL += function(i*h)*h
    return AreaL
    
def ReimannM(function,a,b,N): 
    """ Mid-point Reimann Sum Function. 
        Imbeded with no function!! Requires a defined function as Input!!!
        function = function to be evaluated.Exclude args. when entering function. Ex: If f(x) = e^x + cos(x) + x^x;then input should be 'f'
        a = starting point
        b = ending point
        N = partition number **!!! Must be an int value
        This function uses a loop to calculate the area under a curve
        Calculates in this form----> F = 4*function(i*h)*h """
    h = .5*((b-a)/(N))
    AreaM = 0
    for i in xrange(0,N):
        AreaM += 4*function(i*h)*h
    return AreaM
    
def TrapezoidSum(function,a,b,N):
    """ Trapezoid method for the Reimann Sum of a Function. 
        Imbeded with no function!! Requires a defined function as Input!!!
        function = function to be evaluated.Exclude args. when entering function. Ex: If h(x) = tan(x) + (1/x)^x;then input should be 'h'
        a = starting point
        b = ending point
        N = partition number **!!! Must be an int value
        AR = Reimann Right-hand 
        AL = Reimann Left--hand 
        This function uses a loop to calculate the area under a curve
        Calculates in this form---->  TrapArea = .5*(AR+AL)"""
    AR = ReimannR(function,a,b,N)
    AL = ReimannL(function,a,b,N)    
    TrapArea = .5*(AR+AL)
    return TrapArea

   
def SimpsonSum(function,a,b,N):
    """ Simpson's method for the Reimann Sum of a Function. 
        Imbeded with no function!! Requires a defined function as Input!!!
        function = function to be evaluated.Exclude args. when entering function. Ex: If f(x) = e^(kx) + (ky/x);then input should be 'f'
        a = starting point
        b = ending point
        N = partition number **!!! Must be an int value
        This function uses a loop to calculate the area under a curve
        Calculates in this form---->  Ssum = (h/6.)*(fab+STrap+SMid)"""
    STrap = 0
    SMid  = 0
    Ssum  = 0
    h = (b-a)/(N)
    fab = function(a)+function(b) 
    for i in xrange(0,N+1):
        STrap += 2*(function(a+i*h))
    for k in xrange(0,N):
        SMid  += 4*(function(a+(k+.5)*h))
    Ssum = (h/6.)*(fab+STrap+SMid)
    return Ssum

#### -----------List/Matrix Stuff ------------ ####

def ListSum(A):
    """ 'ListSum' returns the Sum of elements within a list 
        A = sum list
        Calculates in this form---> Asum += float(a) """
    Asum = 0
    for a in A:
        Asum += float(a)
    return Asum
    
def ListAvg(A):
    """'Uses my 'ListSum' function to produce an avgerage value for the list
        A = List of numbers. Non-nested list.
        Calculates in this form---> SumAvg = float(ListSum(A))/float(len(A))"""
    SumAvg = float(ListSum(A))/float(len(A))
    return SumAvg
    
def ListAvg2DDis(A):
    """ Uses my 'ListSum' and 'ListAvg' functions to produce an avgerage value for the list
        A = List of numbers. Nested list. Made for my paired list functions.
        Uses embeded list format from my other functions. Ex: A == (A[0][0],A[0][1])
        Xavg = ListAvg(A[0])
        Yavg = ListAvg(A[1]) 
        Calculates in this form---> r = n.sqrt((Xavg)**2+(Yavg)**2)"""
    Xavg = ListAvg(A[0])
    Yavg = ListAvg(A[1])
    r = nrt((Xavg)**2+(Yavg)**2)
    return r
    
def ListsqSum(A):
    """ 'ListsqSum' returns the Sum of each element squared within a list 
        A = sum list
        Calculates in this form---> AsqSum += float(a**2) """
    AsqSum = 0
    for a in A:
        AsqSum += float(a**2)
    return AsqSum
    
def TwoListProdSum(A,B):
    """ Two List Product Sum 
     'TwoListProdSum' returns the Sum of two elements within a list
     List 'A' and 'B' must be the same size!! 
        A = sum list 1
        B = sum list 2
        Calculates in this form---> ABSum += float(A[i])*float(B[i]) """
    ABSum = 0
    for i in xrange(0,len(A)):
        ABSum += (A[i]/1.)*(B[i]/1.)
    return ABSum 
    
def MatrixMaker(N,rl,k):
    """ MatrixMaker is a function that... makes a "matrix"(list really) of with random integer numbers as there entries.
        N  = Number of list within
        rl = Number of elements within each list
        k  = variance number """
    A = range(N)
    for i in xrange(0,N):
        A[i] = []
        for j in xrange(0,rl):
            K = p.randint(k)
            A[i].append(K)
    return A

def MatrixMakeFlt(N,rl,k):
    """ MatrixMaker is a function that... makes a "matrix"(list really) of with random float numbers as there entries.
        N  = Number of list within
        rl = Number of elements within each list
        k  = variance number """
    A = range(N)
    for i in xrange(0,N):
        A[i] = []
        for j in xrange(0,rl):
            K = p.randint(k)*p.random()
            A[i].append(K)
    return A

####_____________Best Fit Solutions_____________######

def BestFitL(A,B):
    """ Best Fit line gen.
        Only useful 2-D matrix
        Uses 'linalg.solve(A_,r_)'
        A & B must be list!!
        Cannot directly use the return value:
            return is inbeded as such: matrix([[ some num],[ sum number2]])"""
    Suma = ListSum(A)
    Sumb = ListSum(B)
    Ssqb = ListsqSum(B)
    ABSum= TwoListProdSum(A,B)
    N    = len(A)
    A_= n.matrix([[N,Sumb],[Sumb,Ssqb]])
    r_= n.matrix([[Suma],[ABSum]])
    x = n.linalg.solve(A_,r_)
    return x   

####______Signal Processing______####

def signal2Psd(y,a,b,k):
    t = n.arange(a,b,k)
    N = len(t)
    y= y(t)
    Y = n.fft.fft(y)
    indf = n.arange(1,(N/2)+1) 
    f = n.fft.fftfreq(N)
    Psd = (abs(Y[indf]))**2+(abs(Y[-indf]))**2
    return "PSD =",Psd,"&",f
       
#### ------------Template for a Function and it Derivative---------    
    
def newt(f,fd,t):
    """ Newton's Method for root finding
        f  = function
        fd = Derivative of the 'f' function
        t  = Independent variable of the 'f' and fd function 
        calculates in this form --> x = t + ((f)/(fd))"""
    x = t + ((f)/(fd))
    return x 
  

##--------##--Moddeling functions--##--------##

def SleepWalkerSAvgDis(A):
    """ Uses my ListAvg2DDis function to produce an avgerage value for multiple list
        A = List of numbers. Uses Nested list-pair. Made for my paired list functions.
        Uses embeded list format from my other functions. Ex: A == (A[0][0],A[0][1]) == (Ax,Ay)
        Calculates in this form---> Avg[i] = ListAvg2DDis(A[i])"""
    Avg = range(len(A))
    for i in xrange(len(A)):
        Avg[i] = ListAvg2DDis(A[i])
    return Avg
    
    
###_________Random_Walker_Functions__________###


def SleepWalker(T,V):
    """This function produces an array simulation of a drunk walker, 1-D.
        T = time interval
        V = variance number""" 
    r = p.random(T)
    for i in xrange(1,T):
        r[0] = 0
        if r[i] >= .5:
            r[i] = r[i-1] + p.randint(V)*p.random()

        else:
            r[i] = r[i-1] - p.randint(V)*p.random()
    return r
    
def SleepWalker2D(T,V):
    """This function produces an array simulation of a drunk walker in 2-D.
        T = time interval
        V = variance number"""
    x = SleepWalker(T,V)
    y = SleepWalker(T,V)
    return x,y

def SleepWalkerS(T,N,V):
    """This function produces an array simulation of multiple drunk walkers.
        Output is in embeded list format A == A[0],A[1],A[2],A[3]
        T = time interval
        N = Number of sleepwalkers
        V = variance number"""
    A = range(N)
    for i in xrange(0,N):
        A[i] = SleepWalker(T,V)
    return A
    
def SleepWalkerS2D(T,N,V):
    """This function produces an array simulation of multiple drunk walkers in 2-D.
        Output is in embeded list format A == (A[0][0],A[0][1]) == (A[x],A[y]) The x and y of the index '[0]'
        T = time interval
        N = Number of sleepwalkers
        V = variance number"""
    A = range(N)
    for i in xrange(0,N):
        A[i] = SleepWalker2D(T,V)
    return A 
###______Random_Walker_Square_movement_Function__###//
### square steps ###
def SleepWalkersq(T):
    """This function produces an array simulation of a drunk walker, 1-D.
        T = time interval
        V = variance number""" 
    r = p.random(T)
    for i in xrange(1,T):
        r[0] = 0
        if r[i] >= .5:
            r[i] = r[i-1] + 1
        else:
            r[i] = r[i-1] - 1
    return r
    
def SleepWalkersq2D(T):
    """This function produces an array simulation of a drunk walker in 2-D.
        T = time interval
        V = variance number"""
    x = SleepWalkersq(T)
    y = SleepWalkersq(T)
    return x,y

def SleepWalkerSsq(T,N):
    """This function produces an array simulation of multiple drunk walkers.
        Output is in embeded list format A == A[0],A[1],A[2],A[3]
        T = time interval
        N = Number of sleepwalkers
        V = variance number"""
    A = range(N)
    for i in xrange(0,N):
        A[i] = SleepWalkersq(T)
    return A
    
def SleepWalkerSsq2D(T,N):
    """This function produces an array simulation of multiple drunk walkers in 2-D.
        Output is in embeded list format A == (A[0][0],A[0][1]) == (A[x],A[y]) The x and y of the index '[0]'
        T = time interval
        N = Number of sleepwalkers
        V = variance number"""
    A = range(N)
    for i in xrange(0,N):
        A[i] = SleepWalkersq2D(T)
    return A   

###------Mach-Calculation------###

def Mach(A,B):
    '''Mach function. Useful for single-use or iterative calculating for-loops.
    A = (speed of object)
    B = (speed of sound at the location of object)
    Calculates with this form--> M = A/B'''
    M = A/B
    return M  

###------Freefall-velocity-calculation-with-drag------###

def freefallV(M,A,t,g,P):
    '''Exact freefall Velocity function 
    that accounts for the drag Force
    Solves for velocity
    M = (mass of object)
    g = (acceleration due to gravity)
    P = rho = (density of the fluid)
    A = (cross sectional area of object)
    t = t[time value]=(time)
    Calculates with this form--> (n.sqrt((2*M*g)/(p*A)))*(n.tanh(t*n.sqrt((g*p*A)/(2*M))))'''
    V = (nrt((2*M*g)/(P*A)))*(n.tanh(t*nrt((g*P*A)/(2*M))))
    return V
    

###---Exact--decay--function---###
    
def dekay(A,B,C):
    '''Exact decay function
    numpy import needed!!
    A=StartingNucleiValue=(Scalar)
    B=-(DeltaTime)=-(Numerator of exponent)
    C=Tau=(Denomnator of exponent)
    Calculates with this form--> E = A*(n.exp(-B/C))'''
    E = A*(n.exp(-B/C))
    return E
    
##--##---Methods of approxmation---##--##
        
##---Runga Kutta----###

def RungKdekay(A1,A0,B,C):
    '''Runga Kutta method for decay
    A1= N[t] = N[current step]
    A0= N[mid step] = N[t] - N[t](delt/2*tau)
    B =(step size)
    C =tau
    Returns A2=N[next step]
    Define A0 as N[t]-N[t]*(delt/2*tau)
    Calculates with this form-->  A2 = A1 - A0*(B/C)'''
    A2 = A1 - A0*(B/C)
    return A2  
    
###--Euler's-decay-function---###

def Eulerdekay(A,B,C):
    '''Useful for only step.Euler method for nuKe decay
    A=(intial/current value)
    B=-(Delta time)=Numerator
    C=(tau)=Denomnator
    Calculates with this form--> E = A + A*(-B/C)'''
    E = A + A*(-B/C)
    return E

###--Euler's-freefall-function---###

def EulerFreefV(M,A,delt,g,P,V0):
    '''M = (mass of object)
    g = (acceleration due to gravity)
    P = rho = (density of the fluid)
    A = (cross sectional area of object)
    delt = t[time interval] = (delta t)
    V0= V[current step]  
    Calculates with this form--> V = V0 + g*t - ((P*A*(V0**2)*t)/(2*M))'''
    V = V0 + (g*delt) - ((P*A*(V0**2)*delt)/(2*M))
    return V
    
###----Deg2rads converter-----##

def D2Rads(O):
    '''Converts degrees to rads'''
    D = O*(3.14159265359)/(180.)
    return D
    
    
#### Bisection method of root finding #####
""" Detail more detail!!"""
def Bisect(f,a,b,toler):
    if p.sign(f(a)) == p.sign(f(b)):       ###Initial If-Statement###
        print "The signs of f(a) and f(b) are the same. Your two starting 'f(t)' values don't satisfy the intermediate value theorem! Pick new values" 
    else:                                  ### else-statement ####
        while (abs(a-b)) >= toler:         ### while function #####
            c=(a+b)/2.                     ### bisection center point calculation ####
            if p.sign(f(c)) == p.sign(f(a)): ### If-statement for independent value switching ####
                print c                    ### realtime readout for 'c' value #####
                a = c                      ### reasignment of 'a' value ####
            else:
                b = c                      ### reasignment of 'b' value ####
                print c                    ### realtime readout for 'c' value ####
    return c