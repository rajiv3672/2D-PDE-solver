'''
@@@@@@@@@@@@@@@@@
Rajiv Das
2015-132-036
@@@@@@@@@@@@@@@@@
===================================================================
YOU MUST READ EVERY COMMENT FOR UNDERSTANDING THE ALGORITHM.IT WILL
BE HELPFUL FOR YOURSELF IF YOU KEEP THE DIFFERENCE EQUATION OF THE
POISSON EQUATION IN YOUR HAND.

====================================
STEADY STATE GROUND WATER DYNAMICS
IS GOVERNED BY THE POISSON EQUATION.
====================================

===================================================================
'''


# This import registers the 3D projection, but is otherwise unused.
#from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import

import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np


'''
====================================
Boundary condition u_i,0 and u_i,Ny
If Domain,D = [a,b]X[c,d] then the
values of u(x,y), for any value of x
at y=c and y=d is given by this 
condition.
====================================
'''
def boundX(i,j,c,d):
    U_x0=1                            #u(x,y) at y=c:::1st boundary
    U_xNy=1                           #u(x,y) at y=d:::2nd boundary
    if j==c:
        return U_x0
    elif j==d:
        return U_xNy
    else:
        return 0

'''
====================================
Boundary condition u_0,j and u_Nx,j
If Domain,D = [a,b]X[c,d] then the
values of u(x,y), for any value of y
at x=a and x=b is given by this 
condition.
====================================
'''
def boundY(i,j,a,b):
    U_0y=0                              #u(x,y) at x=a:::1st boundary
    U_Nxy=0                           #u(x,y) at x=b:::2nd boundary
    if i==a:
        return U_0y
    elif i==b:
        return U_Nxy
    else:
        return 0

'''
======================================
This is our function for source term.
This can be any function of x and y.
In language of ground water dynamics 
this is called inflitaration.
======================================
'''
def Source(i,j):
    #f can be edited to any finction of x and y
    f=np.sin(i*j)*np.cos(j*i)
    return f


#Taking input of boundaries
a=0
b=1
c=0
d=1

#Here we fix our number of solution points:
Nx=50
Ny=50

#Calculating our step size:
hx=(b-a)/Nx
hy=(d-c)/Ny

#Coefficents of difference equation/sparse matrix
kx=(1/hx)
ky=(1/hy)
kxy=2*(kx+ky)





'''
=================================================
Suppose i=1,j=1 putting this on the difference
equation we see there are u_0,1 and u_1,0 which
are known from boundary conditions. These terms
has to be moved on the right side of the equation.
In the following loop the R.H.S of the  equation 
was constructed as discussed above.
=================================================
'''
R_ij=[]
for i in np.arange(a+hx,b,hx):
    for j in np.arange(c + hy, d, hy):

        '''
        ================================
        This section calculates the 
        boundary terms.
        ================================
        '''
        X_0=kx*boundY(round(i-hx,4),j,a,b)         #Zero if i-hx is not equal to "a"
        Y_0=ky*boundX(i,round(j-hy,4),c,d)         #Zero if j-hy is not equal to "c"
        X_Nx=kx*boundY(round(i+hx,4),j,a,b)        #Zero if i+hx is not equal to "b"
        Y_Ny = ky * boundX(i, round(j + hy,4),c,d) #Zero if j+hy is not equal to "d"
        B=X_0+X_Nx+Y_0+Y_Ny
        R_ij.append((Source(i,j))-B)

#Creating meshgrid for plotting the solution.
x=[]
y=[]
for i in np.arange(a+hx,b,hx):
    x.append(i)
for j in np.arange(c+hy,d,hy):
    y.append(j)

x,y=np.meshgrid(x,y)

N=(Nx-1)*(Ny-1)
a = []
for i in range(Nx-1):
    for j in range(Ny-1):
        for k in range(Nx-1):
            for l in range(Ny-1):
                if k == i and l == j:
                    a.append(-kxy)
                elif k == i and l == j - 1:
                    a.append(ky)
                elif k == i and l == j + 1:
                    a.append(ky)
                elif k == i - 1 and l == j:
                    a.append(kx)
                elif k == i + 1 and l == j:
                    a.append(kx)
                else:
                    a.append(0)


Cof=np.reshape(a,(N,N))
SOL=np.linalg.solve(Cof,R_ij)
z=np.reshape(SOL,(Nx-1,Ny-1))
Z=np.transpose(z)


fig=plt.figure()
ax=fig.add_subplot(111,projection='3d')
ax.plot_wireframe(x,y,Z,color='black')
ax.set_xlabel('X axis')
ax.set_ylabel('Y axis')
ax.set_zlabel('Potential')
ax.set_title('Dirichlet 2D Poisson equation')

plt.show()
