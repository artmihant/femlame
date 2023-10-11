import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve


def D_matrix(Nu):
    return np.array([
        [1, Nu/(1-Nu), 0],
        [Nu/(1-Nu), 1, 0],
        [0, 0, (1-2*Nu)/(2-2*Nu)]
    ])*((1-Nu)/(1+Nu)/(1-2*Nu))


def B_matrix(a,b,c):
    return np.array([
        [b[1]-c[1], 0, c[1]-a[1], 0, a[1]-b[1], 0],
        [0, c[0]-b[0], 0, a[0]-c[0], 0, b[0]-a[0]],
        [c[0]-b[0], b[1]-c[1], a[0]-c[0], c[1]-a[1], b[0]-a[0], a[1]-b[1]],
    ])/(Area(a,b,c))

def solver(task):
    N = 2*len(task['nodes'])
    K_global = sparse.lil_matrix((N,N))
    F_global = np.zeros((N, 1))

    Nu = task['material']['Nu']
    E = task['material']['E']

    D = D_matrix(Nu)

    for el in task['elems']:
        a = task['nodes'][el[0]]
        b = task['nodes'][el[1]]
        c = task['nodes'][el[2]]

        B = B_matrix(a,b,c)

        K_local = B.transpose().dot(D.dot(B))*Area(a,b,c)

        for i in range(6):
            for j in range(6):
                K_global[2*el[i//2]+i%2, 2*el[j//2]+j%2] += K_local[i, j]


    for bound in task['bcs']:
        if bound['type'] == "neumann":
            
            j = bound['nodes'][0]
            for i in bound['nodes'][1:]:
                a = task['nodes'][j]
                b = task['nodes'][i]

                F_global[2*i,0] += bound['value']*(a[1]-b[1])
                F_global[2*i+1,0] += bound['value']*(b[0]-a[0])
                F_global[2*j,0] += bound['value']*(a[1]-b[1])
                F_global[2*j+1,0] += bound['value']*(b[0]-a[0])

                j = i

        if bound['type'] == "dirichlet":
            for node in bound['nodes']:

                n = 2*node+bound['fix']
                for i in range(N):
                    if i != n:
                        K_global[i,n] = 0
                        K_global[n,i] = 0

                K_global[n,n] = 1
                F_global[n,0] = bound['value']

    U = spsolve(K_global, F_global)

    U.shape = (N//2, 2)

    return Solution(task, U)


def Area(a,b,c):
    return (a[0]*b[1]+b[0]*c[1]+c[0]*a[1]-a[1]*b[0]-b[1]*c[0]-c[1]*a[0])

class Solution:
    def __init__(self, task, displacement):
        self.nodes = np.array(task['nodes'])
        self.elems = task['elems']
        self.len = len(task['nodes'])
        self.displacement = displacement
        self.material = task['material']
        # Вычислить страйн и стресс
    def __call__(self,x,y):

        for elem in self.elems:

            #координаты 
            a = self.nodes[elem[0]]
            b = self.nodes[elem[1]]
            c = self.nodes[elem[2]]
            
            area = Area(a,b,c)

            Nu = self.material['Nu']
            E = self.material['E']

            alp = ((b[1]-c[1])*x + (c[0]-b[0])*y + (b[0]*c[1]-b[1]*c[0]))/area
            bet = ((c[1]-a[1])*x + (a[0]-c[0])*y + (c[0]*a[1]-c[1]*a[0]))/area
            gam = ((a[1]-b[1])*x + (b[0]-a[0])*y + (a[0]*b[1]-a[1]*b[0]))/area

            if alp >= 0 and bet >= 0 and gam >= 0:
                U = np.array([self.displacement[elem[i]] for i in range(3)])

                displacement = (U[0]*alp + U[1]*bet + U[2]*gam)/E

                U.shape = (6, 1)

                strain = B_matrix(a,b,c).dot(U)

                stress = D_matrix(Nu).dot(strain)

                strain = strain/E

                return {
                    "displacement":displacement,
                    "stress":stress,
                    "strain":strain
                }

        return {
            "displacement":[0,0],
            "stress":[0,0,0],
            "strain":[0,0,0]
        }