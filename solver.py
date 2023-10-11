import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve


def D_matrix(E, Nu):
    """ Матрица линейного преобразования компонентов тензора напряжения в тензор деформации """
    return np.array([
        [1, Nu/(1-Nu), 0],
        [Nu/(1-Nu), 1, 0],
        [0, 0, (1-2*Nu)/(2-2*Nu)]
    ])*((1-Nu)/(1+Nu)/(1-2*Nu))*E


def N_matrix(a,b,c):
    """ Матрица функций формы элемента ABC в пределах этого элемента """
    area = Area(a,b,c)

    def N(x,y):
        N_a = ((b[1]-c[1])*x + (c[0]-b[0])*y + (b[0]*c[1]-b[1]*c[0]))/(2*area)
        N_b = ((c[1]-a[1])*x + (a[0]-c[0])*y + (c[0]*a[1]-c[1]*a[0]))/(2*area)
        N_c = ((a[1]-b[1])*x + (b[0]-a[0])*y + (a[0]*b[1]-a[1]*b[0]))/(2*area)

        return np.array([
            [N_a, 0, N_b, 0, N_c, 0],
            [0, N_a, 0, N_b, 0, N_c]
        ])
    
    return N


def B_matrix(a,b,c):
    """ Продифференцированная матрица функций формы элемента ABC"""
    area = Area(a,b,c)

    return np.array([
        [b[1]-c[1], 0, c[1]-a[1], 0, a[1]-b[1], 0],
        [0, c[0]-b[0], 0, a[0]-c[0], 0, b[0]-a[0]],
        [c[0]-b[0], b[1]-c[1], a[0]-c[0], c[1]-a[1], b[0]-a[0], a[1]-b[1]],
    ])/(2*area)


def Area(a,b,c):
    """ Площадь треугольника с вершинами в точках a,b,c """
    return (a[0]*b[1]+b[0]*c[1]+c[0]*a[1]-a[1]*b[0]-b[1]*c[0]-c[1]*a[0])/2


def solver(task):
    N = 2*len(task['nodes'])

    K_global = sparse.lil_matrix((N,N))
    F_global = np.zeros((N, 1))

    Nu = task['material']['Nu']
    E = task['material']['E']

    # мы не используем Е, потому что домножение K_global на большое число отвратительно влияет на точность решения, 
    # а E является линейным коэффициентом задачи. Так что мы просто поделим на него в конце вектор перемещений узлов.
    D = D_matrix(1, Nu) 

    for el in task['elems']:
        
        abc = task['nodes'][el]

        B = B_matrix(*abc)

        K_local = B.transpose().dot(D).dot(B)*Area(*abc)*2

        for i in range(6):
            for j in range(6):
                K_global[2*el[i//2]+i%2, 2*el[j//2]+j%2] += K_local[i, j]


    for bound in task['bcs']:
        if bound['type'] == "neumann":
            
            j = bound['nodes'][0]
            for i in bound['nodes'][1:]:
                a = task['nodes'][j]
                b = task['nodes'][i]

                d = np.array([a[1]-b[1], b[0]-a[0]])*bound['value']

                F_global[(2*i, 2*i+1),0] += d
                F_global[(2*j, 2*j+1),0] += d

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

    U = U.reshape((N//2, 2))/E

    return Solution(task, U)


class Solution:

    def __init__(self, task, node_displacements):
        self.nodes = task['nodes']
        self.elems = task['elems']
        self.len = len(task['nodes'])
        self.node_displacements = node_displacements
        self.material = task['material']

    def __call__(self,x,y):
        """Вычислить значение полей напряжения/перемещения/деформации в заданой точке решения"""

        for elem in self.elems:

            #координаты узлов элемента
            abc = self.nodes[elem]

            Nu = self.material['Nu']
            E = self.material['E']

            Shape = N_matrix(*abc)(x,y)

            if Shape[0,0] >= 0 and Shape[0,2] >= 0 and Shape[0,4] >= 0: # если точка внутри элемента

                node_displacement = self.node_displacements[elem].reshape((6, 1)) # вектор перемещений узлов элемента

                displacement = Shape.dot(node_displacement)

                strain = B_matrix(*abc).dot(node_displacement)

                stress = D_matrix(E, Nu).dot(strain)

                return {
                    "displacement":displacement,
                    "stress":stress,
                    "strain":strain
                }

        return {
            "displacement":np.zeros((2,1)),
            "stress":np.zeros((3,1)),
            "strain":np.zeros((3,1))
        }