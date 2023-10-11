import time
from math import sqrt
now = time.time()

def LameTest(solution, config):
    a = config['geometry']['radius_min']
    b = config['geometry']['radius_max']
    p1 = config['bcs']['pressure_inner']
    p2 = config['bcs']['pressure_outer']
    L_an = AnaliticalLameStress(a,b,p1,p2)

    resolution = config['mesh']['resolution']

    def metr(a,b):
        if a!=0 and b!=0:
            return abs(a-b)/sqrt(a**2+b**2)
        return abs(a-b)


    error = [[],[],[]]

    for i in range(11):

        r = a+i/10*(b-a)

        s_an = L_an(r)

        s_cl = [float(solution(r,0)['stress'][i]) for i in range(3)]
    
        for i in range(3):
            error[i].append(metr(s_cl[i],s_an[i]))

        # print(
        #     round(s_cl[0],3),
        #     round(s_an[0],3),'\t\t',
        #     round(s_cl[1],3),
        #     round(s_an[1],3),'\t\t',
        #     round(s_cl[2],3),
        #     round(s_an[2],3),
        # )

    global now

    print('Точность', resolution, "(",solution.len,"точек )", 'время', round(time.time()-now,1))
    print('rr погрешность - ', round(100*sum(error[0])/10,1),'%')
    print('ff погрешность - ', round(100*sum(error[1])/10,1),'%')
    print('rf погрешность - ', round(100*sum(error[2])/10,1),'%')
    now = time.time()


def AnaliticalLameStress(a,b,p1,p2):
    def result(r):
        return [
            p1*a**2/(b**2-a**2)*(1-b**2/r**2) - p2*b**2/(b**2-a**2)*(1-a**2/r**2),
            p1*a**2/(b**2-a**2)*(1+b**2/r**2) - p2*b**2/(b**2-a**2)*(1+a**2/r**2),
            0
        ]
    return result