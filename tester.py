import time
from math import sqrt
now = time.time()

def LameTest(L_calc, param):
    a = param['geometry']['radius_min']
    b = param['geometry']['radius_max']
    p1 = param['bcs']['pressure_inner']
    p2 = param['bcs']['pressure_outer']
    L_an = AnaliticalLameStress(a,b,p1,p2)

    resolution = param['mesh']['resolution']

    def metr(a,b):
        if a!=0 and b!=0:
            return abs(a-b)/sqrt(a**2+b**2)
        return abs(a-b)


    error = [[],[],[]]

    for i in range(11):

        s2 = 1/sqrt(2)
        r = a+i/10*(b-a)

        s_an = L_an(r)

        s_cl = [float(L_calc(r,0)['stress'][i]) for i in range(3)]
    
        # s_cl = [float(L_calc(0.001,r)['stress'][i]) for i in range(3)]
        # s_cl = [s_cl[1],s_cl[0],s_cl[2]]

        # s_cl = [float(L_calc(r*s2,r*s2)['stress'][i]) for i in range(3)]
        # s_cl = [(s_cl[0] + s_cl[1])/2 + s_cl[2], (s_cl[0] + s_cl[1])/2 - s_cl[2] ,(s_cl[1]-s_cl[0])/2]

        for i in range(3):
            error[i].append(metr(s_cl[i],s_an[i]))

        print(
            round(s_cl[0],3),
            round(s_an[0],3),'\t\t',
            round(s_cl[1],3),
            round(s_an[1],3),'\t\t',
            round(s_cl[2],3),
            round(s_an[2],3),
        )

    global now

    print('Точность', resolution, "(",L_calc.len,"точек )", 'время', round(time.time()-now,1))
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