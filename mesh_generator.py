#-*- coding:utf-8 -*-
from math import pi, sin, cos


def polar2сartesian(r, fi):
    """ Перевод из полярной системы координат в декартову"""
    fi = fi*pi/2
    return [r*cos(fi), r*sin(fi)]


def GenLameMesh(config_geometry, config_mesh):

    # разрешение (мелкость) сетки #
    resolution = config_mesh['resolution']

    # максимальный и минимальный радиус трубы
    radius_min, radius_max = config_geometry['radius_min'], config_geometry['radius_max']


    nodes = []
    elems = []

    edges = {
        'inner': [],
        'right': [],
        'outer': [],
        'left': []
    }


    # число слоев сетки вдоль радиуса
    mesh_layers_count = int(2*resolution/pi*(radius_max-radius_min)/radius_min) + 1 

    # типичный размер элемента сетки / расстояние между слоями
    delta = (radius_max-radius_min)/mesh_layers_count

    #генерируем узлы сетки и добавляем их на края

    node_index = 0
    for r in range(mesh_layers_count+1):
        for f in range(resolution+r+1):
            coord = polar2сartesian(radius_min+delta*r, float(f)/(resolution+r))
            nodes.append(coord)

            if r==0:
                edges['inner'].append(node_index)
            if f==0:
                edges['right'].append(node_index)
            if r==mesh_layers_count:
                edges['outer'].append(node_index)
            if f==(resolution+r):
                edges['left'].append(node_index)

            node_index += 1

    edges['inner'].reverse()


    #генерируем треугольные элементы сетки

    for r in range(mesh_layers_count+1):
        for f in range(resolution+r):
            a = (2*resolution+r+1)*r//2 + f
            b = (2*resolution+r+1)*r//2 + f + 1
            c1 = (2*resolution+r+2)*(r+1)//2 + f + 1
            c2 = (2*resolution+r)*(r-1)//2 + f
            if r != 0:
                elems.append([a,b,c2])
            if r != mesh_layers_count:
                elems.append([a,c1,b])


    return {
        'nodes': nodes,
        'elems': elems,
        'edges': edges,
    }

