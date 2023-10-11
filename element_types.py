from typing import TypedDict

import numpy as np

class ElementType(TypedDict):
    type : str
    fc_id : list[int]
    vtk_id: list[int]
    site: int
    order: int
    nodes: int
    structure: dict[int, np.ndarray]
    edges: list[list[int]]
    facets: list[list[int]]
    tetras: list[list[int]]



# https://www.cs.auckland.ac.nz/compsci716s2t/resources/VTK_file-formats.pdf
# https://kitware.github.io/vtk-examples/site/VTKFileFormats/
# vtk_code: [name, count_nodes, edges_struct, trangles_struct]
#file:///home/artem/.bin/CAE-Fidesys-5.1/preprocessor/bin/help/fidesyshelp.htm

ELEMENT_TYPES: list[ElementType] = [
    {
        'type': 'VERTEX1',
        'fc_id': [38, 99],
        'vtk_id': [1],
        'site': 0,
        'order': 1,
        'nodes': 1,
        'edges': [],
        'facets': [],
        'tetras': [],
        'structure': {}
    },
    {
        'type': 'BEAM2',
        'fc_id': [36,89],
        'vtk_id': [3],
        'site': 1,
        'order': 1,
        'nodes': 2,
        'edges': [[0,1]],
        'facets': [],
        'tetras': [],
        'structure': {}
    },
    {
        'type': 'BEAM3',
        'fc_id': [0],
        'vtk_id': [21],
        'site': 1,
        'order': 2,
        'nodes': 3,
        'edges': [[0,2,1]],
        'facets': [],
        'tetras': [],
        'structure': {}
    },
    {
        'type': 'TRI3',
        'fc_id': [29],
        'vtk_id': [5],
        'site': 2,
        'order': 1,
        'nodes': 3,
        'edges': [[0,1,2,0]],
        'facets': [[0,1,2]],
        'tetras': [],
        'structure': {}
    },
    {
        'type': 'TRI6',
        'fc_id': [0],
        'vtk_id': [22],
        'site': 2,
        'order': 2,
        'nodes': 6,
        'edges': [[0,3,1,4,2,5,0]],
        'facets': [[0,3,1,4,2,5]],
        'tetras': [],
        'structure': {}
    },
    {
        'type': 'QUAD4',
        'fc_id': [31],
        'vtk_id': [8,9],
        'site': 2,
        'order': 1,
        'nodes': 4,
        'edges': [[0,1,2,3,0]],
        'facets': [[0,1,2,3]],
        'tetras': [],
        'structure': {}
    },
    {
        'type': 'QUAD8',
        'fc_id': [0],
        'vtk_id': [23],
        'site': 2,
        'order': 2,
        'nodes': 8,
        'edges': [[0,4,1,5,2,6,3,7,0]],
        'facets': [[0,4,1,5,2,6,3,7]],
        'tetras': [],
        'structure': {}
    },
    {
        'type': 'TETRA4',
        'fc_id': [1],
        'vtk_id': [10],
        'site': 3,
        'order': 1,
        'nodes': 4,
        'edges': [[0, 1, 2, 0], [0, 3], [1, 3], [2, 3]],
        'facets': [[0, 2, 1], [0, 1, 3], [1, 2, 3], [2, 0, 3]],
        'tetras': [[0, 1, 2, 3]],
        'structure': {}
    },
    {
        'type': 'TETRA10',
        'fc_id': [0],
        'vtk_id': [24],
        'site': 3,
        'order': 2,
        'nodes': 10,
        'edges': [[0, 4, 1, 5, 2, 6, 0], [0, 7, 3], [1, 8, 3], [2, 9, 3]],
        'facets': [[0, 6, 2, 5, 1, 4], [0, 4, 1, 8, 3, 5], [1, 5, 2, 9, 3, 8], [2, 6, 0, 5, 3, 9]] ,
        'tetras': [],
        'structure': {}
    },
    {
        'type': 'HEX8',
        'fc_id': [3],
        'vtk_id': [11,12],
        'site': 3,
        'order': 1,
        'nodes': 8,
        'edges': [[0, 1, 2, 3, 0], [4, 5, 6, 7, 4], [0, 4], [1, 5], [2, 6], [3, 7]],
        'facets': [[3, 2, 1, 0], [4, 5, 6, 7], [1, 2, 6, 5], [0, 1, 5, 4], [0, 4, 7, 3], [2, 3, 7, 6]],
        'tetras': [[1,3,4,6],[3,1,4,0],[1,3,6,2],[4,1,6,5],[3,4,6,7]],
        'structure': {}
    },
    {
        'type': 'HEX20',
        'fc_id': [0],
        'vtk_id': [25],
        'site': 3,
        'order': 2,
        'nodes': 20,
        'edges': [[0, 8, 1, 9, 2, 10, 3, 11, 0], [4, 12, 5, 13, 6, 14, 7, 15, 4],
            [0, 16, 4], [1, 17, 5], [2, 18, 6], [3, 19, 7]],
        'facets': [[3, 10, 2, 9, 1, 8, 0, 11], [4, 12, 5, 13, 6, 14, 7, 15], [1, 9, 2, 18, 6, 13, 5, 17],
            [0, 8, 1, 17, 5, 12, 4, 16], [0, 16, 4, 15, 7, 19, 3, 11], [2, 10, 3, 19, 7, 14, 6, 18]],
        'tetras': [],
        'structure': {}
    },
    {
        'type': 'WEDGE6',
        'fc_id': [0],
        'vtk_id': [13],
        'site': 3,
        'order': 1,
        'nodes': 5,
        'edges': [[0, 1, 2, 0], [3, 4, 5, 3], [0, 3], [1, 4], [2, 5]],
        'facets': [[0, 1, 2], [5, 4, 3], [0, 2, 5, 3], [0, 3, 4, 1], [1, 4, 5, 2]],
        'tetras': [[0,5,4,3],[0,4,2,1],[0,2,4,5]],
        'structure': {}
    },
    {
        'type': 'WEDGE15',
        'fc_id': [0],
        'vtk_id': [26],
        'site': 3,
        'order': 2,
        'nodes': 15,
        'edges': [[0, 5, 1, 6, 2, 7, 3, 8, 0], [0, 9, 4], [1, 10, 4], [2, 11, 4], [3, 12, 4]],
        'facets': [[3, 7, 2, 6, 1, 5, 0, 8],
            [0, 5, 1, 10, 4, 9], [1, 6, 2, 11, 4, 10], [2, 7, 3, 12, 4, 11], [3, 8, 0, 9, 4, 12]],
        'tetras': [],
        'structure': {}
    },
    {
        'type': 'PYRAMID5',
        'fc_id': [8],
        'vtk_id': [14],
        'site': 3,
        'order': 1,
        'nodes': 5,
        'edges': [[0, 1, 2, 3, 0], [0, 4], [1, 4], [2, 4], [3, 4]],
        'facets': [[3, 2, 1, 0], [0, 1, 4], [1, 2, 4], [2, 3, 4], [3, 0, 4]],
        'tetras': [[1,3,4,0], [3,4,1,2]],
        'structure': {}
    },
    {
        'type': 'PYRAMID13',
        'fc_id': [0],
        'vtk_id': [27],
        'site': 3,
        'order': 2,
        'nodes': 13,
        'edges': [[0, 5, 1, 6, 2, 7, 3, 8, 0], [0, 9, 4], [1, 10, 4], [2, 11, 4], [3, 12, 4]],
        'facets': [[3, 7, 2, 6, 1, 5, 0, 8],
            [0, 5, 1, 10, 4, 9], [1, 6, 2, 11, 4, 10], [2, 7, 3, 12, 4, 11], [3, 8, 0, 9, 4, 12]],
        'tetras': [],
        'structure': {}
    },
]


def split_facet(facet:list[int]) -> list[int]:
    if len(facet) == 3:
        return facet
    if len(facet) < 3:
        return []
    tail = facet[2:]
    tail.append(facet[1])
    tris = [facet[-1],facet[0],facet[1]]
    tris.extend(split_facet(tail))
    return tris


def split_edge(edge:list[int]) -> list[int]:
    if len(edge) == 2:
        return edge
    if len(edge) < 2:
        return []
    tail = edge[1:]
    pairs = [edge[0],edge[1]]
    pairs.extend(split_edge(tail))
    return pairs

def split_polihedron(tetra:list[int]) -> list[int]:
    return tetra


for ELEMENT_TYPE in ELEMENT_TYPES:

    ELEMENT_TYPE['structure'][0] = np.arange(ELEMENT_TYPE['nodes'], dtype=np.int32)

    if ELEMENT_TYPE['site'] > 0:

        pairs = []
        for edge in ELEMENT_TYPE['edges']:
            pairs.extend(split_edge(edge))

        ELEMENT_TYPE['structure'][1] = np.array(pairs, dtype=np.int32)

    if ELEMENT_TYPE['site'] > 1:

        trangles = []
        for facet in ELEMENT_TYPE['facets']:
            trangles.extend(split_facet(facet))

        ELEMENT_TYPE['structure'][2] = np.array(trangles, dtype=np.int32)


    if ELEMENT_TYPE['site'] > 2:

        tetras = []
        for tetra in ELEMENT_TYPE['tetras']:
            tetras.extend(split_polihedron(tetra))

        ELEMENT_TYPE['structure'][3] = np.array(tetras, dtype=np.int32)


