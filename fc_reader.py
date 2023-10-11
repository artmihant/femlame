# -*- coding: utf-8 -*-
import json
from base64 import b64decode

import numpy as np
import pandas as pd

from server.lib.element_types import ELEMENT_TYPES


def tolist(x):
    return list(x.reshape(1, x.size)[0])


def decode(src, dtype):
    data = b64decode(src)
    return np.frombuffer(data, dtype)


FC_ELEMENT_TYPES = {}

for ELEMENT_TYPE in ELEMENT_TYPES:
    for i in ELEMENT_TYPE['fc_id']:
        FC_ELEMENT_TYPES[i] = ELEMENT_TYPE



class FCReader:

    def __init__(self, filepath):

        with open(filepath, "r") as f:
            self.src_data = json.load(f)

        self._decode_header()

        self._decode_blocks()
        self._decode_coordinate_systems()
        self._decode_mesh()
        self._decode_settings()

        self._decode_materials()
        self._decode_loads()
        self._decode_restraints()
        self._decode_sets()


    def __getitem__(self,key):
        return getattr(self,key)


    def _decode_header(self):
        self.header = self.get_from_src('header')
        assert self.header


    def _decode_settings(self):
        self.settings = self.get_from_src('settings')
        assert self.settings


    def _decode_blocks(self):
        self.blocks = self.get_from_src('blocks', [])


    def _decode_coordinate_systems(self):

        self.coordinate_systems = self.get_from_src('coordinate_systems', [])

        for cs in self.coordinate_systems:
            cs['dir1']   = decode(cs['dir1'],   np.float64)
            cs['dir2']   = decode(cs['dir2'],   np.float64)
            cs['origin'] = decode(cs['origin'], np.float64)


    def _decode_mesh(self):
        mesh = self.get_from_src('mesh', {})

        elems = {}
        elems['block'] = decode(mesh['elem_blocks'], np.int32)
        elems['order'] = decode(mesh['elem_orders'], np.int32)
        elems['parent_id'] = decode(mesh['elem_parent_ids'], np.int32)
        elems['type'] = decode(mesh['elem_types'], np.int8)
        elems['id'] = decode(mesh['elemids'], np.int32)
        elems['nodes_list'] = decode(mesh['elems'], np.int32)

        elems['nodes'] = []

        counter = 0

        for elem_type in elems['type']:
            count = FC_ELEMENT_TYPES[elem_type]['nodes']
            elems['nodes'].append(elems['nodes_list'][counter:(counter+count)])
            counter += count

        nodes = {}
        nodes['id'] = decode(mesh['nids'], np.int32)
        nodes['xyz'] = decode(mesh['nodes'], np.float64).reshape(-1,3)

        assert mesh['elems_count'] == len(elems['id'])
        assert mesh['nodes_count'] == len(nodes['id'])

        self.mesh = {}

        self.mesh['nodes'] = pd.DataFrame(data=nodes['xyz'], columns=['x','y','z'], index=nodes['id'])

        self.mesh['elems'] = pd.DataFrame(data={
            'block': elems['block'],
            'order': elems['order'],
            'parent_id': elems['parent_id'],
            'type': elems['type'],
            'nodes': elems['nodes']
        }, index=elems['id'])


    def _decode_sets(self):
        self.sets = self.get_from_src('sets', {
            "nodesets": [],
            "sidesets": []
        })

        for nodeset in self.sets['nodesets']:
            nodeset['apply_to'] = decode(nodeset['apply_to'], np.int32).reshape(sideset['apply_to_size'],-1)

        for sideset in self.sets['sidesets']:
            sideset['apply_to'] = decode(sideset['apply_to'], np.int32).reshape(sideset['apply_to_size'],-1)


    def _decode_materials(self):
        self.materials = self.get_from_src('materials', [])

        for material in self.materials:
            for property_name in material:
                material_property = material[property_name]
                if type(material_property)== list:
                    for cs in material_property:
                        cs['constants'] = [ decode(c, np.float64) for c in cs['constants']]


    def _decode_loads(self):
        self.loads = self.get_from_src('loads', [])

        for load in self.loads:

            if load['apply_to'] != 'all':
                load['apply_to'] = decode(load['apply_to'], np.int32).reshape(load['apply_to_size'],-1)

            load['data'] = [decode(value, np.float64) for value in load['data']]


    def _decode_restraints(self):
        self.restraints = self.get_from_src('restraints', [])

        for restraint in self.restraints:
            if restraint['apply_to'] != 'all':
                restraint['apply_to'] = decode(restraint['apply_to'], np.int32).reshape(restraint['apply_to_size'])
            restraint['data'] = [decode(value, np.float64) for value in restraint['data']]


    def get_from_src(self, key, default=None):
        return self.src_data.get(key, default)


    def stream_fragments(self, dim, rank, index_replace=None):
        index_replace_local = {el:i for i,el in enumerate(self.mesh['nodes'].index.tolist())}

        if index_replace is None:
            index_replace = index_replace_local

        fragments = []

        title = None

        for element in self.mesh['elems'].iloc:
            element_type = FC_ELEMENT_TYPES[element['type']]

            if dim < element_type['site'] or element_type['site'] < rank:
                continue;

            element_structure = element_type['structure'][rank]

            element_parts = np.copy(element_structure)
            nodes = element['nodes']

            for i, el in enumerate(nodes):
                element_parts[element_structure==i] = index_replace[el]

            element_parts = element_parts.reshape((-1,rank+1))

            if title and title[1] == element['block'] and title[2] == element_type['site']:
                title[4] += len(element_parts)
            else:
                title = [dim, element['block'], element_type['site'], rank, len(element_parts)]

                fragments.append(title)

            fragments.extend(element_parts)

        stream = []
        for a in fragments:
            stream.extend(a)

        return stream



if __name__ == '__main__':

    test_path = '/home/artem/ProveDesign/alpha/data/storage/fe4fc46cf25c/body/e89a8a6ae254/stl_surface.fc'

    # src_vol_path = '/home/artem/ProveDesign/alpha/notes/fc_tests/simple_cube_vol.fc'
    # src_surf_path = '/home/artem/ProveDesign/alpha/data/storage/b42d324b8e50/body/2d4dd7ff8c73/stl_surface.fc'
    # src_curve_path = '/home/artem/ProveDesign/alpha/notes/fc_tests/simple_cube_curve.fc'
    # src_vertex_path = '/home/artem/ProveDesign/alpha/notes/fc_tests/simple_cube_vertex.fc'


    test_data = FCReader(test_path)

    # fc_vol_data = FCReader(src_vol_path)
    # fc_surf_data = FCReader(src_surf_path)
    # fc_curve_data = FCReader(src_curve_path)
    # fc_vertex_data = FCReader(src_vertex_path)


    # print(fc_vol_data.stream_tetras())
    # print(test_data.stream_fragments(2))
    # print(fc_curve_data.stream_pairs())
    # print(test_data.mesh['elems']['type'].unique())

    # print(fc_data['loads'])
    # print(fc_data['restraints'])

    # for key in fc_data.mesh:
    #     print(key)
    #     print(fc_data.mesh[key])
    #     print()



    # 1) список треугольников с индексами принадлежности к родителям
    #     - размерность родителя и его индекс
    # 2) список ребер сетки с индексами принадлежности к родителям

