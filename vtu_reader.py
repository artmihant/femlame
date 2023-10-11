from copy import copy

import base64
import zlib

import numpy as np
import xml.etree.cElementTree as ET


class VTUReader:


    vtk_cell_types = {
        # https://www.cs.auckland.ac.nz/compsci716s2t/resources/VTK_file-formats.pdf
        # https://kitware.github.io/vtk-examples/site/VTKFileFormats/
        # vtk_code: (name, count_nodes, edges_struct, trangles_struct)
        1: ("vertex", 1, [], []),
        3: ("edge", 2, [(0, 1)], []),
        5: ("triangle", 3,
            [(0, 1, 2, 0)],
            [(0, 1, 2)]),
        9: ("quad", 4,
            [(0, 1, 2, 3, 0)],
            [(0, 1, 2, 3)]),
        10: ("tetra", 4,
             [(0, 1, 2, 0), (0, 3), (1, 3), (2, 3)],
             [(0, 2, 1), (0, 1, 3), (1, 2, 3), (2, 0, 3)]),
        12: ("hexahedron", 8,
             [(0, 1, 2, 3, 0), (4, 5, 6, 7, 4), (0, 4), (1, 5), (2, 6), (3, 7)],
             [(3, 2, 1, 0), (4, 5, 6, 7), (1, 2, 6, 5), (0, 1, 5, 4), (0, 4, 7, 3), (2, 3, 7, 6)]),
        13: ("wedge", 6,
             [(0, 1, 2, 0), (3, 4, 5, 3), (0, 3), (1, 4), (2, 5)],
             [(0, 1, 2), (5, 4, 3), (0, 2, 5, 3), (0, 3, 4, 1), (1, 4, 5, 2)]),
        14: ("pyramid", 5,
             [(0, 1, 2, 3, 0), (0, 4), (1, 4), (2, 4), (3, 4)],
             [(3, 2, 1, 0), (0, 1, 4), (1, 2, 4), (2, 3, 4), (3, 0, 4)]),
        21: ("edge3", 3, [(0, 2, 1)], []),
        22: ("triangle6", 6,
             [(0, 3, 1, 4, 2, 5, 0)],
             [(0, 3, 1, 4, 2, 5)]),
        23: ("quad8", 8,
             [(0, 4, 1, 5, 2, 6, 3, 7, 0)],
             [(0, 4, 1, 5, 2, 6, 3, 7)]),
        24: ("tetra10", 10,
             [(0, 4, 1, 5, 2, 6, 0), (0, 7, 3), (1, 8, 3), (2, 9, 3)],
             [(0, 6, 2, 5, 1, 4), (0, 4, 1, 8, 3, 5), (1, 5, 2, 9, 3, 8), (2, 6, 0, 5, 3, 9)]),
        25: ("hexahedron20", 20,
             [(0, 8, 1, 9, 2, 10, 3, 11, 0), (4, 12, 5, 13, 6, 14, 7, 15, 4),
              (0, 16, 4), (1, 17, 5), (2, 18, 6), (3, 19, 7)],
             [(3, 10, 2, 9, 1, 8, 0, 11), (4, 12, 5, 13, 6, 14, 7, 15), (1, 9, 2, 18, 6, 13, 5, 17),
              (0, 8, 1, 17, 5, 12, 4, 16), (0, 16, 4, 15, 7, 19, 3, 11), (2, 10, 3, 19, 7, 14, 6, 18)]),
        26: ("wedge15", 15,
             [(0, 6, 1, 7, 2, 8, 0), (3, 9, 4, 10, 5, 11, 3), (0, 12, 3), (1, 13, 4), (2, 14, 5)],
             [(0, 6, 1, 7, 2, 8), (5, 10, 4, 9, 3, 11),
              (0, 8, 2, 14, 5, 11, 3, 12), (0, 12, 3, 9, 4, 13, 1, 6), (1, 13, 4, 10, 5, 14, 2, 7)]),
        27: ("pyramid13", 13,
             [(0, 5, 1, 6, 2, 7, 3, 8, 0), (0, 9, 4), (1, 10, 4), (2, 11, 4), (3, 12, 4)],
             [(3, 7, 2, 6, 1, 5, 0, 8),
              (0, 5, 1, 10, 4, 9), (1, 6, 2, 11, 4, 10), (2, 7, 3, 12, 4, 11), (3, 8, 0, 9, 4, 12)]),
    }


    def __init__(self, filename: str):

        with open(filename, 'r') as file_read:
            self.tree = self.scan_file(file_read)

        PointData = self.tree['Piece']['PointData']
        CellData = self.tree['Piece']['CellData']
        Points = self.tree['Piece']['Points']
        Cells = self.tree['Piece']['Cells']

        point_data = copy(PointData)

        for key in point_data:
            if 'Node' in key:
                del PointData[key]

        self.cells = {
            **CellData, 
            "types":Cells['types'],
            "nodes": np.split(Cells['connectivity'], Cells['offsets'])
        }

        self.points = Points['Points']['data']
        self.point_data = PointData
        for attr in self.point_data:
            attr_data = self.point_data[attr]
            if 'X' in attr_data and 'Y' in attr_data and 'Z' in attr_data and 'Magnitude' not in attr_data:
                self.point_data[attr]['Magnitude'] = (attr_data['X']**2 + attr_data['Y']**2 + attr_data['Z']**2)**(1/2)

        pass

    @property
    def elements_count(self):
        elements = {}
        for element_type_code in self.cells["types"]:
            assert element_type_code in self.vtk_cell_types, f'unknown vkt element {element_type_code}'
            element_type_name = self.vtk_cell_types[element_type_code][0]
            if element_type_name not in elements:
                elements[element_type_name] = 1
            else:
                elements[element_type_name] += 1
        return elements


    def scan_file(self, file_read):

        head = ""
        for line in file_read:
            if '<AppendedData' in line:
                head += '</VTKFile>'
                break
            head += line
        head = ET.fromstring(head)

        assert head.attrib["type"] == "UnstructuredGrid"

        UnstructuredGrid = head[0]

        offset_map = self._scan_offset_map(UnstructuredGrid)

        while file_read.read(1) == ' ':
            pass

        return self._scan_tree(UnstructuredGrid, offset_map, file_read)

    def _scan_offset_map(self, xml):
        offsets = []
        for leaf in xml:
            if 'offset' in leaf.attrib:
                offsets.append(int(leaf.attrib['offset']))
            else:
                offsets.extend(self._scan_offset_map(leaf))
        return offsets

    def _scan_tree(self, xml, offset_map, source_file):
        data = {}
        for leaf in xml:
            if 'offset' in leaf.attrib:
                assert leaf.tag in ['Array', 'DataArray']
                assert 'Name' in leaf.attrib

                offset = int(leaf.attrib['offset'])
                index = offset_map.index(offset)
                if index < len(offset_map) - 1:
                    block = source_file.read(
                        offset_map[offset_map.index(offset) + 1] - offset)
                else:
                    block = source_file.readline()[:-1]

                if leaf.tag == 'DataArray':
                    data[leaf.attrib['Name']] = self._convert_data_array(block, leaf.attrib)
            else:
                data[leaf.tag] = self._scan_tree(leaf, offset_map, source_file)
        return data

    @staticmethod
    def _decode_data_block(data, data_type, header_type='UInt32'):
        # using dark magic
        def vtu_to_np_type(name):
            return np.dtype(getattr(np, name.lower()))

        def num_bytes_to_num_base64_chars(num_bytes):
            return -(-num_bytes // 3) * 4

        dtype = vtu_to_np_type(header_type)
        num_bytes_per_item = np.dtype(dtype).itemsize
        num_chars = num_bytes_to_num_base64_chars(num_bytes_per_item)
        byte_string = base64.b64decode(data[:num_chars])[:num_bytes_per_item]
        num_blocks = np.frombuffer(byte_string, dtype)[0]

        num_header_items = 3 + num_blocks
        num_header_bytes = num_bytes_per_item * num_header_items
        num_header_chars = num_bytes_to_num_base64_chars(num_header_bytes)
        byte_string = base64.b64decode(data[:num_header_chars])
        header = np.frombuffer(byte_string, dtype)

        block_sizes = header[3:]

        byte_array = base64.b64decode(data[num_header_chars:])
        dtype = vtu_to_np_type(data_type)

        byte_offsets = np.empty(block_sizes.shape[0] + 1, dtype=block_sizes.dtype)
        byte_offsets[0] = 0
        np.cumsum(block_sizes, out=byte_offsets[1:])

        # process the compressed data
        block_data = np.concatenate([
            np.frombuffer(
                zlib.decompress(byte_array[byte_offsets[k]: byte_offsets[k + 1]]),
                dtype=dtype,
            )
            for k in range(num_blocks)
        ])

        return block_data

    def _convert_data_array(self, block, attrib):
        data_block = self._decode_data_block(block, attrib['type'], "UInt32")

        if 'NumberOfComponents' in attrib:
            dim = int(attrib['NumberOfComponents'])
            columns = list(range(dim))
            for key in attrib:
                if key.startswith('ComponentName'):
                    columns[int(key[len('ComponentName'):])] = attrib[key]

            data = np.reshape(data_block, (-1, dim))

            if dim == 3 and columns[0] == 0:

                if 'Name' in attrib and attrib['Name'] == 'Points':
                    columns = ['X', 'Y', 'Z']
                else:
                    columns = ['Magnitude', 'X', 'Y', 'Z']
                    m = np.apply_along_axis(
                        lambda x: [(x[0] ** 2 + x[1] ** 2 + x[2] ** 2) ** 0.5],
                        1, data)
                    data = np.hstack([m, data])

            return {'columns':columns, 'data':data}

        return data_block


    @property
    def extremes(self):

        extremes = {}
        attrs = self.attrs
        for attr in attrs:
            extremes[attr] = {}
            for i, axis in enumerate(attrs[attr]):
                column = self.point_data[attr]['data'][:,i]

                argmin = column.argmin()
                argmax = column.argmax()

                extremes[attr][axis] = {
                    'min': {
                        'node': argmin,
                        'point': self.points[argmin].tolist(),
                        'value': column.min()
                    },
                    'max': {
                        'node': argmax,
                        'point': self.points[argmax].tolist(),
                        'value': column.max()
                    },
                }

        return extremes


    @property
    def attrs(self):
        return {
            attr: self.point_data[attr]['columns'] 
            for attr in self.point_data
            if 'columns' in self.point_data[attr]
        }



if __name__ == '__main__':

    vtu_input = './input/calc.vtu'

    result_vtu = VTUReader(vtu_input)

    # points = result_vtu.points
    print(result_vtu.extremes)
    