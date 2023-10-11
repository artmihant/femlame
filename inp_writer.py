
def mesh2inp(mesh):

    file_data = []

    file_data.append("*PART, NAME=Part-Default")
    file_data.append("*NODE")
    
    for i, node in enumerate(mesh['nodes']):
        file_data.append(f"{i}, {node[0]}, {node[1]}")

    file_data.append("*ELEMENT, TYPE=STRI3")

    for i, elem in enumerate(mesh['elems']):
        file_data.append(f"{i}, {elem[0]}, {elem[1]}, {elem[2]}")

    file_data.append("*END PART")
    file_data.append("*ASSEMBLY, NAME=ASSEMBLY1")
    file_data.append("*INSTANCE, NAME=Part-Default_1, PART=Part-Default")
    file_data.append("*END INSTANCE")
    file_data.append("*END ASSEMBLY")

    return '\n'.join(file_data)