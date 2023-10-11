#-*- coding:utf-8 -*-
import json
from mesh_generator import GenLameMesh
from task_generator import GenLameTask
from inp_writer import mesh2inp

from solver import solver
from tester import LameTest

config = {
    "geometry": {
        "radius_min":1,
        "radius_max":2,
    },
    "mesh":{
    	"resolution":4,
    },
    "material":{
        "E":2*10**11,
        "Nu":0.25,
    },
    "bcs":{
        "pressure_inner":0.0,
        "pressure_outer":1.0
    }
}

mesh = GenLameMesh(config["geometry"], config["mesh"])

# with open("task.inp", "w") as write_file:
#     write_file.write(mesh2inp(mesh)) 

task = GenLameTask(mesh, config["material"], config["bcs"])

result = solver(task)
print(result)

LameTest(result, config)

 
# # LameTest(result, param)

# import math

# math.gcd