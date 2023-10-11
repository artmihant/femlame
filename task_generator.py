
def GenLameTask(mesh, config_material, config_bcs):

    task = {
        'nodes': mesh['nodes'],
        'elems': mesh['elems'],
        'bcs': [
            {
                "nodes":mesh['edges']['inner'],
                'type':'neumann',
                'value': config_bcs['pressure_inner']
            },
            {
                "nodes":mesh['edges']['right'],
                'type':'dirichlet',
                'fix':1,
                'value':0,
            },
            {
                "nodes":mesh['edges']['outer'],
                'type':'neumann',
                'value': config_bcs['pressure_outer']
            },
            {
                "nodes":mesh['edges']['left'],
                'type':'dirichlet',
                'fix':0,
                'value':0,
            },
        ],
        'material': config_material,
    }

    return task
