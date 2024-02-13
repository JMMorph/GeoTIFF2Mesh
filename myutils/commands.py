import c4d

    
    
def optimize(obj):
    
    opt = c4d.BaseContainer()
    opt[c4d.MDATA_OPTIMIZE_TOLERANCE] = 0.01
    opt[c4d.MDATA_OPTIMIZE_POINTS] = 1
    opt[c4d.MDATA_OPTIMIZE_POLYGONS] = 1
    opt[c4d.MDATA_OPTIMIZE_UNUSEDPOINTS] = 1

    res = c4d.utils.SendModelingCommand( command = c4d.MCOMMAND_OPTIMIZE,
                                    list = [obj],
                                    mode = c4d.MODELINGCOMMANDMODE_ALL,
                                    bc = opt)
    
    return obj

def triangulate(obj):
    res = c4d.utils.SendModelingCommand( command = c4d.MCOMMAND_TRIANGULATE,
                                    list = [obj],
                                    mode = c4d.MODELINGCOMMANDMODE_ALL)
    return obj

def breakPhong(obj, edges):
    
    print('Not implemented yet')
    return obj

def addPolygonSelectionTag(obj, polys):
    poly_tag = c4d.BaseTag(c4d.Tpolygonselection)
    obj.InsertTag(poly_tag)
    poly_sel = poly_tag.GetBaseSelect()
    poly_sel.SetAll(polys)
    return obj


def addEdgeSelectionTag(obj, edges):
    edge_tag = c4d.BaseTag(c4d.Tedgeselection)
    obj.InsertTag(edge_tag)
    edge_sel = edge_tag.GetBaseSelect()
    edge_sel.SetAll(edges)
    return obj

def reverseNormals(obj):
    c4d.utils.SendModelingCommand(command = c4d.MCOMMAND_REVERSENORMALS,
                                  list = [obj],
                                  mode = c4d.MODELINGCOMMANDMODE_ALL)
    return obj

def get_uvw_from_height(obj, max_m_scale = 1.0):

        
    # Crear un nuevo tag de coordenadas UVW
    UVW = c4d.UVWTag(obj.GetPolygonCount())
    if UVW is None:
        return

    # Obtener las coordenadas de los vértices del objeto
    polys = obj.GetAllPolygons()
    
    vertices = obj.GetAllPoints()
    
    min_x_value = max(vertices, key=lambda v: v.x).x
    max_x_value = min(vertices, key=lambda v: v.x).x
    range_x = max_x_value - min_x_value
    

    # Calcular las coordenadas UVW basadas en la posición en el espacio
    for i, poly in enumerate(polys):
        
        
        pa = c4d.Vector((obj.GetPoint(poly.a).x - min_x_value)/range_x,
                        1 - max(obj.GetPoint(poly.a).y, 0)/max_m_scale, 0)
        pb = c4d.Vector((obj.GetPoint(poly.b).x - min_x_value)/range_x,
                        1 - max(obj.GetPoint(poly.b).y, 0)/max_m_scale, 0)
        pc = c4d.Vector((obj.GetPoint(poly.c).x - min_x_value)/range_x, 
                        1 - max(obj.GetPoint(poly.c).y, 0)/max_m_scale, 0)
        pd = c4d.Vector((obj.GetPoint(poly.d).x - min_x_value)/range_x, 
                        1 - max(obj.GetPoint(poly.d).y, 0)/max_m_scale, 0)

        # Asignar las coordenadas UVW al tag UVW del objeto
        UVW.SetSlow(i, pa, pb, pd, pc)
    
    return UVW