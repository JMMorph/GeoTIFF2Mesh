# ------------------------------------------------------------------------------------------------------
#                                   * SQUARE BAED SHAPE *    
#                                   +-------------------+
#                                   |\                  |
#                                   |  \                |
#                                   |    \              |
#                                   |      \            |
#                                   |        \          |
#                                   |          \        |
#                                   |            \      |
#                                   |              \    |
#                                   |                \  |
#                                   |                  \|
#                                   +-------------------+


import numpy as np
import c4d
import commands

def quadFromIndex(i, j, NX, NZ):
    
    # Define points for a square
    P1 = int(i * NZ + j)
    P2 = int(i * NZ + j + 1)
    P3 = int((i + 1) * NZ + j + 1)
    P4 = int((i + 1) * NZ  + j)
    
    return P1, P2, P3, P4

def createTRN(obj, raster, measureMap, solid_base, uvw_type, step, g_scaler, 
              m_scaler, max_m_scale, center_x, center_z, center_m, surface_thickness, no_data, min_m, arr,
              NX, NZ, all_no_data, surface_polys, bd_points, bd_poly):
    
    '''
    Function to create a Triangular Regular Network (TRN) from a raster file.
    
    Parameters:
    - obj: c4d.PolygonObject
        Object to store the TRN
    - raster: rasterio dataset
        Raster file with the data
    - measureMap: function
        Function to convert the raster data to the desired measure
    - solid_base: bool
        If True, the base of the TRN will be solid
    - uvw_type: str
        Type of UVW mapping for the TRN (top view or height based)
    - step: int
        Step to sample the raster
    - g_scaler: float
        Scaler for all the coordinates
    - m_scaler: float
        Scaler for the measure (height)
    - max_m_scaler: float
        Scaler for the maximum value of the measure
    - center_x: float
        Center of the X axis
    - center_z: float
        Center of the Z axis
    - center_m: float
        Center of the measure (vertical alignment of the TRN)
    - surface_thickness: float
        Bottom surface of the TRN
    - no_data: float
        Value for no data in the raster
    - min_m: float
        Minimum value of the measure
    - max_m_scaler: float
        Scaler for the maximum value of the measure
    - arr: numpy array
        Array with the raster data
    - NX: int
        Number of cells in the X axis
    - NZ: int
        Number of cells in the Z axis
    - all_no_data: int
        Number of polygons with no data
    - surface_polys: int
        Number of polygons with data
    - bd_points: int
        Number of points in the surface
    - bd_poly: int
        Number of polygons in the surface
    
    Returns:
    - obj: c4d.PolygonObject
        Object with the TRN
    - UVW: c4d.UVWTag
        UVW tag for the TRN
    
    '''
        
        
    # ------------------------------------------------------------------------------------------------------
    # Create all points
    # ------------------------------------------------------------------------------------------------------
    
        
    data_points = set()         # Set of points with data
    id_point = 0                # counter for points
    for i in range(0, NX*step, step):
        for j in range(0, NZ*step, step):
            
            xr, zr = raster.transform*(i, j)
            yr = measureMap(arr[i,j], zr, xr)
            
            # Calculate x, y, z scaled and centered
            x, z = (xr - center_x)*g_scaler, (zr - center_z)*g_scaler
            y = (yr - center_m)*g_scaler*m_scaler 

            # Save points
            if yr != no_data:
                
                obj.SetPoint(id_point, c4d.Vector(z, y, x))
                
                if solid_base:
                    obj.SetPoint(id_point + bd_points, 
                                c4d.Vector(z, min_m - center_m - surface_thickness, x))
                    
                data_points.add(id_point)

            # Increase counter
            id_point += 1

    # ------------------------------------------------------------------------------------------------------
    # Create polygons of the surface (and the base if solid_base is True)
    # ------------------------------------------------------------------------------------------------------        

    id_poly = 0 # reset counter
    UVW = c4d.UVWTag(obj.GetPointCount())
    border_edges = set()
    border_poly = []
    data_polygons = set()       # Set of polygons with 4 points with data
    surface_polys_set = set()   # Set of polygons used in the surface
    
    # Create surface polygons from points
    for j in range (0, NZ-1):
        for i in range (0, NX-1):

            # Define points for a square
            P1, P2, P3, P4 = quadFromIndex(i, j, NX, NZ)
            
            # Number of points with data in the square
            poly_data_points = len(set([P1, P2, P3, P4]).intersection(data_points))
            
            # If all points have data, create a polygon and update UVW Map
            if poly_data_points == 4:
                
                # ----------------------------------------------------------------------- Surface polygons
                obj.SetPolygon(id_poly, c4d.CPolygon(P4, P3, P2, P1))
                
                # -------------------------------------------------------------------------- Base polygons
                if solid_base:
                    obj.SetPolygon(id_poly + bd_poly, 
                                c4d.CPolygon(P1 + bd_points, 
                                                P2 + bd_points,
                                                P3 + bd_points,
                                                P4 + bd_points))
                
                data_polygons.add(id_poly)
                

                # Identify edges that are on the border (borders of the raster)
                if i == 0:
                    border_edges.update([4*id_poly+2])   # -X
                if i == NX-2:
                    border_edges.update([4*id_poly])     # +X
                if j == 0:
                    border_edges.update([4*id_poly+3])   # +Z
                if j == NZ-2:
                    border_edges.update([4*id_poly+1])   # -Z
        
            elif poly_data_points > 0 and poly_data_points < 4:

                # Identify edges that are on the border (borders with no data points)
                border_edges.update([4*(id_poly + 1) + 2,           # -X 
                                    4*(id_poly - 1)])               # +X
                border_edges.update([4*(id_poly + (NX-1)) + 3,      # +Z 
                                    4*(id_poly - (NX-1))+1])        # -Z
                
                # Identify polygons that are on the border (shared with no data points
                border_poly.append(id_poly)
            
            
            
            # Update UVW Map
            du0 = float(1)/NX * i
            dv0 = float(1)/NZ * j
            du1 = du0 + float(1)/NX
            dv1 = dv0 + float(1)/NZ
                
            if uvw_type == c4d.GEOMESH_UVW_TOP_VIEW:

                UP0 = c4d.Vector(du0,dv0,0)
                UP1 = c4d.Vector(du0,dv1,0)
                VP0 = c4d.Vector(du1,dv0,0)
                VP1 = c4d.Vector(du1,dv1,0)
                
                UVW.SetSlow(id_poly,UP0,UP1,VP1,VP0)
                UVW.SetSlow(id_poly + bd_poly,UP0,UP1,VP1,VP0)
            
            surface_polys_set.add(id_poly)
            # Increase counter
            id_poly += 1

    
    # ------------------------------------------------------------------------------------------------------
    # Create the polygons for the solid border 
    # ------------------------------------------------------------------------------------------------------
    
    if solid_base:
        
        # ------------------------------------------------------------------------------------ Wall polygons
        wall_min_x = [id_poly for id_poly in range(NX-1) if id_poly in data_polygons]
        wall_max_x = [id_poly for id_poly in range((NZ-2)*(NX-1), (NX-1)*(NZ-1)) if id_poly in data_polygons]
        wall_min_z = [id_poly for id_poly in range(0, (NX-1)*(NZ-1), NX-1) if id_poly in data_polygons]
        wall_max_z = [id_poly for id_poly in range(NX-2, (NX-1)*(NZ-1), NX-1) if id_poly in data_polygons]
        
        count = 0
        for id_poly in wall_min_x + wall_max_x + wall_min_z + wall_max_z:
            i = int(id_poly % (NX - 1))
            j = int(id_poly // (NX - 1))
            P1, P2, P3, P4 = quadFromIndex(i, j, NX, NZ)

            # +Z
            if j == 0:
                obj.SetPolygon(surface_polys*2 + 4*all_no_data + count, 
                            c4d.CPolygon(P4, P1, P1+bd_points, P4+bd_points))
                count += 1
            # -Z
            if j == NZ-2:
                obj.SetPolygon(surface_polys*2 + 4*all_no_data + count , 
                            c4d.CPolygon(P2, P3, P3+bd_points, P2+bd_points))
                count += 1
            # -X
            if i == 0:
                obj.SetPolygon(surface_polys*2 + 4*all_no_data + count , 
                        c4d.CPolygon(P1, P2, P2+bd_points, P1+bd_points))
                count += 1
            # +X
            if i == NX-2:
                obj.SetPolygon(surface_polys*2 + 4*all_no_data + count , 
                        c4d.CPolygon(P3, P4, P4+bd_points, P3+bd_points))
                count += 1
        
            
        # ------------------------------------------------------------------------------------ No data polygons
        if all_no_data > 0:
            
            count = 0
            for id_poly in border_poly:
                i = int(id_poly % (NX - 1))
                j = int(id_poly // (NX - 1))
                P1, P2, P3, P4 = quadFromIndex(i, j, NX, NZ)
                
                # +X
                if id_poly >= 1:
                    edge1 = obj.GetPolygon(id_poly - 1).FindEdge(P1, P2)
                    if edge1 != -1:
                        obj.SetPolygon(surface_polys*2 + count, 
                                    c4d.CPolygon(P2, P1, P1+bd_points, P2+bd_points))
                        
                # +Z
                if id_poly < surface_polys - (NX-1):
                    edge2 = obj.GetPolygon(id_poly + (NX-1)).FindEdge(P2, P3)
                    if edge2 != -1:
                        obj.SetPolygon(surface_polys*2 + all_no_data + count, 
                                    c4d.CPolygon(P3, P2, P2+bd_points, P3+bd_points))
                        
                # -X        
                if id_poly < surface_polys - 1:            
                    edge3 = obj.GetPolygon(id_poly + 1).FindEdge(P3, P4)
                    if edge3 > -1:
                        obj.SetPolygon(surface_polys*2 + 2*all_no_data + count, 
                                    c4d.CPolygon(P4, P3, P3+bd_points, P4+bd_points))
                
                # -Z
                if id_poly > (NX-1):
                    edge4 = obj.GetPolygon(id_poly - (NX-1)).FindEdge(P4, P1)
                    if edge4 != -1:
                        obj.SetPolygon(surface_polys*2 + 3*all_no_data + count, 
                                    c4d.CPolygon(P1, P4, P4+bd_points, P1+bd_points))
                
                count += 1
        
    # ------------------------------------------------------------------------------------------------------
    # Create the Selection tags for the object
    # ------------------------------------------------------------------------------------------------------
    
    # binary list of border edges
    border_edges_list_binary = [1 if i in border_edges else 0 for i in range((NX-1)*(NZ-1)*4)]

    obj = commands.addPolygonSelectionTag(obj, [1 for _ in range(surface_polys)])
    obj = commands.addEdgeSelectionTag(obj, border_edges_list_binary)
    
    # obj = commands.breakPhong(obj, border_edges)
    obj = commands.optimize(obj)
    obj = commands.triangulate(obj)
    obj = commands.reverseNormals(obj)
    
    
    if uvw_type == c4d.GEOMESH_UVW_HEIGHT:
        UVW = commands.get_uvw_from_height(obj, max_m_scale)
    
    return obj, UVW
    