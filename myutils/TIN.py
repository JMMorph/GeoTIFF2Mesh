# ------------------------------------------------------------------------------------------------------
    #                                     * Triangular Shape *
    #                                   +-------------------+
    #                                   | \                 
    #                                    |   \               
    #                                     |     \
    #                                      |       \
    #                                       |         \
    #                                        |           \
    #                                         |             \
    #                                          |               \
    #                                           +-------------------+
    

import numpy as np
import c4d
import commands
from shapely import MultiPoint
import shapely as sp


def dilate_step(image, kernel):
    # Crear una copia de la imagen para almacenar el resultado de la dilatación
    dilated_image = np.zeros_like(image)

    # Obtener las dimensiones de la imagen y del kernel
    rows, cols = image.shape
    k_rows, k_cols = kernel.shape
    k_size = k_rows*k_cols
    threshold = 1

    # Calcular el padding necesario para aplicar el kernel
    pad_rows = k_rows // 2
    pad_cols = k_cols // 2

    # Aplicar el kernel a cada píxel de la imagen
    ind_i = [i for i in range(pad_rows, rows - pad_rows, k_rows)] + [rows - pad_rows-1]
    ind_j = [j for j in range(pad_cols, cols - pad_cols, k_cols)] + [cols - pad_cols-1]
    for i in ind_i:
        for j in ind_j:
            neighborhood = image[i - pad_rows:i + pad_rows + 1, j - pad_cols:j + pad_cols + 1]
            if sum(sum(neighborhood)) > threshold:
                # Aplicar el kernel a la región vecina
                dilated_image[i - pad_rows:i + pad_rows + 1, j - pad_cols:j + pad_cols + 1] = kernel

    return dilated_image


def createTIN(obj, raster,step, measureMap, arr, arr_red, no_data, center_x, center_z, center_m,
              g_scaler, m_scaler, max_m_scale, solid_base, bottom_surface, uvw_type, x_segments, z_segments,
              min_m, range_m, NX, NZ, bd_points, bd_poly, min_size, tin_seed,
              gradient_threshold, uniform_sampling_threshold, gaussian_sampling_threshold):
    '''
    Function to create a TIN from a raster image
    
    Parameters:
    obj: c4d.BaseObject
        Object to store the TIN
    raster: rasterio.io.DatasetReader
        Raster image
    measureMap: function
        Function to map measures
    arr: numpy.ndarray
        Array with the elevation data
    arr_red: numpy.ndarray
        Array with the red channel data
    no_data: int
        No data value
    center_x: float
        Center of the image in x
    center_z: float
        Center of the image in z
    center_m: float
        Center of the measure (vertical alignment of the TIN)
    g_scaler: float
        Scaler for the x and z axis
    m_scaler: float
        Scaler for the y axis
    max_m_scale: float
        Maximum value for the y axis
    solid_base: bool
        If the base is solid
    bottom_surface: float
        Bottom surface
    uvw_type: str
        Type of UVW map
    x_segments: int
        Number of segments in the x axis
    z_segments: int
        Number of segments in the z axis
    min_m: float
        Minimum value for the y axis
    range_m: float
        Range of the y axis
    NX: int
        Number of pixels in the x axis
    NZ: int
        Number of pixels in the z axis
    bd_points: int
        Number of points in the base
    bd_poly: int
        Number of polygons in the base
    min_size: int
        Minimum side of the raster
    tin_seed: int
        Seed for the random number generator
    gradient_threshold: float
        Threshold for the gradient
    uniform_sampling_threshold: float
        Threshold for the uniform sampling
    gaussian_sampling_threshold: float
        Threshold for the gaussian sampling
        
        
    
    
    '''
    np.random.seed(tin_seed)
    
    data_points = set()         # Set of points with data
    id_point = 0                # counter for points

    for i in range(0, NX):
        for j in range(0, NZ):
            
            xr, zr = raster.transform*(i*step, j*step)
            yr = measureMap(arr_red[i,j], zr, xr)
            
            # Calculate x, y, z scaled and centered
            x, z = (xr - center_x)*g_scaler, (zr - center_z)*g_scaler
            y = (yr - center_m)*g_scaler*m_scaler 

            obj.SetPoint(id_point, c4d.Vector(z, y, x))
            
            if solid_base:
                obj.SetPoint(id_point + bd_points, 
                            c4d.Vector(z, center_m - bottom_surface, x))
                
            if yr != no_data:    
                data_points.add(id_point)

            # Increase counter
            id_point += 1

    
    UVW = c4d.UVWTag(obj.GetPointCount())
    data_polygons = set()       # Set of polygons used
    surface_polys_set = set()   # Set of polygons used in the surface
    
    # Avoid no_data in gradient
    substitute = arr.max()
    if no_data != None:
        arrFilled = np.where(arr != no_data, arr, substitute*4)
    else:
        arrFilled = arr
        
    # Vertical and horizontal gradient
    p1, p2 = np.gradient(arrFilled)
    gradient = np.sqrt(p1**2 + p2**2)
    
    # Avoid no_data in gradient
    if no_data != None:
        border = gradient > substitute/4
        gradient = gradient*(gradient <= substitute)
        
        # Improve the border
        kernel = np.ones((4*step + 1,4*step +1),np.uint8)
        border_red = dilate_step(border, kernel)
        
    else:
        border = 0
        
    gradient = (gradient-gradient.min())/(gradient.max()-gradient.min())*100

    # Identify points with data to create the polygons
    gradient_red = (gradient>gradient_threshold)
    gradient_red = gradient_red[::step,::step] + border_red[::step,::step]
    
    # Add frame to ensure triangles on the border
    frame_step = min_size//100
    gradient_red[0,np.linspace(0, NZ-1, frame_step, dtype=int)] = True
    gradient_red[NX-1,np.linspace(0, NZ-1, frame_step, dtype=int)] = True
    gradient_red[np.linspace(0, NX-1, frame_step, dtype=int),0] = True
    gradient_red[np.linspace(0, NX-1, frame_step, dtype=int),NZ-1] = True
    
    # Sub sampling of the gradient with uniform noise
    gradient_red = gradient_red * (np.random.uniform(0, 1, (NX, NZ)) > uniform_sampling_threshold)
    
    
    # Add noise samples for flat areas
    noise_flat = (np.random.random((gradient_red.shape)) > gaussian_sampling_threshold)*(arr_red != no_data)
    gradient_red += noise_flat
    
    
    # Create surface polygons from points, using the points
    # from the gradient to create the triangles
    pixels = np.where(gradient_red)
    points = MultiPoint([(i, j) for i, j in zip(*pixels)])
    delaunay = sp.delaunay_triangles(points)
    
    id_poly = 0
    for poly in delaunay.geoms:

        triangle = poly.exterior.xy
        
        ii = [int(value) for value in triangle[0]]
        jj = [int(value) for value in triangle[1]]
        vertices_data = set([arr_red[i,j] for i, j in zip(ii[0:3], jj[0:3]) if arr_red[i,j] != no_data])
        center_ii = int(sum(ii[0:3])/3)
        center_jj = int(sum(jj[0:3])/3)
        center_trian = arr_red[center_ii, center_jj]
        
        P1, P2, P3, P4 = [int(i * NZ + j) for i, j in zip(ii, jj)]
        
        # If two or more points have data, create a polygon

        if (len(vertices_data) == 2):
            if arr_red[ii[0],jj[0]] == no_data:
                temp = obj.GetPoint(P1)
                obj.SetPoint(P1, c4d.Vector(temp.x, (obj.GetPoint(P2).y + obj.GetPoint(P2).y)/2, temp.z))
            elif arr_red[ii[1],jj[1]] == no_data:
                temp = obj.GetPoint(P2)
                obj.SetPoint(P2, c4d.Vector(temp.x, (obj.GetPoint(P1).y + obj.GetPoint(P3).y)/2, temp.z))
            elif arr_red[ii[2],jj[2]] == no_data:
                temp = obj.GetPoint(P3)
                obj.SetPoint(P3, c4d.Vector(temp.x, (obj.GetPoint(P1).y + obj.GetPoint(P2).y)/2, temp.z))
        
        if center_trian != no_data and len(vertices_data) >= 3:
            # ----------------------------------------------------------------------- Surface polygons
            obj.SetPolygon(id_poly, c4d.CPolygon(P1, P2, P3))
            
            # -------------------------------------------------------------------------- Base polygons
            if solid_base:
                obj.SetPolygon(id_poly + bd_poly, 
                            c4d.CPolygon(P1 + bd_points, 
                                            P2 + bd_points,
                                            P3 + bd_points,))
            
            data_polygons.add(id_poly)
            
        # Update UVW Map if is planar
        ti = np.array(ii)*step/float(x_segments)
        tj = np.array(jj)*step/float(z_segments)
        
        surface_polys_set.add(id_poly)
        
        if uvw_type == c4d.GEOMESH_UVW_TOP_VIEW :
            
            UP0 = c4d.Vector(ti[0],tj[0],(arr[ii[0],jj[0]] - min_m)/range_m)
            UP1 = c4d.Vector(ti[1],tj[1],(arr[ii[1],jj[1]] - min_m)/range_m)
            VP0 = c4d.Vector(ti[2],tj[2],(arr[ii[2],jj[2]] - min_m)/range_m)
            VP1 = c4d.Vector(ti[0],tj[0],(arr[ii[0],jj[0]] - min_m)/range_m)
        
            UVW.SetSlow(id_poly,UP0,UP1,VP1,VP0)

        id_poly += 1
        
    # ------------------------------------------------------------------------------------------------------
    # Create the Selection tags for the object
    # ------------------------------------------------------------------------------------------------------

    obj = commands.addPolygonSelectionTag(obj, [1 if i in data_polygons else 0 for i in range(obj.GetPolygonCount())])
    
    obj = commands.optimize(obj)
    obj = commands.triangulate(obj)
    
    
    
    if uvw_type == c4d.GEOMESH_UVW_HEIGHT:
        UVW = commands.get_uvw_from_height(obj, max_m_scale)
    
    return obj, UVW
    