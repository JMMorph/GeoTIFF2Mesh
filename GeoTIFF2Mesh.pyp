'''
Developed by J. Miguel Medina @JMMorph

ParaSurface is a Cinema 4D plugin that allows you to create parametric surfaces from mathematical expressions.

'''

import c4d
from c4d import plugins, bitmaps, Vector, SplineObject, utils

import math
from math import *
import os
import sys
import ast
import copy

# Third-party libraries
from shapely import GeometryCollection, LineString, MultiPoint, Polygon, affinity
import copy
import numpy as np
import rasterio as rio
import shapely as sp

root = os.path.dirname(__file__)
sys.path.insert(0, os.path.join(root, "myutils"))

import TRN
import TIN

EARTH_RADIUS = 6378137
DEFAULT_FILE = os.path.join(root, "res\GeoData\Example_popocatepetl_EPSG3857.tif")

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



# Class to create the plugin
class GeoTIFF2Mesh(plugins.ObjectData):
    
    def __init__ (self):
        
        # Enable cache
        self.SetOptimizeCache(True)     
        
        # Initialize the attributes
        self.raster = None
        self.arr = None
        self.arr_red = None


    def Init(self, op, isCloneInit=False):
            
        return True


    def Draw(self, op, drawpass, bd, bh):
        
        
        if drawpass != c4d.DRAWPASS_OBJECT:
            return c4d.DRAWRESULT_OK

    
        units = op[c4d.GEOTIFF_UNITS] if op[c4d.GEOTIFF_UNITS] is not None else c4d.GEOTIFF_UNITS_METERS
        no_data = op[c4d.GEOTIFF_NO_DATA] if op[c4d.GEOTIFF_NO_DATA] is not None else 0
        show_coords = op[c4d.GEOTIFF_SHOW_COORDS] if op[c4d.GEOTIFF_SHOW_COORDS] is not None else False
        step = op[c4d.GEOMESH_STEP] if op[c4d.GEOMESH_STEP] is not None else 10
        global_scale = op[c4d.GEOMESH_GLOBAL_SCALE] if op[c4d.GEOMESH_GLOBAL_SCALE] is not None else 0.01
        measure_scale = op[c4d.GEOMESH_MEASURE_SCALE] if op[c4d.GEOMESH_MEASURE_SCALE] is not None else 2
        vertical_align = op[c4d.GEOMESH_VERTICAL_ALIGN] if op[c4d.GEOMESH_VERTICAL_ALIGN] is not None else c4d.GEOMESH_VERTICAL_MIN_VALUE
        solid_base = op[c4d.GEOMESH_SOLID_BASE] if op[c4d.GEOMESH_SOLID_BASE] is not None else True
        base_thickness = op[c4d.GEOMESH_BASE_THICKNESS] if op[c4d.GEOMESH_BASE_THICKNESS] is not None else 100
        
        m = bh.GetMg()
        bd.SetMatrix_Matrix(op, m)
        
        
        # If the file is not None and the raster is not loaded, then load the raster
        if self.raster is not None and self.arr is not None and self.arr_red is not None:
            
            raster = self.raster
            arr = self.arr
            arr_red = self.arr_red
            
            # Size of the original raster "arr"
            x_segments, z_segments = arr.shape
            
            # Bounding box of the raster
            bottom = raster.bounds.bottom
            top = raster.bounds.top
            
            # Center of the raster
            lat_center = (top + bottom)/2

            # Depending on the units, the measure function will be different, by default only meters and degrees are supported
            # Measured data (for example, elevation) will be mapped to proportional degrees based on the latitude of the center of the raster
            if units == 'degrees':
                measureMap = lambda x, lat=0, lon=0: abs(x * 360 / (2 * pi * EARTH_RADIUS * cos(radians(lat))))
            else:
                measureMap = lambda x, lat=0, lon=0: x
            
            center_x, center_z = raster.transform*(int(x_segments/2), int(z_segments/2))
            
            # To find minimum, maximum and range without considering no data values and reduced array
            arr_red = np.where(arr_red == no_data, np.nan, arr_red)     # Replace no_data with nan
            
            center_m = measureMap((np.nanmax(arr_red) - np.nanmin(arr_red))/2, lat_center)    # Range of the array
            min_m = measureMap(np.nanmin(arr_red),lat_center)           # Minimum value of the array
            max_m = measureMap(np.nanmax(arr_red),lat_center)           # Maximum value of the array
            range_m = max_m - min_m                                     # Range of the array
            max_m_scale = range_m*global_scale*measure_scale           # Maximum value of the array scaled
            
            
            # If there is no value for "no data" then the minimum value of the array - 1000 is used
            if no_data == None:
                no_data = int(min_m - 1000)

                
            if vertical_align == c4d.GEOMESH_VERTICAL_MIN_VALUE:
                center_m = min_m
            elif vertical_align == c4d.GEOMESH_VERTICAL_MAX_VALUE:
                center_m = max_m
            elif vertical_align == c4d.GEOMESH_VERTICAL_MEAN_VALUE:
                center_m = measureMap(np.nanmean(arr_red),lat_center) 
            elif vertical_align == c4d.GEOMESH_VERTICAL_CENTER_VALUE:
                center_m = (min_m + max_m)/2
            elif vertical_align == c4d.GEOMESH_VERTICAL_SOLID_BASE:
                center_m = min_m - base_thickness/(measure_scale*global_scale)
                
            
            # Draw the bounding box of the raster, get real coordinates and transformed with scale
            b1x, b1z = raster.transform*(0, 0)
            b1xt, b1zt = (b1x - center_x)*global_scale, (b1z - center_z)*global_scale
            
            b2x, b2z = raster.transform*(x_segments, 0)
            b2xt, b2zt = (b2x - center_x)*global_scale, (b2z - center_z)*global_scale
            
            b3x, b3z = raster.transform*(0, z_segments)
            b3xt, b3zt = (b3x - center_x)*global_scale, (b3z - center_z)*global_scale
            
            b4x, b4z = raster.transform*(x_segments, z_segments)
            b4xt, b4zt = (b4x - center_x)*global_scale, (b4z - center_z)*global_scale
            
            vertical_base = (min_m - center_m)*global_scale*measure_scale
            
            # Real bbox
            p1 = c4d.Vector(b1z, vertical_base, b1x)
            p2 = c4d.Vector(b2z, vertical_base, b2x)
            p3 = c4d.Vector(b3z, vertical_base, b3x)
            p4 = c4d.Vector(b4z, vertical_base, b4x)
            
            # Scaled bbox
            p1t = c4d.Vector(b1zt, vertical_base, b1xt)
            p2t = c4d.Vector(b2zt, vertical_base, b2xt)
            p3t = c4d.Vector(b3zt, vertical_base, b3xt)
            p4t = c4d.Vector(b4zt, vertical_base, b4xt)
            
            # Draw base
            bd.SetPen(c4d.Vector(1, 0, 1))
            bd.DrawLine(p1t, p2t, c4d.NOCLIP_Z)
            bd.DrawLine(p2t, p4t, c4d.NOCLIP_Z)   
            bd.DrawLine(p4t, p3t, c4d.NOCLIP_Z)
            bd.DrawLine(p3t, p1t, c4d.NOCLIP_Z)
            
            
            # Draw sampled points
            step_mini = 15*step
            virtualPoints = []
            pointColors = []
            for i in range(0, arr.shape[0], step_mini):
                for j in range(0, arr.shape[1], step_mini):
                    
                    xr, zr = raster.transform*(i, j)
                    yr = measureMap(arr[i,j], zr, xr)
                    
                    # Calculate x, y, z scaled and centered
                    x, z = (xr - center_x)*global_scale, (zr - center_z)*global_scale
                    y = (yr - center_m)*global_scale*measure_scale

                    if yr != no_data:    
                        virtualPoints.append(c4d.Vector(z, y, x))
                        pointColors += [y]

            pointColors = [c - min(pointColors) for c in pointColors]
            pointColors = [c/max(pointColors) for c in pointColors]
            
            pointColorsRGB = []
            for c in pointColors:
                pointColorsRGB += [1, c, 1-c]
            
            # bd.SetPen(c4d.Vector(1, 1, 1))
            bd.SetPointSize(3)
            bd.DrawPoints(virtualPoints, pointColorsRGB, 3)
            
            
            # Draw bbox displaced to the top (based on max height) ---------------------
            # Scaled bbox top
            p1t_top = c4d.Vector(b1zt, vertical_base + max_m_scale, b1xt)
            p2t_top = c4d.Vector(b2zt, vertical_base + max_m_scale, b2xt)
            p3t_top = c4d.Vector(b3zt, vertical_base + max_m_scale, b3xt)
            p4t_top = c4d.Vector(b4zt, vertical_base + max_m_scale, b4xt)
            
            # Draw base level
            bd.SetPen(c4d.Vector(1, 1, 1))
            bd.DrawLine(p1t_top, p2t_top, c4d.NOCLIP_Z)
            bd.DrawLine(p2t_top, p4t_top, c4d.NOCLIP_Z)
            bd.DrawLine(p4t_top, p3t_top, c4d.NOCLIP_Z)
            bd.DrawLine(p3t_top, p1t_top, c4d.NOCLIP_Z)
            
            # Show the line of the scaled range -------------------------------------------
            bd.SetPen(c4d.Vector(1, 1, 1))
            bd.DrawLine(p1t, p1t_top, c4d.NOCLIP_Z)
                  
            # Draw bbox displaced to the bottom (if solid base is true) ---------------
            
            if solid_base:

                # Scaled bbox bottom
                p1t_bottom = c4d.Vector(b1zt, vertical_base - base_thickness, b1xt)
                p2t_bottom = c4d.Vector(b2zt, vertical_base - base_thickness, b2xt)
                p3t_bottom = c4d.Vector(b3zt, vertical_base - base_thickness, b3xt)
                p4t_bottom = c4d.Vector(b4zt, vertical_base - base_thickness, b4xt)
                
                # Draw base level
                bd.SetPen(c4d.Vector(1, 0, 1))
                bd.DrawLine(p1t_bottom, p2t_bottom, c4d.NOCLIP_Z)
                bd.DrawLine(p2t_bottom, p4t_bottom, c4d.NOCLIP_Z)
                bd.DrawLine(p4t_bottom, p3t_bottom, c4d.NOCLIP_Z)
                bd.DrawLine(p3t_bottom, p1t_bottom, c4d.NOCLIP_Z)
                
                # Draw the sides of the base
                bd.DrawLine(p1t, p1t_bottom, c4d.NOCLIP_Z)
                bd.DrawLine(p2t, p2t_bottom, c4d.NOCLIP_Z)
                bd.DrawLine(p3t, p3t_bottom, c4d.NOCLIP_Z)
                bd.DrawLine(p4t, p4t_bottom, c4d.NOCLIP_Z)
            
            
            # Draw coord of the bbox corners -------------------------------------------
            if show_coords:
                
                # Project scaled bbox to screen
                p1_screen = bd.WS(p1t)
                p2_screen = bd.WS(p2t)
                p3_screen = bd.WS(p3t)
                p4_screen = bd.WS(p4t)
                
                coords_text = lambda s: str(s.z) + ', ' + str(s.x)
                
                
                # Draw text of real values of bbox
                bd.DrawHUDText(int(p1_screen.x), int(p1_screen.y), coords_text(p1))
                bd.DrawHUDText(int(p2_screen.x), int(p2_screen.y), coords_text(p2))
                bd.DrawHUDText(int(p3_screen.x), int(p3_screen.y), coords_text(p3))
                bd.DrawHUDText(int(p4_screen.x), int(p4_screen.y), coords_text(p4))

            # Show measure scaled value
            p1t_top_screen = bd.WS(p1t_top)
            bd.DrawHUDText(int(p1t_top_screen.x), int(p1t_top_screen.y), 'x'+str(round(measure_scale,4)))
      

        return c4d.DRAWRESULT_OK
    
    
    
    # Function to read the attributes from the interface and create the surface object
    def GetVirtualObjects(self, op, hh, restart_gui=False):
        
        file = op[c4d.GEOTIFF_FILE] if op[c4d.GEOTIFF_FILE] is not None else DEFAULT_FILE
        units = op[c4d.GEOTIFF_UNITS] if op[c4d.GEOTIFF_UNITS] is not None else c4d.GEOTIFF_UNITS_METERS
        layer = op[c4d.GEOTIFF_LAYER] if op[c4d.GEOTIFF_LAYER] is not None else 1
        shape = op[c4d.GEOMESH_TYPE] if op[c4d.GEOMESH_TYPE] is not None else c4d.GEOMESH_TYPE_TRN
        step = op[c4d.GEOMESH_STEP] if op[c4d.GEOMESH_STEP] is not None else 10
        global_scale = op[c4d.GEOMESH_GLOBAL_SCALE] if op[c4d.GEOMESH_GLOBAL_SCALE] is not None else 0.01
        measure_scale = op[c4d.GEOMESH_MEASURE_SCALE] if op[c4d.GEOMESH_MEASURE_SCALE] is not None else 2
        vertical_align = op[c4d.GEOMESH_VERTICAL_ALIGN] if op[c4d.GEOMESH_VERTICAL_ALIGN] is not None else c4d.GEOMESH_VERTICAL_MIN_VALUE
        solid_base = op[c4d.GEOMESH_SOLID_BASE] if op[c4d.GEOMESH_SOLID_BASE] is not None else True
        base_thickness = op[c4d.GEOMESH_BASE_THICKNESS] if op[c4d.GEOMESH_BASE_THICKNESS] is not None else 100
        uvw = op[c4d.GEOMESH_UVW] if op[c4d.GEOMESH_UVW] is not None else c4d.GEOMESH_UVW_HEIGHT
        tin_seed = op[c4d.GEOMESH_TYPE_TIN_SEED] if op[c4d.GEOMESH_TYPE_TIN_SEED] is not None else 0
        gaussian_sampling_threshold = op[c4d.GEOMESH_GAUSSIAN_SAMPLING_THRESHOLD] if op[c4d.GEOMESH_GAUSSIAN_SAMPLING_THRESHOLD] is not None else 0.98
        uniform_sampling_threshold = op[c4d.GEOMESH_UNIFORM_SAMPLING_THRESHOLD] if op[c4d.GEOMESH_UNIFORM_SAMPLING_THRESHOLD] is not None else 0.5
        gradient_threshold = op[c4d.GEOMESH_GRADIENT_THRESHOLD] if op[c4d.GEOMESH_GRADIENT_THRESHOLD] is not None else 2
        generate = op[c4d.GEOMESH_GENERATE] if op[c4d.GEOMESH_GENERATE] is not None else False
        no_data = op[c4d.GEOTIFF_NO_DATA] if op[c4d.GEOTIFF_NO_DATA] is not None else int(np.nanmin(self.arr_red) - 1000)

        # If the file is not None and the raster is not loaded, then load the raster
        if self.arr is None and file is not None:
            
            # Open and re sample the raster
            self.raster = rio.open(r'{}'.format(op[c4d.GEOTIFF_FILE]))
            self.arr = self.raster.read(layer)
            self.arr_red = self.arr[::step, ::step]
            no_data = self.raster.nodatavals[layer - 1]              
            self.arr_red = np.where(self.arr_red == no_data, np.nan, self.arr_red)
            
            # Update no data value
            no_data = op[c4d.GEOTIFF_NO_DATA] if op[c4d.GEOTIFF_NO_DATA] is not None else int(np.nanmin(self.arr_red) - 1000)
            op[c4d.GEOTIFF_NO_DATA] = no_data
            
        
        # If the "generate" checkbox is seleted, then create the mesh
        if generate:
            
            raster = self.raster
            arr = self.arr
            arr_red = self.arr_red
            
            if raster is None or arr is None or arr_red is None:
                return None
            
            # Size of the original raster "arr"
            x_segments, z_segments = arr.shape
            max_size = max(x_segments, z_segments)
            min_size = min(x_segments, z_segments)
            
            # Bounding box of the raster
            bottom = raster.bounds.bottom
            top = raster.bounds.top
            left = raster.bounds.left
            right = raster.bounds.right
            
            # Center of the raster
            lat_center = (top + bottom)/2
            lon_center = (left + right)/2

            
            # Depending on the units, the measure function will be different, by default only meters and degrees are supported
            # Measured data (for example, elevation) will be mapped to proportional degrees based on the latitude of the center of the raster
            if units == 'degrees':
                measureMap = lambda x, lat=0, lon=0: abs(x * 360 / (2 * pi * EARTH_RADIUS * cos(radians(lat))))
            else:
                measureMap = lambda x, lat=0, lon=0: x
            
            center_x, center_z = raster.transform*(int(x_segments/2), int(z_segments/2))
            
            center_m = measureMap((np.nanmax(arr_red) - np.nanmin(arr_red))/2, lat_center)    # Range of the array
            min_m = measureMap(np.nanmin(arr_red),lat_center)           # Minimum value of the array
            max_m = measureMap(np.nanmax(arr_red),lat_center)           # Maximum value of the array
            range_m = max_m - min_m                                     # Range of the array
            max_m_scale = range_m*global_scale*measure_scale           # Maximum value of the array scaled
            all_no_data = np.sum(np.isnan(arr_red))                     # All values are no data
            
            # If there is no value for "no data" then the minimum value of the array - 1000 is used
            if no_data == None:
                no_data = int(min_m - 1000)
                op[c4d.GEOTIFF_NO_DATA] = no_data
                
            
            # Replace back the no data values
            arr_red = np.nan_to_num(arr_red, nan=no_data)

            
            centered = True
            if not centered:    
                op.SetRelPos(c4d.Vector(center_x, center_m, center_z))
            
            # Number of points 
            NX, NZ = arr_red.shape
            
            # Compute total number of points and polygons
            total_points = 2*NX*NZ
            surface_polys = (NX-1)*(NZ-1)
            
            # Border division identifier for id's
            bd_points = NX*NZ
            bd_poly = (NX-1)*(NZ-1)
            
            # Create the object
            if solid_base:
                obj = c4d.PolygonObject(total_points + bd_points, 
                                        2*surface_polys + all_no_data*4 + 2*(NX + 1) + 2*(NZ + 1))
            else:
                obj = c4d.PolygonObject(total_points, surface_polys)
                
            # Set the offset based on the vertical alignment
            if vertical_align == c4d.GEOMESH_VERTICAL_MIN_VALUE:
                center_m = min_m
            elif vertical_align == c4d.GEOMESH_VERTICAL_MAX_VALUE:
                center_m = max_m
            elif vertical_align == c4d.GEOMESH_VERTICAL_MEAN_VALUE:
                center_m = measureMap(np.nanmean(arr_red),lat_center) 
            elif vertical_align == c4d.GEOMESH_VERTICAL_CENTER_VALUE:
                center_m = (min_m + max_m)/2
            elif vertical_align == c4d.GEOMESH_VERTICAL_SOLID_BASE:
                center_m = base_thickness*(measure_scale*global_scale)
            
            
            # Call the function to create the mesh
            if shape == c4d.GEOMESH_TYPE_TRN:
                
                # Create the Triangular Regular Network
                obj, uvw_tag = TRN.createTRN(obj, raster, measureMap, solid_base, uvw, step, global_scale, 
                                    measure_scale, max_m_scale, center_x, center_z, center_m, base_thickness, no_data, min_m, arr,
                                    NX, NZ, all_no_data, surface_polys, bd_points, bd_poly)
            
            elif shape == c4d.GEOMESH_TYPE_TIN:
                # Create the Triangular Irregular Network
                obj, uvw_tag = TIN.createTIN(obj, raster,step, measureMap, arr, arr_red, no_data, center_x, center_z, center_m,
                                            global_scale, measure_scale, max_m_scale, solid_base, base_thickness, uvw, x_segments, z_segments,
                                            min_m, range_m, NX, NZ, bd_points, bd_poly, min_size, tin_seed,
                                            gradient_threshold, uniform_sampling_threshold, gaussian_sampling_threshold)
                
            op.KillTag(c4d.Tuvw)  
            op.KillTag(c4d.Tedgeselection)
            op.KillTag(c4d.Tpolygonselection)  
            
            polygon_selection = obj.GetTag(c4d.Tpolygonselection)
            edge_selection = obj.GetTag(c4d.Tedgeselection)
            
            if polygon_selection is not None:
                op.InsertTag(polygon_selection)
            if edge_selection is not None:
                op.InsertTag(edge_selection)    
            if uvw_tag is not None:
                op.InsertTag(uvw_tag)
            
            # Copy tags from the generator to the object
            op.CopyTagsTo(obj, c4d.NOTOK,c4d.NOTOK,c4d.NOTOK)
            
            # Update object
            obj.Message(c4d.MSG_UPDATE)
        

            return obj
        
        else:
            return None
        
    def setDefault(self, node):
        
        node.SetParameter(c4d.GEOTIFF_FILE, DEFAULT_FILE, c4d.DESCFLAGS_SET_0)
            
        node[c4d.GEOTIFF_FILE] = DEFAULT_FILE
        node[c4d.GEOTIFF_UNITS] = c4d.GEOTIFF_UNITS_METERS
        node[c4d.GEOTIFF_LAYER] = 1
        node[c4d.GEOTIFF_NO_DATA] = 0
        node[c4d.GEOTIFF_SHOW_COORDS] = False
        node[c4d.GEOMESH_TYPE] = c4d.GEOMESH_TYPE_TRN
        node[c4d.GEOMESH_STEP] = 10
        node[c4d.GEOMESH_GLOBAL_SCALE] = 0.01
        node[c4d.GEOMESH_MEASURE_SCALE] = 2
        node[c4d.GEOMESH_VERTICAL_ALIGN] = c4d.GEOMESH_VERTICAL_MIN_VALUE
        node[c4d.GEOMESH_SOLID_BASE] = True
        node[c4d.GEOMESH_BASE_THICKNESS] = 100
        node[c4d.GEOMESH_UVW] = c4d.GEOMESH_UVW_HEIGHT
        node[c4d.GEOMESH_TYPE_TIN_SEED] = 0
        node[c4d.GEOMESH_GAUSSIAN_SAMPLING_THRESHOLD] = 0.98
        node[c4d.GEOMESH_UNIFORM_SAMPLING_THRESHOLD] = 0.5
        node[c4d.GEOMESH_GRADIENT_THRESHOLD] = 2
        node[c4d.GEOMESH_GENERATE] = False
        
        if node[c4d.GEOTIFF_FILE] != None:
            self.raster = rio.open(r'{}'.format(node[c4d.GEOTIFF_FILE]))
            self.arr = self.raster.read(node[c4d.GEOTIFF_LAYER])
            
            max_size = max(self.arr.shape)
            step = max(1, int(max_size/250))
            
            self.arr_red = self.arr[::step, ::step]
            node[c4d.GEOMESH_STEP] = step
            
            no_data = self.raster.nodatavals[node[c4d.GEOTIFF_LAYER] - 1]              
            self.arr_red = np.where(self.arr_red == no_data, np.nan, self.arr_red)
            
            if no_data == None:
                no_data = int(np.nanmin(self.arr_red) - 1000)
            node[c4d.GEOTIFF_NO_DATA] = no_data

        
            
        else:
            self.raster = None
            self.arr = None
            self.arr_red = None
        
        
        
    

    # Function to handle the changes in the button
    def Message(self, node, type, data):
        ### Called by Cinema 4D to retrieve the current surface type.
        
        # If one element in the interface has changed
        if type == c4d.MSG_DESCRIPTION_COMMAND:
            
            # If the reset button has been pressed
            if data['id'][0].id == c4d.SURF_RESET:
                self.surf_obj = copy.deepcopy(self.surf_obj_init)              
                self.SetFromSurface(node, self.surf_obj)
                self.GetVirtualObjects(node, None, restart_gui=True)
                # node.Message(c4d.MSG_UPDATE)
        

        elif type == c4d.MSG_DESCRIPTION_POSTSETPARAMETER:

            # If the surface file has changed
            if data['descid'][0].id == c4d.GEOTIFF_FILE:
                
                if node[c4d.GEOTIFF_FILE] != None:
                    
                    try:
                        self.raster = rio.open(r'{}'.format(node[c4d.GEOTIFF_FILE]))
                        self.arr = self.raster.read(node[c4d.GEOTIFF_LAYER])
                        
                        max_size = max(self.arr.shape)
                        step = max(1, int(max_size/250))
                        
                        self.arr_red = self.arr[::step, ::step]
                        node[c4d.GEOMESH_STEP] = step
                        
                        no_data = self.raster.nodatavals[node[c4d.GEOTIFF_LAYER] - 1]              
                        self.arr_red = np.where(self.arr_red == no_data, np.nan, self.arr_red)
                        
                        if no_data == None:
                            no_data = int(np.nanmin(self.arr_red) - 1000)
                        
                        node[c4d.GEOTIFF_NO_DATA] = no_data
                        
                            
                    except:
                        # Not implemented, GUI DIALOG ERROR
                        pass
                        
                pass
            
            # Update array reduced
            if data['descid'][0].id == c4d.GEOMESH_STEP:
                if self.arr is not None and node[c4d.GEOMESH_STEP] != None:
                    step = node[c4d.GEOMESH_STEP]
                    
                    self.raster = rio.open(r'{}'.format(node[c4d.GEOTIFF_FILE]))
                    self.arr = self.raster.read(node[c4d.GEOTIFF_LAYER])
                    
                    self.arr_red = self.arr[::step, ::step]
                    node[c4d.GEOMESH_STEP] = step
                    
                    no_data = self.raster.nodatavals[node[c4d.GEOTIFF_LAYER] - 1]              
                    self.arr_red = np.where(self.arr_red == no_data, np.nan, self.arr_red)
                    
                    if no_data == None:
                        no_data = int(np.nanmin(self.arr_red) - 1000)
                    
                    node[c4d.GEOTIFF_NO_DATA] = no_data
        
        # Mesage received just before the interface is displayed
        
        elif type == c4d.MSG_MENUPREPARE:
            
            # Insert a phong tag by default
            phongTag = c4d.BaseTag(c4d.Tphong)
            phongTag[c4d.PHONGTAG_PHONG_ANGLELIMIT] = True
            phongTag[c4d.PHONGTAG_PHONG_ANGLE] = radians(20)
            node.InsertTag(phongTag)
            
            self.InitAttr(node, str, c4d.GEOTIFF_FILE)
            self.InitAttr(node, int, c4d.GEOTIFF_UNITS)
            self.InitAttr(node, int, c4d.GEOTIFF_LAYER)
            self.InitAttr(node, float, c4d.GEOTIFF_NO_DATA)
            self.InitAttr(node, bool, c4d.GEOTIFF_SHOW_COORDS)
            self.InitAttr(node, int, c4d.GEOMESH_TYPE)
            self.InitAttr(node, float, c4d.GEOMESH_STEP)
            self.InitAttr(node, float, c4d.GEOMESH_GLOBAL_SCALE)
            self.InitAttr(node, float, c4d.GEOMESH_MEASURE_SCALE)
            self.InitAttr(node, int, c4d.GEOMESH_VERTICAL_ALIGN)
            self.InitAttr(node, bool, c4d.GEOMESH_SOLID_BASE)
            self.InitAttr(node, float, c4d.GEOMESH_BASE_THICKNESS)
            self.InitAttr(node, int, c4d.GEOMESH_UVW)
            self.InitAttr(node, int, c4d.GEOMESH_TYPE_TIN_SEED)
            self.InitAttr(node, float, c4d.GEOMESH_GAUSSIAN_SAMPLING_THRESHOLD)
            self.InitAttr(node, float, c4d.GEOMESH_GRADIENT_THRESHOLD)
            self.InitAttr(node, bool, c4d.GEOMESH_GENERATE)
            
            self.setDefault(node)

        return True



    # Function to handle the changes in the interface
    def GetDEnabling(self, node, id, t_data, flags, itemdesc):
        ### "Called  by Cinema 4D to decide which parameters should be enabled or disabled (ghosted).

        # Do not allow to generate the solid base if the mesh is a TIN (not implemented yet)
        triangulation = node[c4d.GEOMESH_TYPE] if node[c4d.GEOMESH_TYPE] is not None else c4d.GEOMESH_TYPE_TRN
        
        if id[0].id == c4d.GEOMESH_SOLID_BASE and triangulation == c4d.GEOMESH_TYPE_TIN:
            node[c4d.GEOMESH_SOLID_BASE] = False
            return False
        
        # No acces to TIN properties if the mesh is a TRN
        tin_props = set([c4d.GEOMESH_TYPE_TIN_SEED, c4d.GEOMESH_GAUSSIAN_SAMPLING_THRESHOLD, c4d.GEOMESH_UNIFORM_SAMPLING_THRESHOLD, c4d.GEOMESH_GRADIENT_THRESHOLD])
        if id[0].id  in tin_props and triangulation == c4d.GEOMESH_TYPE_TRN:
            return False

        return True
    


    
if __name__ == "__main__":
    
    # Load the plugin icon
    icon_absolute_path = os.path.join(os.path.dirname(__file__), 'res/icons', 'icon.png')
    plugin_icon = bitmaps.BaseBitmap()
    plugin_icon.InitWith(icon_absolute_path)

    # Register the plugin
    plugins.RegisterObjectPlugin(
        id = 1062710,
        str = 'GeoTIFF2Mesh',
        g =  GeoTIFF2Mesh,
        description = 'Ogeomesh',
        info = c4d.OBJECT_GENERATOR,
        icon = plugin_icon
    )
