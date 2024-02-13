CONTAINER Ogeomesh {
    INCLUDE Obase;
    NAME Ogeomesh;

    GROUP ID_GEOTIFF_DATA {

        FILENAME GEOTIFF_FILE { }
        GROUP {
            DEFAULT 1;
            LAYOUTGROUP;
            COLUMNS 2;

            GROUP {
            LONG GEOTIFF_UNITS { 
                CYCLE { 
                    GEOTIFF_UNITS_METERS; 
                    GEOTIFF_UNITS_DEGREES; } }
            LONG GEOTIFF_LAYER { MIN 1; }
            }


            GROUP { 
            REAL GEOTIFF_NO_DATA { MIN -99999999; MAX 99999999; STEP 0.1; } 
            BOOL GEOTIFF_SHOW_COORDS { }
            }
        }


        GROUP ID_GEOMESH_PROPERTIES {
            DEFAULT 1;

            LONG GEOMESH_TYPE { 
                CYCLE { 
                GEOMESH_TYPE_TRN; 
                GEOMESH_TYPE_TIN; 
                } 
            }
            LONG GEOMESH_STEP { MIN 1; MAX 10000000; STEP 1; MINSLIDER 0; MAXSLIDER 300; CUSTOMGUI LONGSLIDER; }
            REAL GEOMESH_GLOBAL_SCALE { MIN -999999; MAX 999999; STEP 0.1; MINSLIDER 0; MAXSLIDER 100; CUSTOMGUI REALSLIDER; }
            REAL GEOMESH_MEASURE_SCALE { MIN -999999; MAX 999999; STEP 0.1; MINSLIDER 0; MAXSLIDER 20; CUSTOMGUI REALSLIDER; }
            
            LONG GEOMESH_VERTICAL_ALIGN { 
                CYCLE { 
                    GEOMESH_VERTICAL_MIN_VALUE;
                    GEOMESH_VERTICAL_MAX_VALUE;
                    GEOMESH_VERTICAL_MEAN_VALUE;
                    GEOMESH_VERTICAL_CENTER_VALUE;
                    } 
                }

            BOOL GEOMESH_SOLID_BASE { }

            REAL GEOMESH_BASE_THICKNESS { MIN -999999; MAX 999999; STEP 1; MINSLIDER 0; MAXSLIDER 200; CUSTOMGUI REALSLIDER; }

            LONG GEOMESH_UVW { 
                CYCLE { 
                    GEOMESH_UVW_TOP_VIEW;
                    GEOMESH_UVW_HEIGHT;
                    } 
                }
            
            GROUP ID_GEOMESH_TIN_PROPERTIES {
                DEFAULT 0;
                LONG GEOMESH_TYPE_TIN_SEED { MIN 1; MAX 999999999; STEP 1; }
                REAL GEOMESH_GRADIENT_THRESHOLD { MIN 0.001; MAX 100; STEP 0.01; }
                REAL GEOMESH_UNIFORM_SAMPLING_THRESHOLD { MIN 0.001; MAX 1; STEP 0.01; }
                REAL GEOMESH_GAUSSIAN_SAMPLING_THRESHOLD { MIN 0.001; MAX 1; STEP 0.01; }
                
            }

            BOOL GEOMESH_GENERATE { }

            
        }
    }
}