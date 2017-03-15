#In this file we define map message that will be published
Header header          #standard header information
float64 x_size         #in meters
float64 y_size         #in meters
float64 z_size         #in meters
float64 x_cen          #cneter of the map
float64 y_cen          #cneter of the map
float64 z_cen          #cneter of the map
float64 x_cell_size    #cell size
float64 y_cell_size    #cell size
float64 z_cell_size    #cell size

NDTCellMsg[] cells     #cells
