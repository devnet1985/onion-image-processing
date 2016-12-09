#1)obtain coordinate of corners
#2)obtain nadir images for each left up and right down corners, generate coordinate and image pair 
#3)extract pixels according to coordinates and image 
import scipy.io as sio
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import numpy as np
import math

from os import listdir
from osgeo import gdal
from osgeo import osr

def main():
    global cord_mag,reg_mag,mag,tr_map
    mag=10
    #mag=80
    cord_mag=[10,0]
    tr_map=np.zeros([36,4])
    # rowst,rowend, colst, colend
    mark_coordinate_dir="I:\ImagerySets\USDA mission\MarkersCoordinatesWGS-84.xlsx"
    mark_cord=np.array(pd.read_excel(mark_coordinate_dir))
    #region_left_up=np.array([mark_cord[0:144:4,2],mark_cord[0:144:4,1]])
    region_left_up=np.array([mark_cord[144:288:4,2],mark_cord[144:288:4,1]])
    print(mark_cord[144:288:4,0])
    
    region_left_up=np.transpose(region_left_up)#-dlon*np.array([0,0.000001])
    #region_right_down=np.array([mark_cord[2:144:4,2],mark_cord[2:144:4,1]])
    region_right_down=np.array([mark_cord[146:288:4,2],mark_cord[146:288:4,1]])
    print(mark_cord[146:288:4,0])
    region_right_down=np.transpose(region_right_down)#+np.ones([36,1])*np.array([0,0.000001])*2
    #--1--2---
    #--------- 
    #--4--3---
    # pixelPairs[][0] is hor ; pixelPairs[][1] is ver
    # [lat, lon]    
    
    
 
    img_list=['onion-rgb-nir-tir-06-02-16.tif']
    #date_list=['04-07-16','05-10-16','06-02-16','06-14-16','06-28-16','07-12-16','07-26-16']
    orthophoto_dir_head="I:\\ImagerySets\\USDA mission\\onion-2016-data-blue region\\test\\06-02-16\\"
    for img in img_list:
    #orthophoto_dir="I:\\ImagerySets\\USDA mission\\05-10-16\\onion-nir-orthophotos-05-10-16\\"
    #orthophoto_dir="I:\\ImagerySets\\USDA mission\\04-07-16\\onion-TIR-04-07-16-orthophotos\\"
        orthophoto_dir=orthophoto_dir_head+img
        #region_left_up,region_right_down,nadir_image=obtain_mark_image_pair(mark_lat_lon,orthophoto_dir)
        print(orthophoto_dir)
        region_array=obtain_region_array(region_left_up,region_right_down,orthophoto_dir)
        matrix_name="I:\\ImagerySets\\USDA mission\\onion-2016-data-blue region\\"+img+"-nir.mat"
        sio.savemat(matrix_name,{'img_rgb_nir_tir':region_array})
        matrix_name="I:\\ImagerySets\\USDA mission\\onion-2016-data-blue region\\"+img+"-tr_map.mat"
        sio.savemat(matrix_name,{'tr_map':tr_map})
        plt.imshow(region_array[:,:,0],cmap='gray')
        #plt.imshow(np.transpose(region_array[:,:,1]),cmap='gray')
        plt.show()
    
        """
        matrix_name="I:\\ImagerySets\\USDA mission\\onion-2016-data-blue region\\"+img+"-g.mat"
        sio.savemat(matrix_name,{'img_g':region_array[:,:,1]})
        matrix_name="I:\\ImagerySets\\USDA mission\\onion-2016-data-blue region\\"+img+"-b.mat"
        sio.savemat(matrix_name,{'img_b':region_array[:,:,2]})
        """
        

def obtain_region_array(region_left_up,region_right_down,img_addr):



     
    lat_dif=(region_left_up[:,0]-region_right_down[:,0])
    # latitude of left up is larger than that of right down
    lon_dif=-(region_left_up[:,1]-region_right_down[:,1])
    tune_width=2
    max_gps=np.array([lat_dif.min(),lon_dif.max()*(tune_width+1)/2])
    max_gps.shape=[1,2]
    
    max_gps_array=np.dot(np.multiply(np.ones([36,1]),max_gps),np.array([[1,0],[0,-1]]))
  
    dlon=np.ones([36,1])*0
    # positive move the window left, negative move the window right
    """
    #for 04-07-16
    dlon[0:6]=np.array([-1,0,-1,-1,-1,-1]).reshape((6,1))
    dlon[6:12]=np.array([-1,-1,-2,-1,-1,-1]).reshape((6,1)) 
    dlon[12:18]=np.array([0,0,0,-1,-1,-1]).reshape((6,1))
    dlon[18:24]=np.array([-1,-1,-1,-1,-1,-1]).reshape((6,1))
    dlon[24:30]=np.array([-1,-1,-1,-1,-1,-2]).reshape((6,1))
    dlon[30:36]=np.array([0,0,-1,-1,-1,-1]).reshape((6,1))
    """
    
    #for 06-02-16
    dlon[0:6]=np.array([-1,0,-1,-1,0,-1]).reshape((6,1))
    dlon[6:12]=np.array([0,0,0,0,0,0]).reshape((6,1)) 
    dlon[12:18]=np.array([0,0,0,-1,-1,-1]).reshape((6,1))
    dlon[18:24]=np.array([-1,-1,0,0,0,-1]).reshape((6,1))
    dlon[24:30]=np.array([-1,-1,-1,-1,-1,-2]).reshape((6,1))
    dlon[30:36]=np.array([0,0,0,0,0,-1]).reshape((6,1))
    
    #dlon[2]=
    #dlon[4]=4
    region_left_up_eq=region_left_up-np.dot(np.ones([36,2]),np.array([[0,0],[0,1]]))*lon_dif.max()*(tune_width-1)/2-dlon*np.array([0,0.000001])
    region_right_down_eq=region_left_up_eq-max_gps_array

    # convert gps to image coordinate
    pixelPair_left_up=latLonToPixel(img_addr,region_left_up_eq)
    # right lat+margin(0.000005), left lon-margin(0.000005)
    pixelPair_right_down=latLonToPixel(img_addr,region_right_down_eq)
    region_size=pixelPair_right_down-pixelPair_left_up

    region_size[:,0]=-(pixelPair_left_up[:,0]-pixelPair_right_down[:,0])
    # latitude of left up is larger than that of right down
    region_size[:,1]=-(pixelPair_left_up[:,1]-pixelPair_right_down[:,1])
    max_region=np.max(region_size,axis=0)+[mag*2,mag*2] 

    
    # calculate size of panel
    max_x=(max_region[0]+mag)*12+mag
    max_y=(max_region[1]+mag)*3+mag
    panel =np.zeros((max_y,max_x,5))

    # calculate position for each plot in the panel
    region_loc_xy=np.zeros((3,12,2))
    for j in range(3):
        for i in range(12):
            if i==0:
                xi=mag
            else:
                xi=mag+(mag+max_region[0])*i
            if j==0:
                yj=mag
            else:
                yj=mag+(mag+max_region[1])*j

            region_loc_xy[j,i,:]=[xi,yj]

            
        
    

    # place the pixel in the panel
    ds=gdal.Open(img_addr)
    #myarray2 = np.array(ds.GetRasterBand(2).ReadAsArray())
    #myarray3 = np.array(ds.GetRasterBand(3).ReadAsArray())
    for i in range(36):
        print(i)
        row_lu=math.floor(max(pixelPair_left_up[i,1],0))
        col_lu=math.floor(max(pixelPair_left_up[i,0],0))
        row_rd=math.floor(max(pixelPair_right_down[i,1],0))
        col_rd=math.floor(max(pixelPair_right_down[i,0],0))

        if i<12:
            yj=0
            xi=i
        elif i<24:
            yj=1
            xi=(i-12)
        else:
            yj=2
            xi=(i-24)
      
        [reg_col,reg_row]=region_loc_xy[yj,xi,:]
        print(reg_col,reg_row)
        #col_lu,row_lu,col_rd,row_rd
        tr_map[i,:]=np.array([reg_row,reg_row+row_rd-row_lu,reg_col,reg_col+col_rd-col_lu])
        panel[reg_row:reg_row+row_rd-row_lu,reg_col:reg_col+col_rd-col_lu,0]=ds.GetRasterBand(1).ReadAsArray(col_lu,row_lu,(col_rd-col_lu),(row_rd-row_lu))
        print(row_rd-row_lu,col_rd-col_lu)
        panel[reg_row:reg_row+row_rd-row_lu,reg_col:reg_col+col_rd-col_lu,1]=ds.GetRasterBand(2).ReadAsArray(col_lu,row_lu,(col_rd-col_lu),(row_rd-row_lu))
        panel[reg_row:reg_row+row_rd-row_lu,reg_col:reg_col+col_rd-col_lu,2]=ds.GetRasterBand(3).ReadAsArray(col_lu,row_lu,(col_rd-col_lu),(row_rd-row_lu))
        panel[reg_row:reg_row+row_rd-row_lu,reg_col:reg_col+col_rd-col_lu,3]=ds.GetRasterBand(5).ReadAsArray(col_lu,row_lu,(col_rd-col_lu),(row_rd-row_lu))
        panel[reg_row:reg_row+row_rd-row_lu,reg_col:reg_col+col_rd-col_lu,4]=ds.GetRasterBand(9).ReadAsArray(col_lu,row_lu,(col_rd-col_lu),(row_rd-row_lu))

    ds=None
    return panel

def obtain_max_region_xy(region_left_up,region_right_down,orthophoto_dir,nadir_image):
    """ calculate the size of largest region [x,y] 
    
    Parameters
    ----------
    region_left_up: left_up matrix coordinates [x,y], n*36
    region_right_down:[x,y], n*36
    orthophoto_dir: image path

    Returns:
    --------
    max_region: [x,y]
    """

    pixelPair_left_up=latLonToPixel(img_addr,region_left_up)
    pixelPair_right_down=latLonToPixel(img_addr,region_right_down)
    region_size=pixelPair_right_down-pixelPair_left_up
    max_region=np.max(region_size,axis=0)+[mag*2,mag*2]    
    return max_region
    
def obtain_panel_size_xy(max_region_xy,margin_xy):
    max_x=(max_region_xy[0]+mag*2)*12+margin_xy[0]*13
    max_y=(max_region_xy[1]+mag*2)*3+margin_xy[1]*3
    return max_x,max_y

def obtain_region_loc_xy(max_region_xy,margin_xy):
    region_loc_xy=np.zeros((36,2))
    for num in range(1,37):
        y,x=divmod(num,12)
        if x==0:
            x=12
            y=y-1
        loc_x=margin_xy[0]*(x)+(max_region_xy[0]+mag*2)*(x-1)
        loc_y=margin_xy[1]*(y)+(max_region_xy[1]+mag*2)*(y)
        region_loc_xy[num,:]=[loc_x,loc_y]
        #region_loc_xy.append([loc_x,loc_y])
    return region_loc_xy

# The following method translates given latitude/longitude pairs into pixel locations on a given GEOTIF
# INPUTS: geotifAddr - The file location of the GEOTIF
#         latLonPairs - The decimal lat/lon pairings to be translated in the form [[lat1,lon1],[lat2,lon2]]
# OUTPUT: The pixel translation of the lat/lon pairings in the form [[x1,y1],[x2,y2]]
# NOTE:   This method does not take into account pixel size and assumes a high enough 
#  image resolution for pixel size to be insignificant
def latLonToPixel(geotifAddr, latLonPairs):
    """ convert point in WGS-84 coordinate to pixel matrix coordinate    
    Parameters
    ----------
    geotifAddr: geotif image path in the format  #dir="I:\\xxx\\xxx.tif"
    latLonPairs: [lat,lon] in WGS-84 coordinate, with shape 2*n
    
    Returns
    -------
    pixelPairs: [x,y] in image matrix coordinate, with shape 2*n
    """
    
    # Load the image dataset
    ds = gdal.Open(geotifAddr)
    # Get a geo-transform of the dataset
    gt = ds.GetGeoTransform()
    # Create a spatial reference object for the dataset
    srs = osr.SpatialReference()
    srs.ImportFromWkt(ds.GetProjection())
    # Set up the coordinate transformation object
    srsLatLong = srs.CloneGeogCS()
    ct = osr.CoordinateTransformation(srsLatLong,srs)
    # Go through all the point pairs and translate them to latitude/longitude pairings
    
    pixelPairs=np.zeros(latLonPairs.shape)
    point=[0,0]
    for i in range(0,latLonPairs.shape[0]):
        # Change the point locations into the GeoTransform space
        (point[1],point[0],holder) = ct.TransformPoint(latLonPairs[i,1],latLonPairs[i,0])
        # Translate the x and y coordinates into pixel values
        x = (point[1]-gt[0])/gt[1]
        y = (point[0]-gt[3])/gt[5]
        # Add the point to our return array
        #  ------------x----------
        #  |
        #  y
        #  |
        #  |
        # this is the direction cord in gdal system
        pixelPairs[i,1]=y
        pixelPairs[i,0]=x
        #pixelPairs.append([int(x),int(y)])
    ds=None
    return pixelPairs

if __name__ == "__main__":
    main()
    
