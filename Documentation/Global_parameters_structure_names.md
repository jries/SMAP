# Global parameter structure names

## getPar('name')

* V=getPar(‘name’)

If:

        uicontrol={‘edit’ ,‘checkbox’, ‘togglebutton’} 
        
Then V is numerical, logical or string.

If: 

        uicontrol={‘popupmenu’ ,‘listbox’} 
        
Then:
        
 * V.String: cell array of all entries. 
 * V.Value: selected entry (number). 
 * V.selection: String of selected entry.

## For every channel X

* '',layerX_: all paramters from Layer x GUI:
* 'channels': vector of channels to be plotted
* 'ch_filelist': file to be displayed
* 'frame_max', 'frame_min': range of frames used
* 'groupcheck': if grouping is used
* 'imax': maximum for scaling
* 'imaxtoggle': if Value=0: imax is maximum, for .Value=1: imax is quantile
* 'layercheck': display this layer
* 'lut': look up table
* 'remout': remove localizations outside color range?
* 'renderfield': if render_colormode == field: use this field for color.
* 'rendermode': Which renderer to use: Gauss/Histogram/diffraction limited/other
* 'render_colormode': If to use specific field to define color.
* 'selectedField': Field used for color
* 'shiftxy_max', 'shiftxy_min': shift of specific layer
* 'znm_max','znm_min': z-range

##

* 'layerX_filtertable': Table with min / max settings for filter

##

* 'cam_pixelsize_nm': size of camera pixel
* 'currentfileinfo': structure with info on last loaded sml file.
* 'filelist_localize': file which to localize (tiffs)
* 'filelist_long': list of sml files loaded. Including path
* 'filelist_short': list of sml files loaded. Without path.
* 'group_dt':'group_dx': parameters for grouping
* 'layerson': vector which layers are on
* 'linewidth_roi': Width of linear roi (nm)
* 'locFields': fieldnames of interfaces.LocalizationData.loc: fields of localizations

## Global parameters used by localization workflow modules

* 'loc_cameraSettings': structure with camera settings
* 'loc_currentframe': last frame and frame number opened in tiff reader
* 'loc_fileinfo': info on tiff file
* 'loc_fitOnBackground': if to fit on background
* 'loc_metadatafile': path to metadata file
* 'loc_outputfig': handle to output figure
* 'loc_preview': if this run is a preview
* 'loc_previewframe': which frame to do the preview on
* 'loc_previewmode': display-chooser: what to display
* 'loc_ROIsize': size of ROI for fitting (used by ROI cutter).

##

* 'mainGui': main GUI obj
* 'mainGuihandle': handle of main figure
* 'mainfile': last loaded file. Usually used to define default directory 
* 'numberOfLayers': now many layers are there
* 'ov_axes': handle to axis of overview image
* 'menu_plugins': structure with all plugins including their path
* 'transformationfile': file name for last use transformation


## ROI mangager paramters

* 'se_cellfov'
* 'se_cellpixelsize'
* 'se_drawboxes'
* 'se_imax'
* 'se_imaxcheck'
* 'se_rotate'
* 'se_sitefov'
* 'se_sitepixelsize'
* 'se_siteroi'
* 'se_viewer'

## Reconstruction parameters

* 'sr_axes': handle of axis of SR figure
* 'sr_figurehandle': handle of figure
* 'sr_figurenumber': number of figure
* 'sr_image': image structure
* 'sr_imagehandle': handle to image
* 'sr_imagesize': manually set size of reconstructed image in pixels
* 'sr_imsizecheck': if to use manually set size
* 'sr_layerson': vector which layers are on
* 'sr_pixfactor': if binning is used: factor
* 'sr_pixrec': pixelsize in nm
* 'sr_pos': position of sr image (nm).
* 'sr_roihandle': handle to roi
* 'sr_size':  size of sr image in nm. (divided by 2!): Roi=pos-size:pos+size
* 'sr_sizeRecPix': size of sr image in pixels
* 'status': string to be displayed as status

 
## Results

* 'counting_histogram'

