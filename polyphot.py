#!/usr/bin/python
#=============================================================================#
# NAME:     polyphot.py                                                       #
#                                                                             #
# PURPOSE:  Script to to load a FITS image and allow the user to measure the  #
#           flux of a source by drawing a polygon around the emission.        #
#                                                                             #
# UPDATED:  08-Nov-2010                                                       #
#                                                                             #
# TODO:     Finished?                                                         #
#=============================================================================#
import sys, copy, optparse, os, sys, re, string, math, getopt
from matplotlib.collections import RegularPolyCollection
from matplotlib.collections import PatchCollection
from matplotlib.patches import Polygon
import matplotlib.cm as cm
from matplotlib.widgets import Button, Slider
from matplotlib.pyplot import figure, show
from pylab import *
import pyfits
import pywcs
import numpy as np
import numpy.ma as ma


# Plot the annotation file if it exists?
doplotann=True

#-----------------------------------------------------------------------------#
# Main control function                                                       #
#-----------------------------------------------------------------------------#
def main():
    
    # Get the filename from the command line
    args = sys.argv[1:]
    if len(args)<1: usage()
    else: filename = args[0]
    
#---------------- Load the FITS file and determine coordinates ---------------#

    [header,xydata]=load_fits_image(filename)
    wcs = pywcs.WCS(header)
    
    # Read the max and min of the data array
    zmin = (np.nan_to_num(xydata)).min()
    zmax = (np.nan_to_num(xydata)).max()

#---------------- Load the annotation file if it exists ----------------------#

    suffix = '.fits'
    if filename.endswith(suffix):
        annfilename = filename[:-len(suffix)] + '.ann'
        if os.path.exists(annfilename) and doplotann:
            ann_list = kvis_ann_parse(annfilename)
            ann_list = wcs.wcs_sky2pix(ann_list, 1)

#------------------------ Setup the figure for plotting ----------------------#

    # Define a figure on which to plot the image
    fig = figure()
    axplot = fig.add_subplot(1,1,1)
    axplot.set_title("Click to outline the source emission.")
    subplots_adjust(bottom=0.12)

    # Format the axis labels
    #majorFormatterX = FuncFormatter(label_format_deg) 
    #majorFormatterY = FuncFormatter(label_format_deg)
    #axplot.xaxis.set_major_formatter(majorFormatterX)
    #axplot.yaxis.set_major_formatter(majorFormatterY)
    
    # Plot the image data and colorbar
    cax = axplot.imshow(xydata,interpolation='nearest',origin='lower',
                        cmap=cm.jet,vmin=zmin,vmax=zmax)
    cbar=colorbar(cax,pad=0.0)

    # Plot the fits annotation
    if os.path.exists(annfilename) and doplotann:
        
        # Create a circle for each Gaussian fit
        ann = RegularPolyCollection(
            10,
            rotation=0,
            sizes=(5,),
            facecolors = 'none',
            edgecolors = 'black',
            linewidths = (1,),
            offsets = ann_list,
            transOffset = axplot.transData
            )
        axplot.add_collection(ann)
        
    # Add the buttons to the bottom edge
    axreset = axes([0.49, 0.025, 0.1, 0.04])
    axmeasure = axes([0.60, 0.025, 0.1, 0.04])
    axsave = axes([0.71, 0.025, 0.1, 0.04])
    breset = Button(axreset, 'Reset')
    bmeasure = Button(axmeasure, 'Measure')
    bsave = Button(axsave, 'Save')

#-------------------------- Polygon editor class -----------------------------#

    class ThreePolyEditor:
        """
        Left-click to add a point to the current polygon. Middle-click to delete
        the last point and right-click to close the polygon. Three polygons
        addedin sequence: source-aperture, sky aperture and exclusion aperture.
        """
    
        def __init__(self):
            self.mode = 'src' 
        
            # Lists to store the vertices of the source, sky background and
            # exclusion-zone polygons
            self.offsets_src = []
            self.offsets_sky = []
            self.offsets_exc = []        
            self.offsets = []        # general working array
            self.apertures = []      # completed polygons for plotting
            self.measurements = {}

            # Working polygon collection (points plotted as small polys)
            self.points = RegularPolyCollection(
                10,
                rotation=0,
                sizes=(50,),
                facecolors = 'white',
                edgecolors = 'black',
                linewidths = (1,),
                offsets = self.offsets,
                transOffset = axplot.transData
                )

        # Handle mouse clicks
        def onclick(self,event):

            # Disable click events if using the zoom & pan modes.
            if fig.canvas.widgetlock.locked():
                return
            

            if event.button==1:
            
                # Left-click on plot = add a point
                if axplot.contains(event)[0] and self.mode != 'done':
                    x = event.xdata
                    y = event.ydata
                    self.offsets.append((x,y))
                    self.points.set_offsets(self.offsets)
                    if len(self.offsets)==1 : axplot.add_collection(self.points)
                    self.update()
                
                # Right-click on wedge = halve the lower colour clip
                if cbar.ax.contains(event)[0]:
                    clims = cax.get_clim()
                    cax.set_clim(clims[0]/2.0,clims[1])
                    self.update()
                
            elif event.button==2:
            
                # Middle-click = delete last created point
                if axplot.contains(event)[0] and len(self.offsets)>0 \
                       and self.mode != 'done':
                    self.offsets.pop(-1)
                    if len(self.offsets)==0:
                        self.points.remove()
                    else:
                        self.points.set_offsets(self.offsets)
                    self.update()
                
                # Middle-click on wedge = reset the colour clip
                if cbar.ax.contains(event)[0]:
                    clims = cax.get_clim()
                    cax.set_clim(zmin,zmax)
                    self.update()
                
            if event.button==3:
            
                # Right-click on plot = complete the polygon
                if  axplot.contains(event)[0] and len(self.offsets)>2 \
                       and self.mode != 'done':
                    cpoly = Polygon(self.offsets, animated=False,linewidth=2.0)
                    if self.mode == 'src':
                        cpoly.set_edgecolor('white')
                        cpoly.set_facecolor('none')
                    elif self.mode == 'sky':
                        cpoly.set_edgecolor('yellow')
                        cpoly.set_facecolor('none')
                    elif self.mode == 'exc':
                        cpoly.set_edgecolor('black')
                        cpoly.set_facecolor('none')
                    self.apertures.append(axplot.add_patch(cpoly))
                    self.update('complete')
            
                # Left-click on wedge = halve the upper colour clip
                if cbar.ax.contains(event)[0]:
                    clims = cax.get_clim()
                    cax.set_clim(clims[0],clims[1]/2.0)
                    self.update()
                
        # Store completed polygons in the relevant lists
        def update(self,action=''):

            # When a polygon is complete switch to the next mode
            if action == 'complete':
                if self.mode == 'src':
                    self.offsets_src = self.offsets
                    self.mode = 'sky'
                elif self.mode == 'sky':
                    self.offsets_sky = self.offsets
                    self.mode = 'exc'
                elif self.mode == 'exc':
                    self.offsets_exc = self.offsets
                    self.mode = 'done'
                self.offsets=[]
                
            # When the reset button is complete clear the polygon lists & plot
            elif action == 'reset':
                self.offsets_src = []
                self.offsets_sky = []
                self.offsets_exc = []        
                self.offsets = []
                self.mode = 'src'
                for i in range(len(axplot.collections)):
                    axplot.collections.pop()
                for aperture in self.apertures: aperture.remove()
                self.apertures = []
                self.cent_pt = None
                cax.set_clim(zmin,zmax)
                self.measurements = {}

            # Update the title
            if self.mode == 'src':
                title_str = "Click to outline the source emission."
            elif self.mode == 'sky':
                title_str = 'Click to outline a region of background sky.'
            elif self.mode == 'exc':
                title_str = 'Click to define a source exclusion zone.'
            else:
                title_str = 'Click Reset or Measure when done'
            axplot.set_title(title_str)
            fig.canvas.draw()

        # Reset the plot and clear polygons
        def doreset(self, event):
            self.update('reset')

        # Measure the polygon definitions
        def domeasure(self, event):

            if not self.offsets_src or not self.offsets_sky:
                print "Please complete both the source and sky polygons."
                return
            
            # Measure the source properties
            # [sum,avg,stdev,max,min,npts,centroid_px]
            # [0   1   2     3   4   5    6          ]
            print "\nGetting pixels in source polygon ...",
            sys.stdout.flush() 
            src_indices = get_pix_in_poly(self.offsets_src,xydata)
            print "done." ; sys.stdout.flush() 
            print "Measuring properties of the source ...",
            sys.stdout.flush() 
            src_prop_raw = measure_pix_props(src_indices,xydata,docent=True)
            print "done."; sys.stdout.flush() 
            
            # Calculate the source polygon area and equivalent radius
            pix_scale_x = abs(header['CD1_1'])
            pix_scale_y = abs(header['CD2_2'])
            pix_scale_deg = (pix_scale_x+pix_scale_y)/2.0
            area_px = poly_area(self.offsets_src)
            area_degsq = area_px*(pix_scale_deg**2)
            diameter_px = 2.0*math.sqrt(area_px/math.pi)
            diameter_deg = diameter_px*pix_scale_deg
            
            # Measure the sky properties
            print "Getting pixels in sky polygon ...",
            sys.stdout.flush() 
            sky_indices = get_pix_in_poly(self.offsets_sky,xydata)
            print "done."; sys.stdout.flush() 
            print "Measuring properties of the sky ...",
            sys.stdout.flush() 
            sky_prop_raw  = measure_pix_props(sky_indices,xydata)
            print "done."; sys.stdout.flush() 

            # Calculate the beam area in pixels
            pix_per_beam = calc_beam_area(header)

            # Calculate the integrated flux
            integ_flux = (src_prop_raw[0]-(sky_prop_raw[1]*src_prop_raw[5]))\
                         /pix_per_beam

            # Calculate the uncertainty in the integrated flux
            # F.Masci, IPAC: 'Flux-Uncertainty from Aperture Photometry'
            dinteg_flux = math.sqrt(src_prop_raw[5]*(sky_prop_raw[2]**2)*
                                    (1.0+src_prop_raw[5]/sky_prop_raw[5]))
            
            # Calculate the centroid in world coordinates
            [wcentroid_deg] = wcs.wcs_pix2sky([src_prop_raw[6]], 1)
            [centroid_deg] = wcs.wcs_pix2sky([src_prop_raw[7]], 1)
            [centmax_deg] = wcs.wcs_pix2sky([src_prop_raw[8]], 1)

            # Convert the polygons to world coordinates            
            offsets_src_wld = wcs.wcs_pix2sky(self.offsets_src,1)
            offsets_sky_wld = wcs.wcs_pix2sky(self.offsets_sky,1) 
            try:
                offsets_exc_wld = wcs.wcs_pix2sky(self.offsets_exc,1)
            except Exception:
                offsets_exc_wld = np.array([])

            # Calculate final numbers for the user
            self.measurements['FLUX_JY'] = integ_flux
            self.measurements['dFLUX_JY'] = dinteg_flux
            self.measurements['BGAVG_JYBM'] = sky_prop_raw[1]
            self.measurements['BGSTD_JYBM'] = sky_prop_raw[2]
            self.measurements['MAX_JYBM'] = src_prop_raw[3]
            self.measurements['MIN_JYBM'] = src_prop_raw[4]
            self.measurements['NPIXSRC'] = src_prop_raw[5]
            self.measurements['NPIXSKY'] = sky_prop_raw[5]
            self.measurements['WCENT_DEG'] = wcentroid_deg
            self.measurements['CENT_DEG'] = centroid_deg
            self.measurements['CENTMAX_DEG'] = centmax_deg
            self.measurements['DIAM_DEG'] = diameter_deg
            self.measurements['AREA_DEGSQ'] = area_degsq
            self.measurements['BEAMAREA_PX'] = pix_per_beam
            self.measurements['POLYSRC'] = offsets_src_wld.flatten().tolist()
            self.measurements['POLYSKY'] = offsets_sky_wld.flatten().tolist()
            self.measurements['POLYEXC'] = offsets_exc_wld.flatten().tolist()

            # Feeback to user
            print ""
            for key,val in self.measurements.iteritems():
                if not (key=='POLYSRC' or  key=='POLYSKY' or key=='POLYEXC'):
                    print key, "        \t=\t", val
            print ""

            # Plot the centroid of the pixels
            self.pcentre = RegularPolyCollection(
                10,
                rotation=0,
                sizes=(50,),
                facecolors = 'white',
                edgecolors = 'black',
                linewidths = (1,),
                offsets = [src_prop_raw[6]],
                transOffset = axplot.transData
                )
            axplot.add_collection(self.pcentre)
            self.update()
            
        # Save the measurements to a file
        def dosavefile(self, event):

            # Strip the '.fits' to form the datfile name
            suffix = '.fits'
            if filename.endswith(suffix):
                datfilename = filename[:-len(suffix)] + '.polyflux.dat'
            else:
                datfilename = filename + '.polyflux.dat'

            if not self.measurements:
                print "Measurements empty! Click measure before saving."
                return

            # Open the datfile for appending
            if os.path.exists(datfilename):
                print "\nAppending measurement to file '%s' ... "%datfilename,
            else:
                print "\nSaving measurement to NEW file '%s' ... "%datfilename,
            sys.stdout.flush()
            datHandle = open(datfilename, 'a')
            datHandle.write('START\n')
            for key,val in self.measurements.iteritems():
                datHandle.write("%s=%s\n" % (key,val))
            datHandle.write('END\n\n')
            datHandle.close()
            print "done."
            
#-----------------------------------------------------------------------------#

    # Make an instance of the poly-editor and bind events to methods
    editor = ThreePolyEditor()
    fig.canvas.mpl_connect('button_press_event', editor.onclick)
    breset.on_clicked(editor.doreset)
    bmeasure.on_clicked(editor.domeasure)
    bsave.on_clicked(editor.dosavefile)

    # Draw the plot to the screen
    show()


#=============================================================================#


#-----------------------------------------------------------------------------#
# Calculate the sum, averag and stdev of selected pixels                      #
#-----------------------------------------------------------------------------#
def measure_pix_props(indices,array, docent=False,wcent=True):
    mask = ones_like(array)
    for x,y in indices: mask[y,x]=0.0
    marray = ma.masked_array(array, mask)
    summ = marray.sum()
    avg = marray.mean()
    stdev = marray.std()
    maximum = marray.max()
    minimum = marray.min()
    centmax = list(np.unravel_index(marray.argmax(), marray.shape))
    centmax.reverse()

    # Calculated the intensity weighted centroid
    cent = 0
    wcent = 0
    if docent:
        # Calculate the weighted centroid
        x_sum=0.0; y_sum=0.0; weight_sum=0.0
        for i,j in indices:
            x_sum+=(i*array[j,i])
            y_sum+=(j*array[j,i])
            weight_sum+=array[j,i]
        wcent = ((x_sum/weight_sum),(y_sum/weight_sum))
        
        # Calculate the normal centroid
        xind = np.array(indices,dtype='float')[:,0]
        yind = np.array(indices,dtype='float')[:,1]
        cent = (sum(xind)/len(xind),sum(yind)/len(yind))
    
    return [summ,avg,stdev,maximum,minimum,len(indices),wcent,cent,centmax]
    


#-----------------------------------------------------------------------------#
# Return the indices of the pixels inside a polygon                           #
#-----------------------------------------------------------------------------#
def get_pix_in_poly(poly,array):

    indices = []

    # Get the bounds of the poly vertices
    px = np.array(poly)[:,0]
    py =  np.array(poly)[:,1]
    px_min = int(floor(min(px)))
    px_max = int(ceil(max(px)))
    py_min = int(floor(min(py)))
    py_max = int(ceil(max(py)))
    
    for i in range(px_min,px_max+1):
        for j in range(py_min,py_max+1):
            if point_in_poly(i,j,poly):
                indices.append((i,j))

    return indices
            

#-----------------------------------------------------------------------------#
# Query if a point is inside the polygon                                      #
#-----------------------------------------------------------------------------#
def point_in_poly(px,py,poly):

        cn = 0

        i = -1
        j = 0
        while j < len(poly):
            qx, qy = poly[i]
            rx, ry = poly[j]

            if (px, py) == (qx, qy):
                return True

            if (    ((qy <= py) and (ry > py)) or \
                    ((qy > py) and (ry <= py))    ):
                vt = (py - qy) / (ry - qy)
                if (px < qx + vt * (rx - qx)):
                    cn += 1

            i = j
            j += 1

        return cn%2 


#-----------------------------------------------------------------------------#
# Load a fits image and return the data and dimensions                        #
#-----------------------------------------------------------------------------#
def load_fits_image(filename):
    
    # Read the header and image data from the file
    header=pyfits.getheader(filename)
    data = pyfits.getdata(filename)
    naxis = len(data.shape)

    # Strip unused dimensions from the data array
    if naxis == 2:
        xydata = data.copy()
        del data
    elif naxis == 3:
        xydata = data[0].copy()
        del data
    elif naxis == 4:
        xydata = data[0][0].copy()
        del data
    elif naxis == 5:
        xydata = data[0][0][0].copy()
        del data
    else:
        print "Data array contains %s axes" % naxis 
        print "This script supports up to 5 axes only."
        sys.exit(1)

    # Strip unused dimensions from the header
    stripHeadKeys = ['NAXIS','CRVAL','CRPIX','CDELT','CTYPE','CROTA',
                     'CD1_','CD2_']
    for i in range(3,naxis+1):
        for key in stripHeadKeys:
            del header[key+str(i)]
    header.update('NAXIS',2)

    # Determine the coordinate type and the corresponding keyword index
    # Make a note in the header
    ra_regx = re.compile('^RA')
    dec_regx = re.compile('^DEC')
    glon_regx = re.compile('^GLON')
    glat_regx = re.compile('^GLAT')
    for i in range(int(header['NAXIS'])):
        keyword = "CTYPE"+str(i+1)
        if ra_regx.match(header[keyword]): coord_type="EQU"; x_indx=i+1
        if dec_regx.match(header[keyword]): y_indx=i+1
        if glon_regx.match(header[keyword]): coord_type="GAL"; x_indx=i+1
        if glat_regx.match(header[keyword]): y_indx=i+1
    if not x_indx or not y_indx:
        print "Failed to find Equatorial or Galactic axis coordinate types."
        del data; del header
        sys.exit(1)
    else:
        header.update('XINDEX', x_indx)
        header.update('YINDEX', y_indx)

    # Convert AIPS clean-beam types to standard BMAJ, BMIN, BPA
    try:
        bmaj = header['CLEANBMJ']
        bmin = header['CLEANBMN']
        bpa = header['CLEANBPA']
        header.update('BMAJ',bmaj)
        header.update('BMIN',bmin)
        header.update('BPA',bpa)
    except Exception:
        bmaj,bmin,bpa = get_beam_from_history(header)
        if bmaj:
            header.update('BMAJ',bmaj)
            header.update('BMIN',bmin)
            header.update('BPA',bpa)
        else:
            print "No AIPS style beam keywords found."
        
    # Check for PIXSCAL keywords and write to CDELT standard
    try:
        xdelt=(-1)*(header['PIXSCAL'+str(x_indx)])/3600.0
        ydelt=(header['PIXSCAL'+str(y_indx)])/3600.0
        header.update('CD1_'+str(x_indx), xdelt)
        header.update('CD2_'+str(y_indx), ydelt)
    except Exception:
        pass

    return [header,xydata]


#-----------------------------------------------------------------------------#
# Calculate the effective beam area in px given the Gaussian FWHMs            #
#-----------------------------------------------------------------------------#
def calc_beam_area(header):

    # Read the beam and pixel dimensions
    try:
        bmaj = header['BMAJ']
        bmin = header['BMIN']
        bpa = header['BPA']
    except Exception:
        print "WARINIG: Using hardcoded beam size (1.5'')."
        bmaj = 1.5/3600.0
        bmin = 1.5/3600.0
        bpa = 0.0
        
    pix_scale_x = header['CD1_1']
    pix_scale_y = header['CD2_2']

    # Convert FWHMs to sigma
    gfactor = 2.0*math.sqrt(2*math.log(2))  # FWHM=gfactor*sigma
    sigma_maj = bmaj/gfactor
    sigma_min = bmin/gfactor
    
    # Calculate the beam area in the input units and in pixels
    beam_area_sq = 2.0*math.pi*sigma_maj*sigma_min
    pix_area_sq = abs(pix_scale_x*pix_scale_y)
    beam_area_pix = beam_area_sq/pix_area_sq
    
    return beam_area_pix    



#-----------------------------------------------------------------------------#
# Search the FITS HISTORY fields for the AIPS clean beam string               #
#-----------------------------------------------------------------------------#
def get_beam_from_history(header):

    bmaj=None; bmin=None; bpa=None
    history=header.get_history()
    history.reverse()

    #'AIPS   CLEAN BMAJ=  4.3403E-04 BMIN=  3.1039E-04 BPA= -11.55'
    beamHistStr = 'AIPS\s+CLEAN\sBMAJ=\s+(\S+)\s+BMIN=\s+(\S+)\s+BPA=\s+(\S+)'
    bmHistPat = re.compile(beamHistStr)
    
    for i in range(len(history)):
        
        m=bmHistPat.match(history[i])
        if m:
            bmaj = float(m.group(1))
            bmin = float(m.group(2))
            bpa  = float(m.group(3))
            break

    return [bmaj,bmin,bpa]

            
#-----------------------------------------------------------------------------#
# Calculate the area of a polygon                                             #
#-----------------------------------------------------------------------------#
def poly_area(p):
    if not p[0] == p[-1]:
        p.append(p[0])
    return 0.5 * abs(sum(x0*y1 - x1*y0
                         for ((x0, y0), (x1, y1)) in segments(p)))
def segments(p):
    return zip(p, p[1:] + [p[0]])


#-----------------------------------------------------------------------------#
# Parse the ellipses from a Kis .ann file                                     #
#-----------------------------------------------------------------------------#
def kvis_ann_parse(filename):
    coord_list=[];
    ell = re.compile('^ELLIPSE\s+W\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)')
    
    ANNFILE = open(filename, "r")
    for line in ANNFILE:
        entry={}
        line = line.rstrip("\n\r")
        m = ell.match(line)
        if m:
            entry = (float(m.group(1)),float(m.group(2)))
            coord_list.append(entry)
    return coord_list
    


#-----------------------------------------------------------------------------#
# Format decimal tick labels                                                  #
#-----------------------------------------------------------------------------#
def label_format_deg(deg, pos):
    return "%.0f" % deg


#-----------------------------------------------------------------------------#
# Usage text printed on error                                                 #
#-----------------------------------------------------------------------------#
def usage():
    print "\n\tUSAGE: polyphot.py <file.fits>\n"
    sys.exit(1)


#-----------------------------------------------------------------------------#
main()
#-----------------------------------------------------------------------------#
