import astropy.io.fits as fits
from astropy import wcs as WCS
import numpy as np
from numpy import linalg
from trippy import scamp
from os import path
import os, sys
import pylab as pyl
from astropy.visualization import interval
import numdisplay
from matplotlib.patches import Circle


#sextractor shape cut  -- done
#window size command line option -- done
#zoom on sources -- done
#Bayes information criterion -- is now printed
#parameters into the headers
#instructions


def trimCatalog(cat, minBA=0.85):

    good=[]
    for i in range(len(cat['XWIN_IMAGE'])):
        try:
            ba = cat['BWIN_IMAGE'][i]/cat['AWIN_IMAGE'][i]
            m = cat['MAG_AUTO'][i]
            flag = cat['FLAGS'][i]
        except:
            pass
        #if cat['FLAGS'][i]==0 and m>0 and m<26:
        if m>0 and m<26 and flag == 0 and ba>minBA:
            good.append(i)
    good = np.array(good)

    #(X,Y,A,B) = (cat['XWIN_IMAGE'][good],cat['YWIN_IMAGE'][good],cat['AWIN_IMAGE'][good],cat['BWIN_IMAGE'][good])
    #w = np.where(B/A>minBA)
    #good = good[w]

    outcat = {}
    for i in cat:
        outcat[i] = cat[i][good]
    return outcat


d2r= np.pi/180.0

class matrixWCSSolver(object):

    def __init__(self, imagedata, header, imageSources, refSources, xo = 512.0, yo = 512.0, windowSize = 13):

        if 'B_DEC_7' in header:
            self.b_ra = []
            self.b_dec = []
            for i in range(8):
                self.b_ra.append(header['B_RA_{}'.format(i)])
                self.b_dec.append(header['B_DEC_{}'.format(i)])
            self.b_ra = np.array(self.b_ra)
            self.b_dec = np.array(self.b_dec)

        self.imageSources = imageSources
        self.refSources = refSources
        self._xo = xo
        self._yo = yo
        self.header = header
        self.imagedata = np.copy(imagedata)

        #scratch crap pixel value fix
        mean = np.median(self.imagedata[::3,::3])
        w = np.where(self.imagedata<0)
        self.imagedata[w] = mean

        (self._z1,self._z2) = numdisplay.zscale.zscale(self.imagedata, contrast = 0.4)
        self._normer = interval.ManualInterval(self._z1,self._z2)

        self._wcs  = WCS.WCS(header)
        self._refSourcePix = self._wcs.wcs_world2pix(self.refSources[:,:2],0)

        (a,b) = self.imagedata.shape
        w = np.where((self._refSourcePix[:,0]>-100)&(self._refSourcePix[:,0]<b+100)&(self._refSourcePix[:,1]>-100)&(self._refSourcePix[:,1]<a+100))
        self._refSourcePix = self._refSourcePix[w]
        self.refSources = self.refSources[w]

        self._initSourceSelection = None
        self._initRefSelection = None
        self.matches = []

        self.windowSize = windowSize
        self._lastKilled = []

    def initialMatch(self, maxDeltaMag = 1.5):
        fig = pyl.figure('Full Image', figsize = (self.windowSize, self.windowSize))
        sp = fig.add_subplot(111)
        implot = pyl.imshow(self._normer(self.imagedata))
        implot.set_cmap('hot')
        pyl.scatter(self._refSourcePix[:,0], self._refSourcePix[:,1], c = 'g', s = 15)
        #pyl.scatter(self.imageSources[:,0]-1, self.imageSources[:,1]-1,c='b',alpha = 0.5)


        self.circles = []
        for i in range(len(self.imageSources)):
            circle=Circle(self.imageSources[i,:2]-np.ones(2),20,facecolor="none",edgecolor='b',linestyle='dashed',linewidth=2, alpha=0.75,zorder=10)
            sp.add_patch(circle)
            self.circles.append(circle)
        sp.invert_yaxis()
        sp.set_xlim(0,self.imagedata.shape[1])
        sp.set_ylim(0,self.imagedata.shape[1])

        pyl.connect('button_press_event',self._getStar)

        pyl.show()

    def _getStar(self, event, maxDist =  10, maxMagDiff = 1.0):
        if event.button==1:
            rcx = event.xdata
            rcy = event.ydata
            if rcx==None or rcy==None: return

            distToImageSource = ((self.imageSources[:,0] - (rcx + 1))**2 + (self.imageSources[:,1] - (rcy + 1))**2)**0.5
            argMin = np.argmin(distToImageSource)
            minDist = distToImageSource[argMin]

            self._initSourceSelection = self.imageSources[argMin]
            (x,y,m) = self._initSourceSelection[:3]
            print '   Selected image source at {:.3f}, {:.3f} with magnitude {:.2f}.'.format(x,y,m)

        elif event.button==3:
            rcx = event.xdata
            rcy = event.ydata
            if rcx==None or rcy==None: return

            distToImageSource = ((self._refSourcePix[:,0] - (rcx + 1))**2 + (self._refSourcePix[:,1] - (rcy + 1))**2)**0.5
            argMin = np.argmin(distToImageSource)
            minDist = distToImageSource[argMin]

            self._initRefSelection = argMin
            (x,y) = self._refSourcePix[argMin]
            (ra,dec,r,g) = self.refSources[argMin,:4]
            print '   Selected reference star at {:.3f}, {:.3f} with magnitudes r={:.2f} g={:.2f}.'.format(x,y,r,g)

        if self._initRefSelection is not None and self._initSourceSelection is not None:

            self.matches = []

            off_x = self._initSourceSelection[0] - self._refSourcePix[self._initRefSelection,0]
            off_y = self._initSourceSelection[1] - self._refSourcePix[self._initRefSelection,1]

            deltaMags = []
            for i in range(len(self.refSources)):
                (ox,oy) = self._refSourcePix[i,:2] - np.ones(2)
                (nx,ny) = (ox + off_x, oy + off_y)

                dist = ((nx - self.imageSources[:,0])**2 + (ny - self.imageSources[:,1])**2)**0.5
                argMin = np.argmin(dist)
                if dist[argMin] < maxDist:
                    deltaMags.append(self.refSources[i,2]-self.imageSources[argMin,2])

            deltaMags = np.array(deltaMags)
            correctedMags = self.imageSources[:,2] + np.median(deltaMags)


            for i in range(len(self.refSources)):
                (ox,oy) = self._refSourcePix[i,:2] - np.ones(2)
                (nx,ny) = (ox + off_x, oy + off_y)

                dist = ((nx - self.imageSources[:,0])**2 + (ny - self.imageSources[:,1])**2)**0.5
                argMin = np.argmin(dist)
                if dist[argMin] < maxDist and (correctedMags[argMin] - self.refSources[i,2])<maxMagDiff:

                    pyl.plot([ox,nx], [oy,ny], 'b-', lw=3)

                    self.matches.append([self.refSources[i,0], self.refSources[i,1],self.imageSources[argMin,0], self.imageSources[argMin,1],self.imageSources[argMin,4],self.imageSources[argMin,5]])
            self.matches = np.array(self.matches)
            pyl.draw()

    def _resid4Panel(self):
        fig = pyl.figure(' Fourpanel', figsize = (self.windowSize, self.windowSize))
        fig.subplots_adjust(hspace = 0, wspace = 0)
        self._sp1 = fig.add_subplot(221,xticklabels = '')
        self._sp2 = fig.add_subplot(222,xticklabels = '',yticklabels = '')
        self._sp3 = fig.add_subplot(223)
        self._sp4 = fig.add_subplot(224,yticklabels = '')

        self.dra = (self.RA - self.pra)*3600.0
        self.ddec = (self.DEC - self.pdec)*3600.0

        self._sp1.scatter(self.X,self.dra)
        self._sp3.scatter(self.X,self.ddec)
        self._sp2.scatter(self.Y,self.dra)
        self._sp4.scatter(self.Y,self.ddec)

        self._sp4.set_xlabel('Y')
        self._sp3.set_xlabel('X')
        self._sp3.set_ylabel('delta Dec (")')
        self._sp1.set_ylabel('delta RA (")')

        w = np.where(self.goodMatches)
        self._sp1.set_title('{:.3f}"'.format(np.std(self.dra[w])))
        self._sp2.set_title('{:.3f}"'.format(np.std(self.ddec[w])))

        pyl.connect('button_press_event', self._killResid)
        pyl.connect('key_press_event', self._zoomResid)

        self._fullZooms = [self._sp1.get_xlim(),self._sp1.get_ylim(), self._sp2.get_xlim(),self._sp2.get_ylim(), self._sp3.get_xlim(),self._sp3.get_ylim(), self._sp4.get_xlim(),self._sp4.get_ylim()]
        self._zoomed = False
        pyl.show()

    def _zoomResid(self, event):

        if event.key in ['z','Z'] or event.key in ['r','R']:
            if not self._zoomed and event.key in ['z','Z']:
                w = np.where(self.goodMatches)

                x_min = np.min(self.X[w])
                x_max = np.max(self.X[w])
                y_min = np.min(self.Y[w])
                y_max = np.max(self.Y[w])
                r_min = np.min(self.dra[w])-0.005
                r_max = np.max(self.dra[w])+0.005
                d_min = np.min(self.ddec[w])-0.005
                d_max = np.max(self.ddec[w])+0.005

                self._zoomed = True
            else:
                (r_min,r_max) = self._fullZooms[1]
                (d_min,d_max) = self._fullZooms[3]

                self._zoomed = False

            self._sp1.set_ylim(r_min, r_max)
            self._sp2.set_ylim(r_min, r_max)
            self._sp3.set_ylim(d_min, d_max)
            self._sp4.set_ylim(d_min, d_max)

            pyl.draw()

        if event.key in ['r', 'R']:
            self.goodMatches *= 0
            self.goodMatches += 1
            self._lastKilled = []
            self._fullRedraw()

        elif event.key in ['k','K']:
            w = np.where(self.goodMatches)
            (A, b_ra, b_dec, predRA, predDEC, lnL_RA, lnL_DEC, BICS) = self._solveMatrix()
            std_ra = np.std(predRA[w] - self.RA[w])*3600.0
            std_dec = np.std(predDEC[w] - self.DEC[w])*3600.0
            delta = (std_ra**2 + std_dec**2)**0.5

            k_deltas = []
            for i in w[0]:
                self.goodMatches[i] = 0

                fake_w = np.where(self.goodMatches)
                (A, b_ra, b_dec, predRA, predDEC, lnL_RA, lnL_DEC, BICS) = self._solveMatrix()
                k_std_ra = np.std(predRA[fake_w] - self.RA[fake_w])*3600.0
                k_std_dec = np.std(predDEC[fake_w] - self.DEC[fake_w])*3600.0
                k_delta = (k_std_ra**2 + k_std_dec**2)**0.5
                k_deltas.append(k_delta)

                self.goodMatches[i] = 1
            k_deltas = np.array(k_deltas)
            argmin = np.argmin(k_deltas)
            self._lastKilled.append(w[0][argmin])
            self.goodMatches[w[0][argmin]] = 0
            self._fullRedraw()

        elif event.key in ['j','J'] and len(self._lastKilled) > 0:
            if self.goodMatches[self._lastKilled[-1]] == 0:
                self.goodMatches[self._lastKilled[-1]] = 1
                self._fullRedraw()
            self._lastKilled = self._lastKilled[:-1]

        elif event.key == '?':
            print 'HELP!!!'
            print 'Closing the window will accept the current fit as is.'
            print
            print " -press z to zoom in on the selected good stars or out to all"
            print " -press r to reset all stars back to good and try again"
            print " -press k to kill the most discrepant point"
            print " -press j to restore the most recent point you nuked with k"
            print " -press ? to see this message"

    def _fullRedraw(self):
        (A, b_ra, b_dec, predRA, predDEC, lnL_RA, lnL_DEC, BICS) = self._solveMatrix()


        self.pra = np.copy(predRA)
        self.pdec = np.copy(predDEC)
        self.dra = (self.RA - self.pra)*3600.0
        self.ddec = (self.DEC - self.pdec)*3600.0
        self.b_ra = b_ra[:]
        self.b_dec = b_dec[:]

        self.lnL_RA = lnL_RA
        self.lnL_DEC = lnL_DEC
        self.BIC_RA = BICS[0]
        self.BIC_DEC = BICS[1]
        self.BIC_RA_low = BICS[2]
        self.BIC_DEC_low = BICS[3]

        self._sp1.cla()
        self._sp2.cla()
        self._sp3.cla()
        self._sp4.cla()

        colours = []
        for i in range(len(self.goodMatches)):
            if self.goodMatches[i]:
                colours.append('b')
            else:
                colours.append('r')
        self._sp1.scatter(self.X,self.dra,c=colours)
        self._sp2.scatter(self.Y,self.dra,c=colours)
        self._sp3.scatter(self.X,self.ddec,c=colours)
        self._sp4.scatter(self.Y,self.ddec,c=colours)

        self._sp4.set_xlabel('Y')
        self._sp3.set_xlabel('X')
        self._sp3.set_ylabel('delta Dec (")')
        self._sp1.set_ylabel('delta RA (")')

        w = np.where(self.goodMatches)
        self.std_ra = np.std(self.dra[w])
        self.std_dec = np.std(self.ddec[w])
        self._sp1.set_title('RA residual {:.3f}"'.format(np.std(self.dra[w])))
        self._sp2.set_title('Dec residual {:.3f}"'.format(np.std(self.ddec[w])))

        self._fullZooms = [self._sp1.get_xlim(),self._sp1.get_ylim(), self._sp2.get_xlim(),self._sp2.get_ylim(), self._sp3.get_xlim(),self._sp3.get_ylim(), self._sp4.get_xlim(),self._sp4.get_ylim()]
        self._zoomed = False

        pyl.draw()



    def _killResid(self,event):
        rcx = event.xdata
        rcy = event.ydata
        if rcx==None or rcy==None: return

        if event.inaxes == self._sp1:
            dist = ((rcx - self.X)**2 + (rcy - self.dra)**2)**0.5
        elif event.inaxes == self._sp2:
            dist = ((rcx - self.Y)**2 + (rcy - self.dra)**2)**0.5
        elif event.inaxes == self._sp3:
            dist = ((rcx - self.X)**2 + (rcy - self.ddec)**2)**0.5
        elif event.inaxes == self._sp4:
            dist = ((rcx - self.Y)**2 + (rcy - self.ddec)**2)**0.5

        arg = np.argmin(dist)

        if event.button == 1:
            self.goodMatches[arg] = 0
        elif event.button == 3:
            self.goodMatches[arg] = 1

            if arg in self._lastKilled:
                del self._lastKilled[self._lastKilled.index(arg)]

        self._fullRedraw()



    def solveMatrix(self):
        if self.matches == []:
            print 'Must run _getStar first!'
            return
        self.X = self.matches[:,2] - self._xo
        self.Y = self.matches[:,3] - self._yo
        self.RA = self.matches[:,0]
        self.DEC = self.matches[:,1]
        self.dX2 = self.matches[:,4]
        self.dY2 = self.matches[:,5]




        self.goodMatches = np.ones(len(self.matches))

        (A, b_ra, b_dec, predRA, predDEC, lnL_RA, lnL_DEC, BICS) = self._solveMatrix()

        self.pra = np.copy(predRA)
        self.pdec = np.copy(predDEC)
        self.dra = (self.RA - self.pra)*3600.0
        self.ddec = (self.DEC - self.pdec)*3600.0
        self.lnL_RA = lnL_RA
        self.lnL_DEC = lnL_DEC
        self.BIC_RA = BICS[0]
        self.BIC_DEC = BICS[1]
        self.BIC_RA_low = BICS[2]
        self.BIC_DEC_low = BICS[3]
        self.b_ra = b_ra[:]
        self.b_dec = b_dec[:]

        w = np.where(self.goodMatches)
        self.std_ra = np.std(self.dra[w])
        self.std_dec = np.std(self.ddec[w])


        self._resid4Panel()

    def _solveMatrix(self):

        w = np.where(self.goodMatches)

        #solve for the solution coefficients
        X = np.copy(self.X[w])
        Y = np.copy(self.Y[w])
        RA = np.copy(self.RA[w])
        DEC = np.copy(self.DEC[w])

        X2 = X*X
        X3 = X2*X
        Y2 = Y*Y
        Y3 = Y2*Y
        XY = X*Y

        #cubed order
        A = np.zeros( (8, len(RA)) ).astype('float64')
        A[0,:] = 1.0
        A[1,:] = X
        A[2,:] = Y
        A[3,:] = X2
        A[4,:] = XY
        A[5,:] = Y2
        A[6,:] = X3
        A[7,:] = Y3
        At = A
        A = At.T

        AtA = np.dot(At,A)
        AtAi = linalg.inv(AtA)
        AtAiAt = np.dot(AtAi,At)
        b_ra = np.dot(AtAiAt,RA)
        b_dec = np.dot(AtAiAt,DEC)


        #squared order
        A = np.zeros( (6, len(RA)) ).astype('float64')
        A[0,:] = 1.0
        A[1,:] = X
        A[2,:] = Y
        A[3,:] = X2
        A[4,:] = XY
        A[5,:] = Y2
        At = A
        A = At.T

        AtA = np.dot(At,A)
        AtAi = linalg.inv(AtA)
        AtAiAt = np.dot(AtAi,At)
        b_ra_low = np.dot(AtAiAt,RA)
        b_dec_low = np.dot(AtAiAt,DEC)



        #now get the ra/dec predictions
        X = np.copy(self.X)
        Y = np.copy(self.Y)
        RA = np.copy(self.RA)
        DEC = np.copy(self.DEC)

        X2 = X*X
        X3 = X2*X
        Y2 = Y*Y
        Y3 = Y2*Y
        XY = X*Y

        #cubed order
        A = np.zeros( (8, len(RA)) ).astype('float64')
        A[0,:] = 1.0
        A[1,:] = X
        A[2,:] = Y
        A[3,:] = X2
        A[4,:] = XY
        A[5,:] = Y2
        A[6,:] = X3
        A[7,:] = Y3
        A = A.T

        predRA = np.dot(A,b_ra)
        predDEC = np.dot(A,b_dec)


        #squared order
        A = np.zeros( (6, len(RA)) ).astype('float64')
        A[0,:] = 1.0
        A[1,:] = X
        A[2,:] = Y
        A[3,:] = X2
        A[4,:] = XY
        A[5,:] = Y2
        A = A.T

        predRA_low = np.dot(A,b_ra_low)
        predDEC_low = np.dot(A,b_dec_low)



        #these seem in error
        lnL_RA = np.sum( -3600**2*(RA[w]-predRA[w])**2/(2*self.dX2[w])) - np.sum(0.5*np.log((2*np.pi*self.dX2[w])))
        lnL_DEC = np.sum( - 3600**2*(DEC[w] - predDEC[w])**2/(2*self.dY2[w]))  - np.sum(0.5*np.log((2*np.pi*self.dY2[w])))

        lnL_RA_low = np.sum( -3600**2*(RA[w]-predRA_low[w])**2/(2*self.dX2[w])) - np.sum(0.5*np.log((2*np.pi*self.dX2[w])))
        lnL_DEC_low = np.sum( - 3600**2*(DEC[w] - predDEC_low[w])**2/(2*self.dY2[w]))  - np.sum(0.5*np.log((2*np.pi*self.dY2[w])))

        BIC_RA = -2*lnL_RA + np.log(len(w[0]))*8
        BIC_RA_low = -2*lnL_RA_low + np.log(len(w[0]))*6
        BIC_DEC = -2*lnL_DEC + np.log(len(w[0]))*8
        BIC_DEC_low = -2*lnL_DEC_low + np.log(len(w[0]))*6
        #print lnL_RA,lnL_DEC
        #print lnL_RA_low,lnL_DEC_low

        return (A, b_ra, b_dec, predRA, predDEC, (lnL_RA,lnL_RA_low), (lnL_DEC,lnL_DEC_low), (BIC_RA,BIC_DEC,BIC_RA_low,BIC_DEC_low))


    def xy2sky(self,xy):
        x = xy[:,0]
        y = xy[:,1]

        x2 = x*x
        x3 = x2*x
        y2 = y*y
        y3 = y2*y
        xy = x*y

        A = np.ones((8, len(x2))).astype('float64')
        A[1,:] = x
        A[2,:] = y
        A[3,:] = x2
        A[4,:] = xy
        A[5,:] = y2
        A[6,:] = x3
        A[7,:] = y3
        A=A.T

        ra = np.dot(A,self.b_ra)
        dec = np.dot(A,self.b_dec)
        coords = np.zeros((len(ra),2)).astype('float64')
        coords [:,0] = ra
        coords [:,1] = dec
        return coords









if __name__ == "__main__":

    import pickle
    from optparse import OptionParser

    parser = OptionParser()
    parser.add_option('--minBA', default = 0.85,
                      type = float, dest = 'minBA', action = 'store',
                      help = 'Minimum B/A sextractor roundness shape parameter to call a star a star. Lower this a small amount if you have trouble finding stars. DEFAULT = %default')
    parser.add_option('--windowSize', default = 13,
                      type = int, dest = 'window', action = 'store',
                      help = 'Window size for plots. DEFAULT = %default')
    parser.add_option('--tsvPath', default = '/Users/fraserw/Documents/Queens Documents/4th Year Projects/Harkness-1617/MatchFiles/',
                      dest = 'tsvPath', type = 'str', action = 'store',
                      help ='Path to the Orcus PS1 tsv file. DEFAULT = %default')
    (opt,args) = parser.parse_args()


    if len(args)>0:
        imagefn = args[0]
    else:
        imagefn = 'EFOSC_Image010_0111_corr.fits'

    with fits.open(imagefn) as han:
        header = han[0].header
        imdata = han[0].data
    header.set('RADESYSa','FK5')

    #scratch crap pixel value fix
    mean = np.median(imdata[::3,::3])
    w = np.where(imdata<0)
    imdata[w] = mean
    HDU = fits.PrimaryHDU(imdata,header)
    List = fits.HDUList([HDU])
    List.writeto('sex.fits',clobber = True)


    #run sextractor
    overwriteSexFiles = True
    if overwriteSexFiles:
        os.system('rm OV.sex default.conv def.param')
    if not path.isfile('OV.sex') or overwriteSexFiles:
        scamp.makeParFiles.writeSex('OV.sex',
                                    minArea=4.,
                                    threshold=5.,
                                    zpt=26.2,
                                    aperture=8.,
                                    min_radius=2.0,
                                    catalogType='FITS_LDAC',
                                    saturate=55000)
    if not path.isfile('default.conv') or overwriteSexFiles:
        scamp.makeParFiles.writeConv()
    if not path.isfile('def.param') or overwriteSexFiles:
        scamp.makeParFiles.writeParam(numAps=1) #numAps is thenumber of apertures that you want to use. Here we use 1

    scamp.runSex('OV.sex', 'sex.fits' ,options={'CATALOG_NAME':'OV.cat'}, verbose=True)
    catalog = trimCatalog(scamp.getCatalog('OV.cat',paramFile='def.param'), minBA = opt.minBA)

    """
    for i in range(len(catalog['XWIN_IMAGE'])):
        print catalog['XWIN_IMAGE'][i],catalog['YWIN_IMAGE'][i]
        if abs(catalog['XWIN_IMAGE'][i]-63)<5 and abs(catalog['YWIN_IMAGE'][i] - 934)<5:
            for k in catalog.keys():
                print k, catalog[k][i]
    sys.exit()

    """


        #


    imSources = []
    for i in range(len(catalog['XWIN_IMAGE'])):
        imSources.append([catalog['XWIN_IMAGE'][i],
                          catalog['YWIN_IMAGE'][i],
                          catalog['MAG_AUTO'][i]+2.5*np.log10(header['EXPTIME']),
                          catalog['FLUX_AUTO'][i],
                          catalog['ERRX2WIN_IMAGE'][i]*0.12**2,
                          catalog['ERRY2WIN_IMAGE'][i]*0.12**2])
    imSources = np.array(imSources)


    #stellar g-r colour range for wcs matching
    #orcus g-r ~0.44
    starColourRange = [0.0,1.0]


    #load up stellar colours and mjds.
    with open('{}/Orcus_PS.tsv'.format(opt.tsvPath)) as han:
        tsv = han.readlines()

    tsvStars = []
    for i in range(1,len(tsv)):
        s = tsv[i].split()
        tsvStars.append([float(s[1]),
                         float(s[2]),
                         float(s[3]),
                         float(s[5])])
    tsvStars = np.array(tsvStars)
    w = np.where((tsvStars[:,2]-tsvStars[:,3]>starColourRange[0])&(tsvStars[:,2]-tsvStars[:,3]<starColourRange[1]))
    tsvStars = tsvStars[w]


    ms_wcs = matrixWCSSolver(imdata, header, imSources, tsvStars, windowSize = opt.window)
    ms_wcs.initialMatch()
    #with open('test.pickle','w+') as han:
    #    pickle.dump(ms_wcs.matches,han)

    #with open('test.pickle') as han:
    #    matches = pickle.load(han)
    #    ms_wcs.matches = matches[:]

    ms_wcs.solveMatrix()

    for i in range(len(ms_wcs.b_ra)):
        header.update('B_RA_{}'.format(i),ms_wcs.b_ra[i])
        header.update('B_DEC_{}'.format(i),ms_wcs.b_dec[i])

    header.update('STD_RA',ms_wcs.std_ra)
    header.update('STD_DEC',ms_wcs.std_dec)
    header.update('NMATCH',len(np.where(ms_wcs.goodMatches==1)[0]))

    HDU = fits.PrimaryHDU(imdata, header)
    List = fits.HDUList([HDU])
    List.writeto('{}_wcs.fits'.format(imagefn.replace('.fits','')), clobber = True)

    #print ms_wcs.xy2sky(np.array([[336.81,  595.00]]))
