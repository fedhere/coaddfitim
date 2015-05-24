import os,sys,glob
import optparse
import pylab as pl
import numpy as np

import pyfits as PF
import Image


from coaddutils import *
from coaddimconfig import *
print keys

SAFE=1
colors={0:'R',1:'G',2:'B'}


if __name__=='__main__':
    parser = optparse.OptionParser(usage="coaddim.py <path to images>", conflict_handler="resolve")
    parser.add_option('-i','--interactive', default=False, action="store_true",
                      help='interactive plotting (with ion) and option to reject the last coadded image')
    parser.add_option('--noalign', default=False, action="store_true",
                      help='just stack images without aligning')
    parser.add_option('--norescale', default=False, action="store_true",
                      help='do not rescale image: otherwise the BITSCALE is set to 16')
    parser.add_option('--nosatmask', default=False, action="store_true",
                      help='do not mask saturation')
    parser.add_option('--mask', default=None, type="string",
                      help='selecting a restricted region of the image to be used for cross correlation: x0,y0,x1,y1')
    parser.add_option('-s','--saturation', default=None, type="float",
                      help='the saturation value')
    parser.add_option('--saturationkey', default=None, type="float",
                      help='the saturation header keyword')
    parser.add_option('--exposure', default=None, type="float",
                      help='the exposure')
    parser.add_option('--exposurekey', default=None, type="string",
                      help='the exposure header keyword')
    parser.add_option('--objectkey', default=None, type="string",
                      help='the target object header keyword')
    parser.add_option('--object', default=None, type="string",
                      help='the target object')
    parser.add_option('--date', default=None, type="string",
                      help='the date (in format like : 2000-01-01)')
    parser.add_option('--datekey', default=None, type="string",
                      help='the date header keyword')
    parser.add_option('--filter', default=None, type="string",
                      help='the filter')
    parser.add_option('--filterkey', default=None, type="string",
                      help='the filter header keyword')
    parser.add_option('--mjd', default=None, type="float",
                      help='the mjd')
    parser.add_option('--mjdkey', default=None, type="string",
                      help='the mjd header keyword')
    parser.add_option('--readnoise', default=None, type="float",
                      help='the readnoise')
    parser.add_option('--readnoisekey', default=None, type="string",
                      help='the readnoise header keyword')
    parser.add_option('--gain', default=None, type="float",
                      help='the gain')
    parser.add_option('--gainkey', default=None, type="string",
                      help='the gain header keyword')
    parser.add_option('--telescope', default=None, type="string",
                      help='the telescope')
    parser.add_option('--telescopekey', default=None, type="string",
                      help='the telescope name header keyword')
    parser.add_option('--impath', default='./', type="string",
                      help='path to images')
    parser.add_option('--jpg',default=False,action="store_true",
                       help="using jpg's (or png, tif...)")
    parser.add_option('--color',default=1,type=int,
                       help="RGB color channel for color images")
    parser.add_option('--swarp', default=False, action="store_true",
                      help='run swarp command as well')
    parser.add_option('--onlyshow', default=False, action="store_true",
                      help='just check keywords and show output names')
    parser.add_option('-v','--verbose', default=False, action="store_true",
                      help='verbose mode')
    
    options,  args = parser.parse_args()
    print options,args

    if len(args)==0:
        options,  args = parser.parse_args(args=['--help'])
        sys.exit(0)
    SHOW=False
    pl.ion()
    if options.interactive:
        pl.ion()
        SHOW=True
    pixmask=[0,0,-1,-1]
    if options.mask:
        pixmask=[int(m) for m in options.mask.split(',')]

    mysaturate= options.saturation    
    if not mysaturate and options.saturationkey:
        keys['saturate']=options.saturationkey
    myobject= options.object
    if not myobject and options.objectkey:
        keys['object']=options.objectkey
    mydate= options.date
    if not mydate and options.datekey:
        keys['date']=options.datekey
    myscope= options.telescope
    if not myscope and options.telescopekey:
        keys['telescope']=options.telescopekey
    myfilter= options.filter
    if not myfilter and options.filterkey:
        keys['filter']=options.filterkey
    mymjd= options.mjd
    if not mymjd and options.mjdkey:
        keys['mjd']=options.mjdkey
    mygain= options.gain
    if not mygain and options.gainkey:
        keys['gain']=options.gainkey
    myrdnoise= options.readnoise
    if not myrdnoise and options.readnoisekey:
        keys['rdnoise']=options.readnoisekey
    myexposure= options.exposure
    if not myexposure and options.exposurekey:
        keys['exposure']=options.exposurekey


    impath=options.impath
#    print [glob.glob(impath+'/'+arg) for arg in args]
    allimgs=[glob.glob(impath+'/'+arg) for arg in args]
    print allimgs
    ##trick to flatten the list of images and remove possible duplicates
    try :allimgs=[val for sublist in allimgs for val in sublist]
#sum(allimgs, []) 
    except TypeError: allimgs=allimgs[0]
    print allimgs
    allimgs=sorted(set(allimgs),key=allimgs.index) 
    inpath=allimgs[0]    
    
    if options.jpg:
        im=Image.open(inpath)
        tmp=np.array(im)
        tmp=tmp[:,:,options.color]
        hdu =PF.PrimaryHDU(tmp)
        header0=hdu.header
        hdu.writeto(inpath+'_'+colors[options.color]+'.fits',clobber=True)
    else:
        image0=PF.open(inpath)
        tmp=PF.getdata(allimgs[0])
        header0=image0[0].header
    nfiles=len(allimgs)
#    allimgs=np.array(allimgs)[-1:0:-1]

    for k,value in keys.items():
        print "looking for :",value
        if k=='exposure':
            if not myexposure:
                try:
                    print header0[value]
                except KeyError:
                    print "missing required keyword:",k,value
                    sys.exit()
                myexposure = float(header0[value])
            print "exposure :",myexposure

        elif k=='saturate' :
            if not mysaturate:
                try:
                    print header0[value]
                    mysaturate = header0[value]
                except KeyError:
                    print "missing required keyword ",k,value
                    sys.exit()
            print "saturate :",mysaturate

        elif k == 'gain':
            if not mygain:
                try:
                    print header0[value]
                    mygain = header0[value]
                except KeyError:
                    print "missing required keyword ",k,value
                    sys.exit()
            print "gain :",mygain

        elif k == 'rdnoise':
            if not myrdnoise:
                try:
                    print header0[value]
                    myrdnoise = header0[value]
                except KeyError:
                    print "missing required keyword ",k,value
                    sys.exit()
            print "read noise :",myrdnoise

        elif k=='object':
            if not myobject:
                try:
                    print header0[value]
                    myobject = header0[value]
                except KeyError:
                    print "missing required keyword ",k,value
                    myobject='MYSTERYOBJECT'
            print "object:",myobject

        elif k=='mjd':
            if not mymjd:
                try:
                    print header0[value]
                    mymjd = header0[value]
                except KeyError:
                    print "missing keyword ",k,value
                    mymjd=0.0
            print "mjd :",mymjd

        elif k=='date':
            if not mydate:
                try:
                    print header0[value]
                    mydate = header0[value]
                except KeyError:
                    print "missing keyword ",k,value
                    import datetime
                    mydate=str(datetime.date.today())
            print "date:",mydate

        elif k == 'telescope':
            if not myscope:
                try:
                    print header0[value]
                    myscope=header0[value]                       
                except KeyError:
                    print "missing keyword ",k,value
                    myscope='SOMESCOPE'
            print "telescope: ",myscope

        elif k == 'filter':
            if not myfilter:
                try:
                    print header0[value]
                    myfilter = header0[value]
                except KeyError:
                    print "missing keyword ",k,value
                    myfilter='UNKNOWN'
            print "filter: ",myfilter



        

    ##check for date format:
    if 'T' in mydate or 'UT' in mydate:
        mydate=mydate[:10]
    nodate=False
    assert( mydate[4]=='-' and mydate[7]=='-'), "unknown format for date. proper format is YYYY-MM-DD"
    nodate=True

    print "number of files:", nfiles
    exposure = (nfiles)*myexposure
    print "first image:", inpath, args
    print "stacked exposure: ",exposure
    print '''WARNING: for now all files must have the same exposure 
for this value to be correct. 
the value will enter the name of the output file, but not the header'''





    outfile=impath+myobject.replace(' ','_')+'_'+mydate.replace('-','')+'_'+myscope+'_%d'%int(exposure)+'sec_'+myfilter+'_%02dstack.fits'%nfiles 
    print "outfie:", outfile
    print options.onlyshow
    if not options.onlyshow:
        PROCESS=True
    else:
        sys.exit()
    if PROCESS:
        imtype='COADDED'
        if options.noalign:
            align='UNALIGN'
            ALIGN=False
        else:
            align='ALIGN'
            ALIGN=True
        print "Generating ", imtype
        sys.stdout.flush()
        
        if options.swarp:
            cmd = "swarp "
            for i in args:
                cmd = cmd + " " + i
                print cmd
            outfile=impath+myobject.replace(' ','_')+'_'+mydate.replace('-','')+'_'+myscope+'_%d'%int(exposure)+'sec_'+myfilter+'_swarp.fits' 
            
            cmd = "mv coadd.fits "+outfile
            
            os.system(cmd)
            
            cmd = '~/etadata/FITS2jpeg/fits2jpeg -fits '+ outfile +' -jpeg '+outfile.replace('.fits','.jpg')+' -nonLinear'
            print cmd
            os.system(cmd)
        if SHOW:
            print "original image here we go!"
            pl.figure(1)
            pl.imshow(tmp)
            pl.title("original")
            pl.draw()
            inrange=raw_input("input new range as x0,y0,x1,y1\n")
            if not inrange:
                x0,y0,x1,y1=0,0,-1,-1
            else:
                x0,y0,x1,y1=[int(ir) for ir in inrange.split(',')]

        nimgs=1
        acceptedlist=[allimgs[0].split("/")[-1]] 
        if options.mask:
            print pixmask[0],pixmask[2],pixmask[1],pixmask[3]
            tmp=tmp[pixmask[0]:pixmask[2],pixmask[1]:pixmask[3]]
        if ALIGN:
            for i in range(1,nfiles):
                print "working on",allimgs[i]
                if options.jpg:
                    image=Image.open(allimgs[i])
                    tmp2=np.array(image)
                    tmp2=tmp2[:,:,options.color][pixmask[0]:pixmask[2],pixmask[1]:pixmask[3]]
                    hdu = PF.PrimaryHDU(tmp2)
                    header=hdu.header
                    hdu.writeto(allimgs[i]+'_'+colors[options.color]+'.fits',clobber=True)
                else:
                    image=PF.open(allimgs[i])
                    tmp2= PF.getdata(allimgs[i])[pixmask[0]:pixmask[2],pixmask[1]:pixmask[3]]
                    header=image[0].header

                (t1,t2)=correlate([tmp,tmp2],mysaturate, showme=SHOW)

                if  t1==None and t2==None:
                    print "image correlation failed"
                    continue
                print t1.shape
                print t2.shape

                #         s2= header['STDDEV']
                #w2=1.0/s1#sigclippedstd(t2)
                #print w2
                #w1=1.0/s2#sigclippedstd(t1)
                #print w1
                #w1=w1*w1
                #w2=w2*w2	
                #ratio = (w1+w2)

                w1,w2,ratio=1,1,1
                newtmp=(t1*w1+t2*w2)/ratio
                acc=''
                if SHOW:
                    print "image cut:", x0,y0,x1,y1
                    pl.figure(3)
                    pl.clf()
                    pl.imshow(newtmp[x0:x1,y0:y1])#+80:750+160,750+160:750+260])
                    pl.title("stack %d"%i)
                    pl.draw()
                    acc=raw_input("accept this image? Y/n")
                    print acc
                    if not acc.lower().startswith('n'):
                        print "accepting this image and phase"
                        tmp=newtmp
                        acceptedlist.append(allimgs[i].split("/")[-1])
                        try:
                            myexposure += float(header[keys['exposure']])
                        except KeyError:
                            myexposure+=exposure
                        nimgs+=1
                    else:
                        pl.imshow(tmp[x0:x1,y0:y1])
                        pl.title("stack %d"%i)
                        pl.draw()
                else:
                    tmp=newtmp
                    acceptedlist.append(allimgs[i].split("/")[-1])
                    try:
                        myexposure += float(header[keys['exposure']])
                    except KeyError:
                        myexposure+=exposure
                    nimgs+=1

                #print tmp[0,0:100]
                    #         s1=sqrt(s1*s1+s2*s2)
    else:
        for i in range(1,nfiles):
            image=PF.open(allimgs[i])
            header=image[0].header
            tmp+=PF.getdata(allimgs[i])
            try:
                myexposure += float(header[keys['exposure']])
            except KeyError:
                myexposure=exposure
            nimgs+=1
            acceptedlist.append(allimgs[i].split("/")[-1])
    pl.imshow(tmp)
    pl.title("final image")
    print mymjd,type(mymjd), int(float(mymjd))
    mymjd=float(mymjd)
    imid=int((mymjd-int(mymjd))*1e4)
    print imid
    tmp/=nimgs
    effgain=float(mygain)*nimgs
    if not options.nosatmask:
        tmp=masksaturated(tmp, header0,mysaturate)
    if len(tmp)> 0:
        median = np.median(tmp)
    else:
        median = 0
            #print 'median flux:',median
            #tmp = where(tmp<median*0.7,SATURATION,tmp)
            #tmp = where(tmp>SATURATION,SATURATION,tmp)#
            
    if len(maskcenters)>0:
        tmp=fixim(tmp,maskcenters)
    if myexposure != exposure:
        print "some image must have failed or something.... total exposure shorter than single exposure x number of files..."
        
        #####     UPDATE HEADER     #####

    header0[keys['exposure']]=myexposure/nimgs
    header0.update('MJD-OBS','%f' % mymjd)
    readnoise=float(myrdnoise)#/sqrt(nimgs)
    readnoise*=np.sqrt(nimgs)
    readnoise=float(readnoise)
    header0.update('IMAGTYPE', '%s' % imtype, 'Coadded (summed stack)')
    header0.update('NCOADDED', '%d' % nfiles)
    header0.update(keys['gain'], '%f' % effgain, 'Gain, corrected for stacking')
    header0.update(keys['rdnoise'], '%f' % readnoise, '[e] Readout noise corrected for stacking')
    header0.update('SATURATE',mysaturate)
    header0.update('OBSTYPE','EXPOSE')
    if ALIGN:
        header0.update('ALIGN', '%s' % imtype, 'Aligned')
        for ii,img in enumerate(acceptedlist):
            print ii,img
            header0.update('INIMG%02d'%ii, 'input image %d'%ii)
            header0.comments['INIMG%02d'%ii]=img
    else:
        header0.update('ALIGN', '%s' % imtype, 'Unaligned')

    base=myobject.replace(' ','_')+'_'+mydate.replace('-','')+'_'+myscope+'_%d'%int(myexposure)+'sec_'+myfilter+'_coadd'
    outfile=impath+base+'.fits' 

    if SAFE:
        if os.path.isfile(outfile):
            rm = raw_input("output file exists: remove? (Y/n)\n")
            if rm == '' or rm.lower().startswith('y'):
                os.remove(outfile)
            else:
                if base+'_1.fits' not in os.listdir(impath) :
                    outfile = impath+base+'_1.fits' 
                else :
                    existing=sorted([f for f in os.listdir(impath) if base in f and f.endswith('.fits')])
                    try: int(existing.split('_')[-1].split('.')[0])
                    except: 
                        existing=existing[1:]
                    outfile = base+'_' + str(np.max([int(ii.split('_')[-1].split('.')[0]) for ii in existing])  + 1) + '.fits'


    out_fits = PF.PrimaryHDU(header=header0,data=tmp)
    if not options.norescale: out_fits.scale(type=out_fits.NumCode[16],bzero=32768.0,bscale=1.0)
    if options.jpg:
        outfile=outfile[:-5]+'_'+colors[options.color]+'.fits'
    print 'writing',outfile
    if SAFE:
        out_fits.writeto(outfile, clobber=False)
    else: 
        out_fits.writeto(outfile, clobber=True)
#PF.writeto(outfile, tmp, header0)

    print 'Composite image written to %s' %  outfile


