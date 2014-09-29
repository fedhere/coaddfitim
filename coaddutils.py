
###############################################################
###############written by fbb 2011########################
#coadds images correlating in fft space (if ALIGN is set to 1)
####
#f.b.b
#last modified 10/25/2011
################################################################

from scipy.fftpack import *
import sys,os,time
import pyfits as PF
from numpy import *
import numpy as np 
from scipy import *
from scipy import cluster
from scipy.stats import mode
from scipy.interpolate import interp1d
from pylab import imshow,title,figure,show
import time

#import ds9
SHIFT=1
CUT=0
#CUT=1
RAD=20
'''
SATURATION=27000
ALIGN=1
PROCESS = 1
THRESHOLD=10
SAFE = 0
FIXHEADER=1
FIXSATURATION=0
MASKIM=0

'''
def read2list(filename):

    listfile = open(filename)
    lst = []
    for line in listfile.readlines():
        lst.append(line.strip())

    return lst

def increment_writename(outpath, writename):

    if os.path.isfile(outpath+'/'+writename):run = True
    else: run = False

    num = int((writename[-6:])[:-5])+1
    while run==True:
        writename = writename[:-6]+str(num)+'.fits'
        if os.path.isfile(outpath+'/'+writename):
            num += 1
        else:
            run = False

    return writename
    
def fixim(data):
   #if regions of the ccd need to be masked 
   #because of a damaged chip, or a stadd that is difficult to correct, 
   #they get set to 0 here
   
   print data.shape
   print maskcenters
   n1=data.shape[0]
   n2=data.shape[1]
   x=array(range(n1*n2))
   y=array(range(n1*n2))
   
   if len(maskcenters)==0:
      return data
   x=x%n1
   y=y/n2

   x=reshape(x,(n2,n1))
   y=reshape(y,(n1,n2))

   for c in maskcenters:
      print c
      print x
      print y
      maxsize=min(n1,n2)
      w=where(((x[:maxsize,:maxsize]-c[0][0])**2 + (y[:maxsize,:maxsize]-c[0][1])**2)<c[1]*c[1])
      data[w]=0
      
      
#   imshow(data[maskcenters[0][0][0]-maskcenters[0][1]*5:maskcenters[0][0][0]+maskcenters[0][1]*5,maskcenters[0][0][1]-maskcenters[0][1]*5:maskcenters[0][0][1]+maskcenters[0][1]*5])
#   show()
#   sys.exit()
   return(data)






def dist(pair):
   return abs(pair[0][0]-pair[1][0])+abs(pair[0][1]-pair[0][1])


def sigclippedstd(arr):
   iter = 0
   oldarr=zeros(len(arr.flatten())+1,float)
   while (iter < 20 and mean(arr)>median(arr)): 
#len(arr.flatten())<len(oldarr.flatten())):
    	oldarr=arr
        m=mean(arr)
        s=std(arr)
#	print s, m
	arr=arr.clip(min=m-1*s,max=m+1*s)
	iter+=1
   return(s)	


def mymode(data):
    mydata=data#[where(data<5000)]
    counts = {}
    for x in mydata.flatten():
       counts[x] = counts.get(x,0) + 1
#       print counts[x]
    maxcount = max(counts.values())
    modelist = []
    for x in counts:
       if counts[x] == maxcount:
              modelist.append(x)
    return modelist[0]


def masksaturated(data, header0, saturation):
   #this function looks for saturated pixels 
   # and also fixes the pixels around and inside of saturation regious
   # where saturation occcurs the core of stars counts may be artificially
   # low in some telescopes. that must be fixed 
   #since saturated pixels really mess up the correlation and 
   #someetime saturation occurs BELOW the saturation limit
   n1=data.shape[0]
   n2=data.shape[1]
   index = array(range(data.shape[0]*data.shape[1]))
   print index
   index1=index % data.shape[1]
   index2=index / data.shape[1]
   index1=index1.reshape(data.shape[0],data.shape[1])
   index2=index2.reshape(data.shape[0],data.shape[1])
   
   index=[index1,index2]
   mode= mymode(data)
   wh=where(data<(mode)-mode*0.5)
   
   allsat=[]
   centers = []
   for i in range(20,3,-1):
      mocdata=zeros((data.shape[0],data.shape[1]),float)
      data1=data[:-i,:]
      data2=data[i:,:]
      k1=i/2
      k2=data2.shape[0]
      data3=data[k1:k2+k1,:]

      sat= where ((data1>saturation) & (data2 > saturation) & (data3 <saturation))
      sat=array(zip(*sat))
      
      print "sat" , sat
      if len(sat)<=1:
         continue

      for s in sat:
         print s[1],s[0]
         allsat.append([s[0],s[1]])
         mocdata[s[0],s[1]]=saturation
#         imshow(mocdata)
#         show()
#      print sat


      outfilehere="tmp_%d.fits"%i
      out_fits = PF.PrimaryHDU(header=header0,data=mocdata)
      out_fits.scale(type=out_fits.NumCode[16],bzero=32768.0,bscale=1.0)
      print 'writing',outfilehere
      out_fits.writeto(outfilehere, clobber=True)

   allsat = array(allsat)
#   print allsat
   if len(allsat)==0:
	wh=where(data>saturation)
   	data[wh[0],wh[1]]=saturation
 
	return data

   clusters= cluster.hierarchy.fclusterdata(allsat, 10, criterion='distance', metric='euclidean', depth=2, method='single', R=None)

#   print clusters
   for i in range(max(clusters)): 
      if len(where(clusters == i+1)[0])>1:
         centers.append(mean(allsat[where(clusters==i+1)[0]], axis=0))
         print "centers", mean(allsat[where(clusters==i+1)[0]], axis=0), len(where(clusters==i)[0])
   
   print  data1.shape,data2.shape
   for x in centers:
      print "circle(",x[1],",",x[0],",10)"
      ind=where(sqrt((index1-x[1])**2+(index2-x[0])**2)<RAD)
#      print ind[0],ind[1]
      data[ind[0],ind[1]]=saturation
 
   try:
      data[wh[0],wh[1]]=saturation
   except:
      pass
   try:
      data[wh[0]-1,wh[1]-1]=saturation
   except:
      pass
   try:
      data[wh[0]+1,wh[1]+1]=saturation
   except:
      pass
   try:
      data[wh[0]+1,wh[1]]=saturation
   except:
      pass
   try:
      data[wh[0]+1,wh[1]]=saturation
   except:
      pass
   try:
      data[wh[0]-1,wh[1]]=saturation
   except:
      pass
   try:
      data[wh[0],wh[1]-1]=saturation
   except:
      pass   
   wh=where(data>saturation)
   data[wh[0],wh[1]]=saturation
 
   return data


def correlate(fits,saturation):

    if fits[0].shape[0] >= fits[1].shape[0]:
        larger = 0
        smaller = 1
    else:
        larger = 1
        smaller = 0

    naxis1=fits[smaller].shape[0]
    naxis2=min(fits[smaller].shape[1],fits[larger].shape[1])

#    print fits[0],fits[1]
### passing data to dummy array
    arr = zeros((naxis1,naxis2),float)
    for i in range(arr.shape[0]):
        for j in range(arr.shape[1]):
            arr[i,j] = fits[smaller][i,j]

    fits[smaller] = arr

### cropping larger down to smaller size
    arr = zeros((naxis1,naxis2),float)
    for i in range(arr.shape[0]):
        for j in range(arr.shape[1]):
            arr[i,j] = fits[larger][i,j]
    fits[larger] = arr 

    print "size before finding phase ",fits[smaller].shape
    print "size before finding phase ",fits[larger].shape

    #imshow(fits[0])    
    #figure()
    #imshow(fits[1])
    #show()

####   FINDING PHASE  ####
    if SHIFT:
      print 'Finding phase...'
      if CUT:    
         mask = zeros(shape(fits[larger]))
         mask[1275:1295,470:490]=+1#400:750,1000:1350]+=1
         fft1 = fftn(fits[larger]*mask)
         ifft2 = ifftn(fits[smaller]*mask)
         R = ifftn(fft1*ifft2).real
      else:
         fft1 = fftn(fits[larger])
         ifft2 = ifftn(fits[smaller])
         R = ifftn(fft1*ifft2).real
#         print R
#         imshow(R)
#         show()
      phase = np.where(R == np.max(R))
      print "phase = " + str(phase)
      ### Checks if img_small has negative shift ###
    
      if phase[0] > naxis2/2:
         axis2_shift =  phase[0] - naxis2
      else:
         axis2_shift = phase[0]
         
      if phase[1] > naxis1/2:
         axis1_shift = phase[1] - naxis1
      else:
         axis1_shift = phase[1]

      if axis2_shift >= 0:
         if axis2_shift == 0: axis2_shift = -naxis1
         if axis1_shift >= 0:
            if axis1_shift == 0: axis1_shift = -naxis2
            return(  fits[larger][axis2_shift:,axis1_shift:],
                     fits[smaller][:-axis2_shift,:-axis1_shift])
            #          stack[axis2_shift:,axis1_shift:] += w*fitspad[:-axis2_shift,:-axis1_shift]
         else: #axis1_shift < 0

            return(  fits[larger][axis2_shift:,:-abs(axis1_shift)],
                     fits[smaller][:-axis2_shift,abs(axis1_shift):])
            #         stack[axis2_shift:,:-abs(axis1_shift)] += w*fitspad[:-axis2_shift,abs(axis1_shift):]

      else: #axis2_shift < 0
         if axis1_shift >= 0:
            if axis1_shift == 0: axis1_shift = -naxis1
            
            return(  fits[larger][:-abs(axis2_shift),axis1_shift:],
                     fits[smaller][abs(axis2_shift):,:-axis1_shift])
            #stack[:-abs(axis2_shift),axis1_shift:] += w*fitspad[abs(axis2_shift):,:-axis1_shift]

         else: #axis1_shift < 0
            return(  fits[larger][:-abs(axis2_shift),:-abs(axis1_shift)],
                     fits[smaller][abs(axis2_shift):,abs(axis1_shift):])

            #          stack[:-abs(axis2_shift),:-abs(axis1_shift)] += w*fitspad[abs(axis2_shift):,abs(axis1_shift):]

    else:
      return(fits[larger],fits[smaller])

    print "...shift in NAXIS1 = %d"  %axis2_shift
    print "...shift in NAXIS2 = %d"  %axis1_shift

    print 'img_sml shape... '+ str(fits[smaller][abs(axis2_shift):,abs(axis1_shift):].shape)
    print 'img_lrg shape... '+ str(fits[larger][:-abs(axis2_shift):,:-abs(axis1_shift)].shape)
    

    #coadded = fits[larger][:-abs(axis2_shift),:-abs(axis1_shift)] + fits[smaller][abs(axis2_shift):,abs(axis1_shift):]
#    print slicex1,slicey1,slicex2,slicey2

#    print "larger cut" ,fits[larger][int(slicey1[0]):int(slicey1[1]), int(slicex1[0]):int(slicex1[1])].shape
#    print "smaller cut", fits[smaller][int(slicey2[0]):int(slicey2[1]), int(slicex2[0]):int(slicex2[1])].shape
    
    #coadded= 
    return (fits[larger][int(slicex1[0]):int(slicex1[1]), int(slicey1[0]):int(slicey1[1])],fits[smaller][int(slicex2[0]):int(slicex2[1]), int(slicey2[0]):int(slicey2[1])])
