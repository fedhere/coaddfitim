coaddfitim
==========

this code stacks fits image files. it is able to shift the imagies to achieve best alignment but correlating the images (in fourier space)

it deals with saturated pixels and articats of saturation, it can remove bad or unwanted regions of the fits file, and of course appropriately corrects the exposure, gain, readnoise values in the header. 

run it as 

$python coaddim.py <wild cards for the fits files> 

and to get the help run the command line without arguments:

$python coaddim.py

