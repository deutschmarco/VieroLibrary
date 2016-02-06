import numpy as np
import pdb
import os
from astropy.wcs import WCS
from shift import shift_twod
from zero_pad import zero_pad

def cut_thumbnail(image, xpos, ypos, side_pix, error=None):

	if (side_pix % 2 == 0): side_shift = side_pix / 2# even 
	else: side_shift = (side_pix-1) / 2 #odd

	mapsize = np.shape(image)
	if mapsize[0] >= mapsize[1]:
		l2 = mapsize[0]+side_pix * 2.0
	else: 
		l2 = mapsize[1]+side_pix * 2.0

	padded_image = shift_twod(zero_pad(image,l2 = l2),side_shift,side_shift)
	if error != None:
		padded_error = shift_twod(zero_pad(error,l2 = l2),side_shift,side_shift)
		padded_error[np.where(padded_error == 0)] = 1e10

	ngal = len(xpos)
	cube = np.zeros([ngal,side_pix,side_pix])
	cuberr = np.ones([ngal,side_pix,side_pix])
	for i in range(ngal):
		xlo = xpos[i]
		xhi = xpos[i] + side_pix #2.0 * side_pix
		ylo = ypos[i]
		yhi = ypos[i] + side_pix #2.0 * side_pix
		thumbnail = padded_image[xlo:xhi,ylo:yhi]
		if error != None:
			thumberr = padded_error[xlo:xhi,ylo:yhi]
		else: thumberr = 1.0
		cube[i,:,:] = thumbnail
		if error != None:
			cuberr[i,:,:] = thumberr

	return [cube,cuberr]

def stack_radec(image, hd, radec_positions, side=100, error=None, verb=None, zlist=None, zlo=0., zhi=1000.):
	#if error == None: 
	#	error2 = None
	#else: error2 = error

	#if side_arcmin != 100:
	#	side = 

	w = WCS(hd)
	mapsize = np.shape(image)

	ty,tx = w.wcs_world2pix(radec_positions[0], radec_positions[1], 0) 
	if zlist != None:
		ind_keep = np.where((tx >= 0) & (np.floor(tx) < mapsize[0]) & (ty >= 0) & (np.floor(ty) < mapsize[1]) & (zlist >= zlo) & (zlist < zhi))
	else:
		ind_keep = np.where((tx >= 0) & (np.floor(tx) < mapsize[0]) & (ty >= 0) & (np.floor(ty) < mapsize[1]))
	nt0 = np.shape(ind_keep)[1]
	real_x=np.floor(tx[ind_keep]).astype(int)
	real_y=np.floor(ty[ind_keep]).astype(int)

	if verb != None:
		print str(len(real_x)) + ' galaxies in region'
	thumbs_errors = cut_thumbnail(image,real_x,real_y,side,error=error) 

	wthumb = np.sum(thumbs_errors[0]/thumbs_errors[1],axis=0) / np.sum(1.0/thumbs_errors[1],axis=0) 

	return [wthumb, thumbs_errors,ind_keep]