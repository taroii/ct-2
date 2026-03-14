#########################################
#FDK for flat detector
#########################################
import numpy as n
import sys
#sys.path.append('/home/seunghee/utils.Utils.Python')
import setFreqrange
pi = n.pi

def fdk_sino_weight3D(radius, center2detector, ns, vu, vv, sinomat):
	##########################################################
	#weighting
	#where F is a distance from source to center of rotation and
	#D is a distance from detector to center of rotation	
	#vu = vector u
	#vv = vector v
	#ns = number of projection angles s
	##########################################################
	F = radius
	D = center2detector
	weight = F/n.sqrt((F+D)**2 + vv**2 + vu**2)
	w_sinomat = n.tile(weight,[ns,1,1]) * sinomat
	return w_sinomat

def fdk_sino_filter3D(dp, dim, sinomat, opt_window='', L=1):
	##############################################################
	#apply ramp filter using discrete implementation described in Kak & Slaney
	#dp = sampling interval along the chosen dimension
	#dim = dimension of the sinogram to be filtered
	#L = scaling parameter to set cutoff frequency when additional window is used,
	#    typically, 0.5, 0.75, or 1.0
	#############################################################
	#Zeropadding
	origsize = sinomat.shape[dim]
	#npad = (2**(n.ceil(n.log2(2*sinomat.shape[dim] -1)))).astype(n.integer)
	npad = (2**(n.round(n.log2(2*sinomat.shape[dim] -1)))).astype(n.integer)
	
	#Ramp filter
	N = npad
	if N % 2 ==0:
		nn=n.arange(-N/2,N/2, dtype=float)
	else:
		nn=n.arange(-(N-1)/2,(N-1)/2+1, dtype=float)
	hn=n.zeros(nn.size)
	hn[npad/2] = 1./(4*dp**2) #h[0]
	for ii in range(hn.size):
		if abs(nn[ii]) % 2 == 1: #h[odd]
			hn[ii]=-1./(pi*nn[ii]*dp)**2
	hn = n.fft.ifftshift(hn)
	Hk = n.real(n.fft.fft(hn))
	Hk = dp * Hk; # dp from discrete implementation of convolution integral
	ramp = Hk
	
	#Hanning window: L changes the cutoff frequency of the window
	if opt_window == 'Hann':		
		kappa = setFreqrange.kappa1D(npad,dp)
		cutoff_freq = 1/(2*dp)*L
		
		if L < 1.0:
			#find kappa closest to cutoff_freq
			diff = n.abs(n.abs(kappa) - cutoff_freq)
			mindiff = diff.min()
			kappaindex = n.nonzero(diff == mindiff)[0]
			startindex = kappaindex[0]
			endindex = kappaindex[-1]
			print 'Hanning window starts from freq=%g to freq=%g'%(kappa[startindex],kappa[endindex])
			HannFilt = n.zeros(kappa.shape)
			HannFilt[startindex:endindex+1] = 0.5 + \
			                0.5*n.cos(n.pi*kappa[startindex:endindex+1]/cutoff_freq)
			HannFilt = n.fft.ifftshift(HannFilt)
		else: #L == 1.0
			HannFilt = n.fft.ifftshift(0.5 + 0.5 *n.cos(n.pi*kappa/cutoff_freq))
		Filt = ramp * HannFilt
	else:
		Filt=ramp
	
	#initialize
	if dim == 1: #filter along u axis
		sinofiltmat = n.zeros((sinomat.shape[0], npad, sinomat.shape[2]))
		zeropad = n.zeros(npad-sinomat.shape[1])
	else: #filter along v axis
		sinofiltmat = n.zeros((sinomat.shape[0], sinomat.shape[1], npad))
		zeropad = n.zeros(npad-sinomat.shape[2])
	
	#compute filtered sinogram
	for ii in range(sinofiltmat.shape[0]): #ns
		if dim == 1: #filtering along u axis
			for jj in range(sinofiltmat.shape[2]):  #nv
				Fline = n.fft.fft(n.hstack((sinomat[ii,:,jj], zeropad)))
				Fline = Filt*(Fline)
				sinofiltmat[ii,:,jj]=n.real(n.fft.ifft(Fline))
		else: #filtering along v axis
			for jj in range(sinofiltmat.shape[1]):  #nu
				Fline = n.fft.fft(n.hstack((sinomat[ii,jj,:], zeropad)))
				Fline = Filt*(Fline)
				sinofiltmat[ii,jj,:]=n.real(n.fft.ifft(Fline))
	
	if dim == 1: #along u axis
		sinofiltmat = sinofiltmat[:,0:origsize,:]
	else: #along v axis
		sinofiltmat = sinofiltmat[:,:, 0:origsize]
	
	return sinofiltmat

def setcoordinates(nx, ny, nz, x0, y0, z0, \
	dx, dy, dz, xoffset, yoffset, zoffset, mask_flag, \
	umax, vmax, radius, source2detector):
	
	cubevox = n.mgrid[0:nz, 0:nx, 0:ny]
	vz = z0 + (cubevox[0,:,:] + zoffset) * dz
	vx = x0 + (cubevox[1,:,:] + xoffset) * dx
	vy = y0 + (cubevox[2,:,:] + yoffset) * dy
	
	#masking
	mask = n.ones(vx.shape)
	if mask_flag == 1:
		gammamax = n.arctan2(n.sqrt(umax**2 + vmax**2),source2detector)
		rmax = n.sin(gammamax) * radius
		rr = n.sqrt(vx**2 + vy**2 + vz **2)
	mask = mask * (rr < rmax)
	#output of numpy.nonzero: row comes first
	vx = vx[n.nonzero(mask)]
	vy = vy[n.nonzero(mask)]
	vz = vz[n.nonzero(mask)]
	
	xyz = n.vstack((vx,vy,vz))
	
	return xyz

def fdk_backproj_linear3D(ns, nu, nv, du, dv, u0, v0, \
	radius, source2detector, \
	uoffset, voffset, \
	vs, vx, vy, vz, sinomat):
	##################################################################
	#compute backprojection
        #ns = number of projection angles s
	#nu = number of pixels along the u axis
	#nv = number of pixels along the v axis
	#du = sampling interval along the u axis in the sinogram
	#du = sampling interval along the v axis in the sinogram
	#u0 = value of u corresponding to the first pixel
	#v0 = value of v corresponding to the first pixel
	#vs = vector s of the sinogram
	#vx = vector x of the original image
	#vy = vector y of the original image
	#sino = sinogram	
	#zoomfac = zoom factor based on the geometry where F is a distance 
	#	from source to center of rotation and
	#	D is a distance from detector to center of rotation
	##################################################################	
	F = radius
	D = source2detector - radius
	
	zoomfac = n.zeros(vx.shape)
	
	map_backproj_leftu = n.zeros(vx.shape).astype(n.integer) #left index of u
	map_backproj_leftv = n.zeros(vx.shape).astype(n.integer) #left index of v
	uinlen = n.zeros(vx.shape) # u in length unit: needed for calculating weights in interpolation
	vinlen = n.zeros(vx.shape) # v in length unit: needed for calculating weights in interpolation
		
	#back projection based on linear interpolation
	recon = n.zeros((vx.size, vz.size))
	for jj in range(vz.size):
		for ii in range(ns):
			#compute magnification ratio based on geometry
			zoomfac = (F+D)/(F -vx*n.cos(vs[ii]) - vy*n.sin(vs[ii]))
			uinlen = (-vx*n.sin(vs[ii]) + vy*n.cos(vs[ii])) * zoomfac
			vinlen = vz[jj] * zoomfac
			
			map_backproj_leftu = n.floor((uinlen - u0)/du - uoffset).astype(n.integer) #left bin along u axis
			map_backproj_leftv = n.floor((vinlen - v0)/dv - voffset).astype(n.integer) #left bin along v axis
			
			tmpsino = sinomat[ii,:,:]
			
			#mapping index >=0 based on the left index and index < nu based on the right index
			mapping = n.ones(vx.shape)
			mapping = mapping * (map_backproj_leftu+1 < nu ) * (map_backproj_leftu >= 0) * \
			  (map_backproj_leftv+1 < nv ) * (map_backproj_leftv >= 0)
			mindex = n.nonzero(mapping)
				
			#along u axis at v=left index
			tmpsino_xy_u1v1 = tmpsino[map_backproj_leftu[mindex],map_backproj_leftv[mindex]]
			tmpsino_xy_u2v1= tmpsino[map_backproj_leftu[mindex]+1,map_backproj_leftv[mindex]]
			#along u axis at v=right index
			tmpsino_xy_u1v2 = tmpsino[map_backproj_leftu[mindex],map_backproj_leftv[mindex]+1]
			tmpsino_xy_u2v2= tmpsino[map_backproj_leftu[mindex]+1,map_backproj_leftv[mindex]+1]
			
			#weights along u axis
			weight_rightu = (uinlen[mindex] - u0)/du - (map_backproj_leftu[mindex] + uoffset)
			weight_leftu = 1 - weight_rightu
			
			#weights along v axis
			weight_rightv = (vinlen[mindex] - v0)/dv - (map_backproj_leftv[mindex] + voffset)
			weight_leftv = 1 - weight_rightv
			
			#interpolate along u axis
			tmpsino_xy_v1 = weight_leftu * tmpsino_xy_u1v1 + weight_rightu * tmpsino_xy_u2v1
			tmpsino_xy_v2 = weight_leftu * tmpsino_xy_u1v2 + weight_rightu * tmpsino_xy_u2v2
			
			#then, interpolate along v axis
			tmpsino_xy = weight_leftv * tmpsino_xy_v1 + weight_rightv * tmpsino_xy_v2
			
			recon[mindex,jj] = recon[mindex,jj] + tmpsino_xy * zoomfac[mindex] * zoomfac[mindex]
		print 'nz: %d'%(jj)
		
	recon = pi/ns*recon; #don't forget the constant from the integral over angle!
	
	return recon
