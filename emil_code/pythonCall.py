
imnormlib = numpy.ctypeslib.load_library('/home/emil/python/PyTomo3Diterative/imagenorms','.')
clpgrad = imnormlib._Z7lp_gradffPfS_iiijj
ctpvgrad = imnormlib._Z16total_p_var_gradffPfS_iiijj
ctvgrad = imnormlib._Z14total_var_gradfPfS_iiijj


   def tv_grad_cuda(self,nblocks,blocksize,eps=0.000000001):
      '''Compute tv grad on GPU.
'''
        
      nx = ctypes.c_int()
      nx.value=self.mat.shape[0]
      ny = ctypes.c_int()
      ny.value=self.mat.shape[1]
      nz = ctypes.c_int()
      nz.value=self.mat.shape[2]

      peps = ctypes.c_float()
      peps.value=eps

      nblk = ctypes.c_uint()
      nblk.value=nblocks
      blksz = ctypes.c_uint()
      blksz.value=blocksize

      temp = self.mat*1.
      temp.fill(0.)

      ctvgrad(peps,\
                ctypes.c_void_p(temp.ctypes.data),\
                ctypes.c_void_p(self.mat.ctypes.data),\
                nx,ny,nz,\
                nblk,blksz)

      total = sqrt(sum(temp.flatten()**2.,dtype=float64))
      self.mat = 1.*temp/total

      return total

