totalvar:
	f2py -c totalvar.f -m totalvar
tomo3Dfort:
	f2py -c tomo3DfortranAux.f -m tomo3Dfort
ray_proj3D:
	f2py -c ray_proj3D.f -m ray_proj3D
view_proj3D:
	f2py -c view_proj3D.f -m view_proj3D
view_proj3D_tomo:
	f2py -c view_proj3D_tomo.f -m view_proj3D_tomo
radon3D:
	f2py -c radon3D.f -m radon3D
forts: totalvar tomo3Dfort ray_proj3D view_proj3D view_proj3D_tomo radon3D

clean:
	rm *.so *.pyc phantoms/*.pyc
