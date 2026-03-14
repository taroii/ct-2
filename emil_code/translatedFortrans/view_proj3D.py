# translated with Google Gemini, not checked yet

import numpy as np
from numba import jit, prange
import math

# ==========================================
# Constants
# ==========================================
ALARGE = 100000000.0
EPS = 0.00000001

# ==========================================
# Geometric Root Finders
# ==========================================

@jit(nopython=True, cache=True)
def root_plane(zs, zb):
    """
    Calculates intersection roots with a plane (z=0).
    """
    troot = np.empty(2, dtype=np.float32)
    
    if abs(zs - zb) < EPS:
        if zs >= 0.0:
            troot[0] = -ALARGE
            troot[1] = ALARGE
        else:
            troot[0] = 1.0
            troot[1] = -1.0
    else:
        root1 = zs / (zs - zb)
        if zs < zb:
            troot[0] = root1
            troot[1] = ALARGE
        else:
            troot[0] = -ALARGE
            troot[1] = root1
    return troot

@jit(nopython=True, cache=True)
def root_ellipsoid(xs, ys, zs, xb, yb, zb, ax, ay, az):
    """
    Calculates intersection roots with an ellipsoid.
    """
    troot = np.empty(2, dtype=np.float32)
    
    ax2, ay2, az2 = ax**2, ay**2, az**2
    
    denom = (ax2*az2*(yb - ys)**2 + 
             ay2*(az2*(xb - xs)**2 + ax2*(zb - zs)**2))

    term1 = 4 * ax2 * ay2 * az2
    term2 = -(az2*(xs*yb - xb*ys)**2)
    term3 = ay2*(az2*(xb - xs)**2 + ax2*(zb - zs)**2 - (xs*zb - xb*zs)**2)
    term4 = ax2*(az2*(yb - ys)**2 - (ys*zb - yb*zs)**2)
    
    numerator_dist = term1 * (term2 + term3 + term4)
    dist2 = numerator_dist / (denom**2)

    avet_num = (ax2*az2*ys*(-yb + ys) + 
                ay2*(az2*xs*(-xb + xs) + ax2*zs*(-zb + zs)))
    avet = avet_num / denom

    if dist2 < 0.0:
        troot[0] = 1.0
        troot[1] = -1.0
    else:
        sqrt_dist = math.sqrt(dist2)
        troot[0] = avet - 0.5 * sqrt_dist
        troot[1] = avet + 0.5 * sqrt_dist
        
    return troot

@jit(nopython=True, cache=True)
def root_cylinder(xs, ys, xb, yb, ax, ay):
    """
    Calculates intersection roots with a cylinder along Z.
    """
    troot = np.empty(2, dtype=np.float32)
    
    ax2, ay2 = ax**2, ay**2
    
    num = 4 * ax2 * ay2 * (ay2 * (xb - xs)**2 + ax2 * (yb - ys)**2 - (xs*yb - xb*ys)**2)
    denom = (ay2 * (xb - xs)**2 + ax2 * (yb - ys)**2)**2
    
    dist2 = num / denom
    
    avet = (ay2 * xs * (xs - xb) + ax2 * ys * (ys - yb)) / \
           (ay2 * (xb - xs)**2 + ax2 * (yb - ys)**2)

    if dist2 < 0.0:
        troot[0] = 1.0
        troot[1] = -1.0
    else:
        sqrt_dist = math.sqrt(dist2)
        troot[0] = avet - 0.5 * sqrt_dist
        troot[1] = avet + 0.5 * sqrt_dist
    return troot

@jit(nopython=True, cache=True)
def root_cone(xs, ys, zs, xb, yb, zb, ax, ay):
    """
    Calculates intersection roots with a cone.
    """
    troot = np.empty(2, dtype=np.float32)
    
    term1 = xs**2 * (-yb**2 + ay**2 * zb**2)
    term2 = ax**2 * (ys * zb - yb * zs)**2
    term3 = 2 * xb * xs * (yb * ys - ay**2 * zb * zs)
    term4 = xb**2 * (-ys**2 + ay**2 * zs**2)
    
    numerator = 4 * ax**2 * ay**2 * (term1 + term2 + term3 + term4)
    
    denom_term = (ax**2 * (yb - ys)**2 + 
                  ay**2 * (xb**2 - 2*xb*xs + xs**2 - ax**2 * (zb - zs)**2))
    
    # Avoid division by zero
    if denom_term == 0.0:
        troot[0] = 1.0
        troot[1] = -1.0
        return troot

    dist2 = numerator / (denom_term**2)
    
    avet_num = (ax**2 * ys * (ys - yb) + 
                ay**2 * (-(xb * xs) + xs**2 + ax**2 * (zb - zs) * zs))
    avet = avet_num / denom_term

    sd2 = math.sqrt(max(0.0, dist2)) * 0.5

    if dist2 < 0.0:
        troot[0] = 1.0
        troot[1] = -1.0
    else:
        tr1 = avet - sd2
        tr2 = avet + sd2
        
        zr1 = zs * (1.0 - tr1) + zb * tr1
        zr2 = zs * (1.0 - tr2) + zb * tr2

        # Logic for cone halves (z > 0 vs z < 0)
        if zr1 <= 0.0 and zr2 <= 0.0:
            troot[0] = 1.0
            troot[1] = -1.0
        elif zr1 > 0.0 and zr2 > 0.0:
            troot[0] = tr1
            troot[1] = tr2
        else:
            # One root positive, one negative
            if zr1 >= 0.0:
                if zs < zb:
                    troot[0] = tr1
                    troot[1] = ALARGE
                else:
                    troot[0] = -ALARGE
                    troot[1] = tr1
            if zr2 >= 0.0:
                if zs < zb:
                    troot[0] = tr2
                    troot[1] = ALARGE
                else:
                    troot[0] = -ALARGE
                    troot[1] = tr2
    return troot

# ==========================================
# Analytical Projectors
# ==========================================

@jit(nopython=True, cache=True)
def gen_shape_proj(sinomat, indsino, frame_vectors, 
                   nu, nv, du, dv, u0, v0, att,
                   nsurf, xu1, xu2, xu3, yu1, yu2, yu3, zu1, zu2, zu3, x0, y0, z0,
                   nplane, nellipsoid, ax_ell, ay_ell, az_ell):
    
    xSource, ySource, zSource = frame_vectors[0], frame_vectors[1], frame_vectors[2]
    xDetCenter, yDetCenter, zDetCenter = frame_vectors[6], frame_vectors[7], frame_vectors[8]
    eux, euy, euz = frame_vectors[9], frame_vectors[10], frame_vectors[11]
    evx, evy, evz = frame_vectors[12], frame_vectors[13], frame_vectors[14]

    for jp in range(nu):
        for kp in range(nv):
            u = u0 + ((jp + 1) - 0.5) * du 
            v = v0 + ((kp + 1) - 0.5) * dv 

            if indsino[jp, kp] == 1:
                xbin = xDetCenter + eux*u + evx*v
                ybin = yDetCenter + euy*u + evy*v
                zbin = zDetCenter + euz*u + evz*v
                
                sbdist = math.sqrt((xSource-xbin)**2 + (ySource-ybin)**2 + (zSource-zbin)**2)
                
                froot0 = -ALARGE
                froot1 = ALARGE
                
                isurf = 0 
                
                # Planes
                for iplane in range(1, nplane):
                    isurf += 1
                    xt = xSource - x0[isurf]
                    yt = ySource - y0[isurf]
                    zt = zSource - z0[isurf]
                    
                    zs = zu1[isurf]*xt + zu2[isurf]*yt + zu3[isurf]*zt
                    
                    xt = xbin - x0[isurf]
                    yt = ybin - y0[isurf]
                    zt = zbin - z0[isurf]
                    
                    zb = zu1[isurf]*xt + zu2[isurf]*yt + zu3[isurf]*zt
                    
                    troot = root_plane(zs, zb)
                    froot0 = max(froot0, troot[0])
                    froot1 = min(froot1, troot[1])
                    if froot0 >= froot1: break 

                if froot0 < froot1:
                    # Ellipsoids
                    for iell in range(1, nellipsoid): 
                        isurf += 1
                        xt = xSource - x0[isurf]
                        yt = ySource - y0[isurf]
                        zt = zSource - z0[isurf]

                        xs = xu1[isurf]*xt + xu2[isurf]*yt + xu3[isurf]*zt
                        ys = yu1[isurf]*xt + yu2[isurf]*yt + yu3[isurf]*zt
                        zs = zu1[isurf]*xt + zu2[isurf]*yt + zu3[isurf]*zt

                        xt = xbin - x0[isurf]
                        yt = ybin - y0[isurf]
                        zt = zbin - z0[isurf]

                        xb = xu1[isurf]*xt + xu2[isurf]*yt + xu3[isurf]*zt
                        yb = yu1[isurf]*xt + yu2[isurf]*yt + yu3[isurf]*zt
                        zb = zu1[isurf]*xt + zu2[isurf]*yt + zu3[isurf]*zt
                        
                        troot = root_ellipsoid(xs, ys, zs, xb, yb, zb, 
                                             ax_ell[iell], ay_ell[iell], az_ell[iell])
                        froot0 = max(froot0, troot[0])
                        froot1 = min(froot1, troot[1])
                        if froot0 >= froot1: break

                tlen = 0.0
                if froot0 < froot1:
                    tlen = froot1 - froot0
                
                sinomat[jp, kp] = sbdist * tlen * att

@jit(nopython=True, cache=True)
def ellipsoid_proj(sinomat, indsino, frame_vectors, 
                   nu, nv, du, dv, u0, v0, 
                   ax, ay, az, x0_pos, y0_pos, z0_pos, alpha, beta, att):
    
    xSource, ySource, zSource = frame_vectors[0], frame_vectors[1], frame_vectors[2]
    xDetCenter, yDetCenter, zDetCenter = frame_vectors[6], frame_vectors[7], frame_vectors[8]
    eux, euy, euz = frame_vectors[9], frame_vectors[10], frame_vectors[11]
    evx, evy, evz = frame_vectors[12], frame_vectors[13], frame_vectors[14]

    for jp in range(nu):
        for kp in range(nv):
            u = u0 + ((jp + 1) - 0.5) * du 
            v = v0 + ((kp + 1) - 0.5) * dv 

            if indsino[jp, kp] == 1:
                xbin = xDetCenter + eux*u + evx*v
                ybin = yDetCenter + euy*u + evy*v
                zbin = zDetCenter + euz*u + evz*v

                xspp = (-x0_pos + xSource) * math.cos(beta) - (y0_pos - ySource) * math.sin(beta)
                yspp = (-y0_pos + ySource) * math.cos(beta) + (x0_pos - xSource) * math.sin(beta)
                zspp = zSource - z0_pos
                
                xsp = xspp * math.cos(alpha) + zspp * math.sin(alpha)
                ysp = yspp
                zsp = -xspp * math.sin(alpha) + zspp * math.cos(alpha)

                xbpp = (-x0_pos + xbin) * math.cos(beta) - (y0_pos - ybin) * math.sin(beta)
                ybpp = (-y0_pos + ybin) * math.cos(beta) + (x0_pos - xbin) * math.sin(beta)
                zbpp = zbin - z0_pos
                
                xbp = xbpp * math.cos(alpha) + zbpp * math.sin(alpha)
                ybp = ybpp
                zbp = -xbpp * math.sin(alpha) + zbpp * math.cos(alpha)

                term_dist = ((xbp - xsp)**2 + (ybp - ysp)**2 + (zbp - zsp)**2)
                term_num_1 = -(az**2 * (xsp*ybp - xbp*ysp)**2)
                term_num_2 = ay**2 * (az**2 * (xbp - xsp)**2 + ax**2 * (zbp - zsp)**2 - (xsp*zbp - xbp*zsp)**2)
                term_num_3 = ax**2 * (az**2 * (ybp - ysp)**2 - (ysp*zbp - ybp*zsp)**2)
                
                numerator = 4 * ax**2 * ay**2 * az**2 * term_dist * (term_num_1 + term_num_2 + term_num_3)
                denominator = (ax**2 * az**2 * (ybp - ysp)**2 + 
                               ay**2 * (az**2 * (xbp - xsp)**2 + ax**2 * (zbp - zsp)**2))**2

                pathLength2 = 0.0
                if denominator != 0.0:
                    pathLength2 = numerator / denominator

                if pathLength2 < 0.0: pathLength2 = 0.0
                
                sinomat[jp, kp] = math.sqrt(pathLength2) * att

# ==========================================
# Linear Interpolation Routines
# ==========================================

@jit(nopython=True, cache=True)
def sinoproj_linear(sinomat, indsino, frame_vectors, 
                    ns, nu, nv, du, dv, u0, v0, 
                    smat, dx, dy, dz, x0, y0, z0, nx, ny, nz):
    """
    Forward projection (Ray Tracing) using Linear Interpolation.
    """
    for ip in range(ns):
        xSource = frame_vectors[ip, 0]
        ySource = frame_vectors[ip, 1]
        zSource = frame_vectors[ip, 2]
        xDetCenter = frame_vectors[ip, 6]
        yDetCenter = frame_vectors[ip, 7]
        zDetCenter = frame_vectors[ip, 8]
        eux, euy, euz = frame_vectors[ip, 9], frame_vectors[ip, 10], frame_vectors[ip, 11]
        evx, evy, evz = frame_vectors[ip, 12], frame_vectors[ip, 13], frame_vectors[ip, 14]

        for jp in range(nu):
            for kp in range(nv):
                u = u0 + ((jp + 1) - 0.5) * du 
                v = v0 + ((kp + 1) - 0.5) * dv 

                if indsino[ip, jp, kp] == 1:
                    xbin = xDetCenter + eux*u + evx*v
                    ybin = yDetCenter + euy*u + evy*v
                    zbin = zDetCenter + euz*u + evz*v
                    
                    xl, yl, zl = x0, y0, z0
                    xdiff, ydiff, zdiff = xbin - xSource, ybin - ySource, zbin - zSource
                    xad, yad, zad = abs(xdiff) / dx, abs(ydiff) / dy, abs(zdiff) / dz
                    
                    total = 0.0
                    
                    if xad > yad and xad > zad:
                        yox, zox = ydiff / xdiff, zdiff / xdiff
                        travVoxlen = dx * math.sqrt(1.0 + yox**2 + zox**2)
                        for ix in range(nx):
                            x = xl + dx * (ix + 0.5)
                            y = ySource + yox * (x - xSource)
                            z = zSource + zox * (x - xSource)
                            iy = int(math.floor((y - y0) / dy - 0.5))
                            iz = int(math.floor((z - z0) / dz - 0.5))
                            if (0 <= iy < ny - 1) and (0 <= iz < nz - 1):
                                wy2 = (y - (dy * (iy + 0.5) + y0)) / dy
                                wy1 = 1.0 - wy2
                                wz2 = (z - (dz * (iz + 0.5) + z0)) / dz
                                wz1 = 1.0 - wz2
                                vi1 = wy1 * smat[ix, iy, iz] + wy2 * smat[ix, iy + 1, iz]
                                vi2 = wy1 * smat[ix, iy, iz + 1] + wy2 * smat[ix, iy + 1, iz + 1]
                                total += travVoxlen * (wz1 * vi1 + wz2 * vi2)

                    elif yad > zad:
                        xoy, zoy = xdiff / ydiff, zdiff / ydiff
                        travVoxlen = dy * math.sqrt(1.0 + xoy**2 + zoy**2)
                        for iy in range(ny):
                            y = yl + dy * (iy + 0.5)
                            x = xSource + xoy * (y - ySource)
                            z = zSource + zoy * (y - ySource)
                            ix = int(math.floor((x - x0) / dx - 0.5))
                            iz = int(math.floor((z - z0) / dz - 0.5))
                            if (0 <= ix < nx - 1) and (0 <= iz < nz - 1):
                                wx2 = (x - (dx * (ix + 0.5) + x0)) / dx
                                wx1 = 1.0 - wx2
                                wz2 = (z - (dz * (iz + 0.5) + z0)) / dz
                                wz1 = 1.0 - wz2
                                vi1 = wx1 * smat[ix, iy, iz] + wx2 * smat[ix + 1, iy, iz]
                                vi2 = wx1 * smat[ix, iy, iz + 1] + wx2 * smat[ix + 1, iy, iz + 1]
                                total += travVoxlen * (wz1 * vi1 + wz2 * vi2)

                    else:
                        xoz, yoz = xdiff / zdiff, ydiff / zdiff
                        travVoxlen = dz * math.sqrt(1.0 + xoz**2 + yoz**2)
                        for iz in range(nz):
                            z = zl + dz * (iz + 0.5)
                            x = xSource + xoz * (z - zSource)
                            y = ySource + yoz * (z - zSource)
                            ix = int(math.floor((x - x0) / dx - 0.5))
                            iy = int(math.floor((y - y0) / dy - 0.5))
                            if (0 <= ix < nx - 1) and (0 <= iy < ny - 1):
                                wx2 = (x - (dx * (ix + 0.5) + x0)) / dx
                                wx1 = 1.0 - wx2
                                wy2 = (y - (dy * (iy + 0.5) + y0)) / dy
                                wy1 = 1.0 - wy2
                                vi1 = wx1 * smat[ix, iy, iz] + wx2 * smat[ix + 1, iy, iz]
                                vi2 = wx1 * smat[ix, iy + 1, iz] + wx2 * smat[ix + 1, iy + 1, iz]
                                total += travVoxlen * (wy1 * vi1 + wy2 * vi2)

                    sinomat[ip, jp, kp] = total

@jit(nopython=True, cache=True)
def sinobackproject_linear(sinomat, indsino, frame_vectors,
                           ns, nu, nv, du, dv, u0, v0,
                           smat, dx, dy, dz, x0, y0, z0, nx, ny, nz):
    """
    Backprojection using Linear Interpolation.
    """
    for ip in range(ns):
        xSource = frame_vectors[ip, 0]
        ySource = frame_vectors[ip, 1]
        zSource = frame_vectors[ip, 2]
        xDetCenter = frame_vectors[ip, 6]
        yDetCenter = frame_vectors[ip, 7]
        zDetCenter = frame_vectors[ip, 8]
        eux, euy, euz = frame_vectors[ip, 9], frame_vectors[ip, 10], frame_vectors[ip, 11]
        evx, evy, evz = frame_vectors[ip, 12], frame_vectors[ip, 13], frame_vectors[ip, 14]

        for jp in range(nu):
            for kp in range(nv):
                u = u0 + ((jp + 1) - 0.5) * du
                v = v0 + ((kp + 1) - 0.5) * dv
                
                if indsino[ip, jp, kp] == 1:
                    dataval = sinomat[ip, jp, kp]
                    xbin = xDetCenter + eux*u + evx*v
                    ybin = yDetCenter + euy*u + evy*v
                    zbin = zDetCenter + euz*u + evz*v
                    
                    xl, yl, zl = x0, y0, z0
                    xdiff, ydiff, zdiff = xbin - xSource, ybin - ySource, zbin - zSource
                    xad, yad, zad = abs(xdiff)/dx, abs(ydiff)/dy, abs(zdiff)/dz
                    
                    if xad > yad and xad > zad:
                        yox, zox = ydiff / xdiff, zdiff / xdiff
                        travVoxlen = dx * math.sqrt(1.0 + yox**2 + zox**2)
                        for ix in range(nx):
                            x = xl + dx * (ix + 0.5)
                            y = ySource + yox * (x - xSource)
                            z = zSource + zox * (x - xSource)
                            iy = int(math.floor((y - y0) / dy - 0.5))
                            iz = int(math.floor((z - z0) / dz - 0.5))
                            if (0 <= iy < ny - 1) and (0 <= iz < nz - 1):
                                wy2 = (y - (dy * (iy + 0.5) + y0)) / dy
                                wy1 = 1.0 - wy2
                                wz2 = (z - (dz * (iz + 0.5) + z0)) / dz
                                wz1 = 1.0 - wz2
                                
                                val = dataval * travVoxlen
                                smat[ix, iy, iz]     += val * wy1 * wz1
                                smat[ix, iy+1, iz]   += val * wy2 * wz1
                                smat[ix, iy, iz+1]   += val * wy1 * wz2
                                smat[ix, iy+1, iz+1] += val * wy2 * wz2
                    
                    elif yad > zad:
                        xoy, zoy = xdiff / ydiff, zdiff / ydiff
                        travVoxlen = dy * math.sqrt(1.0 + xoy**2 + zoy**2)
                        for iy in range(ny):
                            y = yl + dy * (iy + 0.5)
                            x = xSource + xoy * (y - ySource)
                            z = zSource + zoy * (y - ySource)
                            ix = int(math.floor((x - x0) / dx - 0.5))
                            iz = int(math.floor((z - z0) / dz - 0.5))
                            if (0 <= ix < nx - 1) and (0 <= iz < nz - 1):
                                wx2 = (x - (dx * (ix + 0.5) + x0)) / dx
                                wx1 = 1.0 - wx2
                                wz2 = (z - (dz * (iz + 0.5) + z0)) / dz
                                wz1 = 1.0 - wz2
                                val = dataval * travVoxlen
                                smat[ix, iy, iz]     += val * wx1 * wz1
                                smat[ix+1, iy, iz]   += val * wx2 * wz1
                                smat[ix, iy, iz+1]   += val * wx1 * wz2
                                smat[ix+1, iy, iz+1] += val * wx2 * wz2

                    else:
                        xoz, yoz = xdiff / zdiff, ydiff / zdiff
                        travVoxlen = dz * math.sqrt(1.0 + xoz**2 + yoz**2)
                        for iz in range(nz):
                            z = zl + dz * (iz + 0.5)
                            x = xSource + xoz * (z - zSource)
                            y = ySource + yoz * (z - zSource)
                            ix = int(math.floor((x - x0) / dx - 0.5))
                            iy = int(math.floor((y - y0) / dy - 0.5))
                            if (0 <= ix < nx - 1) and (0 <= iy < ny - 1):
                                wx2 = (x - (dx * (ix + 0.5) + x0)) / dx
                                wx1 = 1.0 - wx2
                                wy2 = (y - (dy * (iy + 0.5) + y0)) / dy
                                wy1 = 1.0 - wy2
                                val = dataval * travVoxlen
                                smat[ix, iy, iz]     += val * wx1 * wy1
                                smat[ix+1, iy, iz]   += val * wx2 * wy1
                                smat[ix, iy+1, iz]   += val * wx1 * wy2
                                smat[ix+1, iy+1, iz] += val * wx2 * wy2

# ==========================================
# Nearest Neighbor Routines
# ==========================================

@jit(nopython=True, cache=True)
def sinoproj_nearest(sinomat, indsino, frame_vectors, 
                     ns, nu, nv, du, dv, u0, v0, 
                     smat, dx, dy, dz, x0, y0, z0, nx, ny, nz):
    """
    Forward projection using Nearest Neighbor interpolation.
    """
    for ip in range(ns):
        xSource, ySource, zSource = frame_vectors[ip, 0], frame_vectors[ip, 1], frame_vectors[ip, 2]
        xDetCenter, yDetCenter, zDetCenter = frame_vectors[ip, 6], frame_vectors[ip, 7], frame_vectors[ip, 8]
        eux, euy, euz = frame_vectors[ip, 9], frame_vectors[ip, 10], frame_vectors[ip, 11]
        evx, evy, evz = frame_vectors[ip, 12], frame_vectors[ip, 13], frame_vectors[ip, 14]

        for jp in range(nu):
            for kp in range(nv):
                u = u0 + ((jp + 1) - 0.5) * du 
                v = v0 + ((kp + 1) - 0.5) * dv 

                if indsino[ip, jp, kp] == 1:
                    xbin = xDetCenter + eux*u + evx*v
                    ybin = yDetCenter + euy*u + evy*v
                    zbin = zDetCenter + euz*u + evz*v
                    
                    xl, yl, zl = x0, y0, z0
                    xdiff, ydiff, zdiff = xbin - xSource, ybin - ySource, zbin - zSource
                    xad, yad, zad = abs(xdiff)/dx, abs(ydiff)/dy, abs(zdiff)/dz
                    
                    total = 0.0
                    
                    if xad > yad and xad > zad:
                        yox, zox = ydiff / xdiff, zdiff / xdiff
                        travVoxlen = dx * math.sqrt(1.0 + yox**2 + zox**2)
                        for ix in range(nx):
                            x = xl + dx * (ix + 0.5)
                            y = ySource + yox * (x - xSource)
                            z = zSource + zox * (x - xSource)
                            iy = int(math.floor((y - y0) / dy))
                            iz = int(math.floor((z - z0) / dz))
                            if (0 <= iy < ny) and (0 <= iz < nz):
                                total += travVoxlen * smat[ix, iy, iz]

                    elif yad > zad:
                        xoy, zoy = xdiff / ydiff, zdiff / ydiff
                        travVoxlen = dy * math.sqrt(1.0 + xoy**2 + zoy**2)
                        for iy in range(ny):
                            y = yl + dy * (iy + 0.5)
                            x = xSource + xoy * (y - ySource)
                            z = zSource + zoy * (y - ySource)
                            ix = int(math.floor((x - x0) / dx))
                            iz = int(math.floor((z - z0) / dz))
                            if (0 <= ix < nx) and (0 <= iz < nz):
                                total += travVoxlen * smat[ix, iy, iz]

                    else:
                        xoz, yoz = xdiff / zdiff, ydiff / zdiff
                        travVoxlen = dz * math.sqrt(1.0 + xoz**2 + yoz**2)
                        for iz in range(nz):
                            z = zl + dz * (iz + 0.5)
                            x = xSource + xoz * (z - zSource)
                            y = ySource + yoz * (z - zSource)
                            ix = int(math.floor((x - x0) / dx))
                            iy = int(math.floor((y - y0) / dy))
                            if (0 <= ix < nx) and (0 <= iy < ny):
                                total += travVoxlen * smat[ix, iy, iz]

                    sinomat[ip, jp, kp] = total

@jit(nopython=True, cache=True)
def sinobackproject_nearest(sinomat, indsino, frame_vectors,
                            ns, nu, nv, du, dv, u0, v0,
                            smat, dx, dy, dz, x0, y0, z0, nx, ny, nz):
    """
    Backprojection using Nearest Neighbor Interpolation.
    """
    for ip in range(ns):
        xSource, ySource, zSource = frame_vectors[ip, 0], frame_vectors[ip, 1], frame_vectors[ip, 2]
        xDetCenter, yDetCenter, zDetCenter = frame_vectors[ip, 6], frame_vectors[ip, 7], frame_vectors[ip, 8]
        eux, euy, euz = frame_vectors[ip, 9], frame_vectors[ip, 10], frame_vectors[ip, 11]
        evx, evy, evz = frame_vectors[ip, 12], frame_vectors[ip, 13], frame_vectors[ip, 14]

        for jp in range(nu):
            for kp in range(nv):
                u = u0 + ((jp + 1) - 0.5) * du
                v = v0 + ((kp + 1) - 0.5) * dv
                
                if indsino[ip, jp, kp] == 1:
                    dataval = sinomat[ip, jp, kp]
                    xbin = xDetCenter + eux*u + evx*v
                    ybin = yDetCenter + euy*u + evy*v
                    zbin = zDetCenter + euz*u + evz*v
                    
                    xl, yl, zl = x0, y0, z0
                    xdiff, ydiff, zdiff = xbin - xSource, ybin - ySource, zbin - zSource
                    xad, yad, zad = abs(xdiff)/dx, abs(ydiff)/dy, abs(zdiff)/dz
                    
                    if xad > yad and xad > zad:
                        yox, zox = ydiff / xdiff, zdiff / xdiff
                        travVoxlen = dx * math.sqrt(1.0 + yox**2 + zox**2)
                        for ix in range(nx):
                            x = xl + dx * (ix + 0.5)
                            y = ySource + yox * (x - xSource)
                            z = zSource + zox * (x - xSource)
                            iy = int(math.floor((y - y0) / dy))
                            iz = int(math.floor((z - z0) / dz))
                            if (0 <= iy < ny) and (0 <= iz < nz):
                                smat[ix, iy, iz] += dataval * travVoxlen
                    
                    elif yad > zad:
                        xoy, zoy = xdiff / ydiff, zdiff / ydiff
                        travVoxlen = dy * math.sqrt(1.0 + xoy**2 + zoy**2)
                        for iy in range(ny):
                            y = yl + dy * (iy + 0.5)
                            x = xSource + xoy * (y - ySource)
                            z = zSource + zoy * (y - ySource)
                            ix = int(math.floor((x - x0) / dx))
                            iz = int(math.floor((z - z0) / dz))
                            if (0 <= ix < nx) and (0 <= iz < nz):
                                smat[ix, iy, iz] += dataval * travVoxlen

                    else:
                        xoz, yoz = xdiff / zdiff, ydiff / zdiff
                        travVoxlen = dz * math.sqrt(1.0 + xoz**2 + yoz**2)
                        for iz in range(nz):
                            z = zl + dz * (iz + 0.5)
                            x = xSource + xoz * (z - zSource)
                            y = ySource + yoz * (z - zSource)
                            ix = int(math.floor((x - x0) / dx))
                            iy = int(math.floor((y - y0) / dy))
                            if (0 <= ix < nx) and (0 <= iy < ny):
                                smat[ix, iy, iz] += dataval * travVoxlen

# ==========================================
# Matrix Projectors & Adjustment Tools
# ==========================================

@jit(nopython=True, parallel=True, cache=True)
def matproj_linear(sinomat, projmats, 
                   ns, nu, nv, du, dv, u0, v0, 
                   smat, dx, dy, dz, x0, y0, z0, nx, ny, nz):
    """
    Voxel-driven forward projection using projection matrices.
    """
    for ip in prange(ns):
        p11, p12, p13, p14 = projmats[ip, 0, 0], projmats[ip, 0, 1], projmats[ip, 0, 2], projmats[ip, 0, 3]
        p21, p22, p23, p24 = projmats[ip, 1, 0], projmats[ip, 1, 1], projmats[ip, 1, 2], projmats[ip, 1, 3]
        p31, p32, p33, p34 = projmats[ip, 2, 0], projmats[ip, 2, 1], projmats[ip, 2, 2], projmats[ip, 2, 3]

        for ipi in range(nx):
            for jpi in range(ny):
                for kpi in range(nz):
                    # Coordinate scaling preserved from original Fortran
                    x = 10.0 * (x0 + (ipi + 0.5) * dx)
                    y = 10.0 * (y0 + (jpi + 0.5) * dy)
                    z = 10.0 * (z0 + (kpi + 0.5) * dz)

                    up = p11*x + p12*y + p13*z + p14
                    vp = p21*x + p22*y + p23*z + p24
                    sp = p31*x + p32*y + p33*z + p34
                    
                    if sp != 0:
                        u = up / sp
                        v = vp / sp
                        
                        iu = math.floor(u)
                        iv = math.floor(v)
                        
                        w2u, w2v = u - iu, v - iv
                        w1u, w1v = 1.0 - w2u, 1.0 - w2v
                        
                        w11, w12 = w1u * w1v, w1u * w2v
                        w21, w22 = w2u * w1v, w2u * w2v
                        
                        idx_u, idx_v = int(iu) - 1, int(iv) - 1
                        val = smat[ipi, jpi, kpi] * dx * dy * dz

                        # Splatting with bounds check
                        if 0 <= idx_u < nu and 0 <= idx_v < nv:
                            sinomat[ip, idx_u, idx_v] += w11 * val
                        if 0 <= idx_u + 1 < nu and 0 <= idx_v < nv:
                            sinomat[ip, idx_u + 1, idx_v] += w21 * val
                        if 0 <= idx_u < nu and 0 <= idx_v + 1 < nv:
                            sinomat[ip, idx_u, idx_v + 1] += w12 * val
                        if 0 <= idx_u + 1 < nu and 0 <= idx_v + 1 < nv:
                            sinomat[ip, idx_u + 1, idx_v + 1] += w22 * val

@jit(nopython=True, cache=True)
def sumadj(sinomat, indsino, ns, nu, nv, du, dv, u0, v0, 
           nps, npu, npv, coeffs):
    """
    Adjusts sinogram based on polynomial coefficients.
    """
    for ip in range(ns):
        s = 2.0 * ((ip + 1) - 1.0) / (ns - 1.0) - 1.0
        
        for jp in range(nu):
            if indsino[ip, jp, 0] == 1:
                u = 2.0 * (u0 + ((jp + 1) - 0.5) * du) / (nu * du)
                
                for kp in range(nv):
                    v = 2.0 * (v0 + ((kp + 1) - 0.5) * dv) / (nv * dv)
                    
                    val_sum = 0.0
                    for ipp in range(nps):
                        t1 = s**ipp
                        for jpp in range(npu):
                            t2 = u**jpp
                            for kpp in range(npv):
                                t3 = v**kpp
                                val_sum += coeffs[ipp, jpp, kpp] * t1 * t2 * t3
                    
                    sinomat[ip, jp, kp] += val_sum

@jit(nopython=True, cache=True)
def backsumadj(sinomat, indsino, ns, nu, nv, ds, du, dv, u0, v0, 
               nps, npu, npv, coeffs):
    """
    Back-calculates coefficients from sinogram.
    """
    # Initialize coeffs to 0
    coeffs[:] = 0.0

    for ip in range(ns):
        s = 2.0 * ((ip + 1) - 1.0) / (ns - 1.0) - 1.0
        
        for jp in range(nu):
            if indsino[ip, jp, 0] == 1:
                u = 2.0 * (u0 + ((jp + 1) - 0.5) * du) / (nu * du)
                
                for kp in range(nv):
                    v = 2.0 * (v0 + ((kp + 1) - 0.5) * dv) / (nv * dv)
                    
                    val = sinomat[ip, jp, kp]
                    
                    for ipp in range(nps):
                        t1 = s**ipp
                        for jpp in range(npu):
                            t2 = u**jpp
                            for kpp in range(npv):
                                t3 = v**kpp
                                coeffs[ipp, jpp, kpp] += val * t1 * t2 * t3
