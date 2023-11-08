
#execute this script within the pyQSC root folder
def convert_pyQSC_to_GVEC_ncfile(qsc_stel_command='Qsc.from_paper("r1 section 5.1")',ncout_file ="axisNB_r1_section_5p1",boundary_m_max=10,r_bound=0.1):
    '''
        export axis and boundary data to netcdf format, that can then be read by GVEC
    '''
    from qsc import Qsc
    import numpy as np
    import netCDF4 as nc
    import os

    ###########################
    #stel = Qsc(rc=[1, 0.09], zs=[0, -0.09], nfp=2, etabar=0.95, I2=0.9, order='r2', B2c=-0.7, p2=-600000.)
    #stel = Qsc(rc=[1, 0.], zs=[0, ], nfp=3, etabar=0.95, I2=0., order='r1')
    #stel = Qsc.from_paper("2022 QH nfp7")
    #ncout_file ="axisNB_2022_QH_nfp7"
    ############################
    
    print(('QSC CASE: %s\n - netCDF output file: "%s" \n - boundary_m_max=%d \n - r_bound=%f'%(qsc_stel_command,ncout_file,boundary_m_max,r_bound)))
    stel=eval(qsc_stel_command)




    ###export axis data to netcdf format.

    nfp=stel.nfp
    axis_n_max=2*(stel.nfourier+1)

    axis_nzeta = 2*axis_n_max+1
    axis_nzetaFull = axis_nzeta*nfp
    zeta = np.linspace(0,2*np.pi,axis_nzetaFull,endpoint=False)
    zeta = zeta+0.5*zeta[1]  # half grid!

    xyz = np.zeros((3,axis_nzetaFull))

    #exact evaluation or R,Z
    r = sum([stel.rc[i]*np.cos(i*stel.nfp*zeta)+stel.rs[i]*np.sin(i*stel.nfp*zeta)  for i in range(len(stel.rc))])
    z = sum([stel.zc[i]*np.cos(i*stel.nfp*zeta)+stel.zs[i]*np.sin(i*stel.nfp*zeta)  for i in range(len(stel.zs))])
    xyz[0,:]= r*np.cos(zeta)
    xyz[1,:]= r*np.sin(zeta)
    xyz[2,:]= z

    #exact evaluation or R',Z'
    rp = sum([i*stel.nfp*(-stel.rc[i]*np.sin(i*stel.nfp*zeta)+stel.rs[i]*np.cos(i*stel.nfp*zeta))  for i in range(len(stel.rc))])
    zp = sum([i*stel.nfp*(-stel.zc[i]*np.sin(i*stel.nfp*zeta)+stel.zs[i]*np.cos(i*stel.nfp*zeta))  for i in range(len(stel.zs))])

    xyzp = np.zeros((3,axis_nzetaFull))
    xyzp[0,:]=-r*np.sin(zeta)+rp*np.cos(zeta)
    xyzp[1,:]= r*np.cos(zeta)+rp*np.sin(zeta)
    xyzp[2,:]= zp

    #exact evaluation or R'',Z''
    rpp = sum([-(i*stel.nfp)**2*(stel.rc[i]*np.cos(i*stel.nfp*zeta)+stel.rs[i]*np.sin(i*stel.nfp*zeta))  for i in range(len(stel.rc))])
    zpp = sum([-(i*stel.nfp)**2*(stel.zc[i]*np.cos(i*stel.nfp*zeta)+stel.zs[i]*np.sin(i*stel.nfp*zeta))  for i in range(len(stel.zs))])

    xyzpp = np.zeros((3,axis_nzetaFull))
    xyzpp[0,:]=-r*np.cos(zeta)-2*rp*np.sin(zeta)+rpp*np.cos(zeta)
    xyzpp[1,:]=-r*np.sin(zeta)+2*rp*np.cos(zeta)+rpp*np.sin(zeta)
    xyzpp[2,:]= zpp

    # compute TNB Frenet frame (without sign changes!)
    Txyz = np.zeros((3,axis_nzetaFull))
    Nxyz = np.zeros((3,axis_nzetaFull))
    Bxyz = np.zeros((3,axis_nzetaFull))

    lp=np.linalg.norm(xyzp,axis=0)
    Txyz[0,:]=xyzp[0,:]/lp
    Txyz[1,:]=xyzp[1,:]/lp
    Txyz[2,:]=xyzp[2,:]/lp

    Bxyz=np.cross(xyzp,xyzpp,axis=0)
    absB=np.linalg.norm(Bxyz,axis=0)
    Bxyz[0,:]=Bxyz[0,:]/absB
    Bxyz[1,:]=Bxyz[1,:]/absB
    Bxyz[2,:]=Bxyz[2,:]/absB

    Nxyz=np.cross(Bxyz,Txyz,axis=0)


    os.system("rm -f "+ncout_file+".nc")
    ncfile = nc.Dataset(ncout_file+'.nc', 'w') 
    vec_dim = ncfile.createDimension('vec',3)
    axis_nzeta_dim = ncfile.createDimension('nzeta_axis',axis_nzeta)
    #axis_nzeta_dim = ncfile.createDimension('axis/nzeta',axis_nzeta)
    axis_nzetaFull_dim = ncfile.createDimension('nzetaFull_axis',axis_nzetaFull)
    version=300
    for ivar,ival in zip(["NFP","VERSION","axis/n_max","axis/nzeta"],
                        ["nfp","version","axis_n_max","axis_nzeta"]):
        exec(ival+'_var = ncfile.createVariable("'+ivar+'","i8")')
        exec(ival+'_var.assignValue('+ival+')')

    zeta_var = ncfile.createVariable("axis/zeta(:)","double",("nzeta_axis"))
    zeta_var[:] = zeta[0:axis_nzeta]
    for vecvar,vecval in zip(["axis/xyz","axis/Nxyz","axis/Bxyz"],["xyz","Nxyz","Bxyz"]):
        exec(vecval+'_var = ncfile.createVariable("'+vecvar+'(::)","f8",("vec","nzetaFull_axis"))')

        exec(vecval+'_var[:,:] = '+vecval)




    boundary_ntheta=2*(boundary_m_max+1)

    boundary_nzeta=(stel.nphi)   # must be the same resolution (for now)
    boundary_n_max=(boundary_nzeta-1)//2

    if(stel.lasym):
        boundary_lasym=1
    else:
        boundary_lasym=0


    for ivar,ival in zip(["boundary/ntheta","boundary/nzeta","boundary/m_max","boundary/n_max","boundary/lasym"],
                        ["boundary_ntheta","boundary_nzeta","boundary_m_max","boundary_n_max","boundary_lasym"]):
        exec(ival+'_var = ncfile.createVariable("'+ivar+'","i8")')
        exec(ival+'_var.assignValue('+ival+')')

    boundary_ntheta_dim = ncfile.createDimension('ntheta_boundary',boundary_ntheta)

    theta=np.linspace(0,2*np.pi,boundary_ntheta,endpoint=False)    # NO HALF GRID HERE!
    theta_var = ncfile.createVariable("boundary/theta(:)","double",("ntheta_boundary"))
    theta_var[:] = theta

    boundary_nzeta_dim = ncfile.createDimension('nzeta_boundary',boundary_nzeta)

    zeta=np.linspace(0,2*np.pi/nfp,boundary_nzeta,endpoint=False)  # NO HALF GRID HERE!
    zeta_var = ncfile.createVariable("boundary/zeta(:)","double",("nzeta_boundary"))
    zeta_var[:] = zeta

    from qsc.Frenet_to_Frenet_NB import Frenet_to_Frenet_NB

    boundX,boundY,phi2d = Frenet_to_Frenet_NB(stel,r_bound, boundary_ntheta)

    for vecvar,vecval in zip(["boundary/X","boundary/Y"],["boundX","boundY"]):
        exec(vecval+'_var = ncfile.createVariable("'+vecvar+'(::)","f8",("ntheta_boundary","nzeta_boundary"))')
        exec(vecval+'_var[:,:] = '+vecval)


    ncfile.title = "== File that containts axis and boundary information, used in GVEC with the hmap_axisNB module"
    hdr=  "======= HEADER OF THE NETCDF FILE VERSION 3.0 ==================================="
    hdr+= "\n    Note: This file was generated from pyQSC data."
    hdr+= "\n=== FILE DESCRIPTION:"
    hdr+= "\n  * axis, normal and binormal of the frame are given in cartesian coordinates along the curve parameter zeta [0,2pi]."
    hdr+= "\n  * The curve is allowed to have a field periodicity NFP, but the curve must be provided on a full turn."
    hdr+= "\n  * The adata is given in real space, sampled along equidistant zeta point positions:"
    hdr+= "\n      zeta(i)=(i+0.5)/nzeta * (2pi/NFP), i=0,...,nzeta-1"
    hdr+= "\n    always shifted by (2pi/NFP) for the next field period."
    hdr+= "\n    Thus the number of points along the axis for a full turn is NFP*nzeta"
    hdr+= "\n  * definition of the axis-following frame in cartesian coordinates ( boundary surface at rho=1):"
    hdr+= "\n"
    hdr+= "\n     {x,y,z}(rho,theta,zeta)=axis_{x,y,z}(zeta) + X(rho,theta,zeta)*N_{x,y,z}(zeta)+Y(rho,theta,zeta)*B_{x,y,z}(zeta)  "
    hdr+= "\n"
    hdr+= "\n=== DATA DESCRIPTION"
    hdr+= "\n- general data"
    hdr+= "\n  * NFP: number of field periods"
    hdr+= "\n  * VERSION: version number as integer: V3.0 => 300"
    hdr+= "\n- 'axis' data group:"
    hdr+= "\n  * 'axis/n_max'   : maximum mode number in zeta (in one field period)"
    hdr+= "\n  * 'axis/nzeta'   : number of points along the axis, in one field period (>=2*n_max+1)"
    hdr+= "\n  * 'axis/zeta(:)' : zeta positions, 1D array of size 'axis/nzeta', for one field period. zeta[i]=zeta[1] + (i-1)/nzeta*(2pi/nfp). starting value arbitrary"
    hdr+= "\n  * 'axis/xyz(::)' : cartesian positions along the axis for ONE FULL TURN, 2D array of size (3,NFP* nzeta ), sampled at zeta positions,"
    hdr+= "\n                     xyz[:,j+fp*nzeta]=axis(zeta[j]+fp*2pi/NFP), for j=0,..nzeta-1 and  fp=0,...,NFP-1"
    hdr+= "\n  * 'axis/Nxyz(::)': cartesian components of the normal vector of the axis frame, 2D array of size (3, NFP* nzeta), evaluated analogously to the axis"
    hdr+= "\n  * 'axis/Bxyz(::)': cartesian components of the bi-normal vector of the axis frame, 2D array of size (3, NFP*nzeta), evaluated analogously to the axis"
    hdr+= "\n- 'boundary' data group:"
    hdr+= "\n  * 'boundary/m_max'    : maximum mode number in theta "
    hdr+= "\n  * 'boundary/n_max'    : maximum mode number in zeta (in one field period)"
    hdr+= "\n  * 'boundary/lasym'    : asymmetry, logical. "
    hdr+= "\n                           if lasym=0, boundary surface position X,Y in the N-B plane of the axis frame can be represented only with"
    hdr+= "\n                             X(theta,zeta)=sum X_mn*cos(m*theta-n*NFP*zeta), with {m=0,n=0...n_max},{m=1...m_max,n=-n_max...n_max}"
    hdr+= "\n                             Y(theta,zeta)=sum Y_mn*sin(m*theta-n*NFP*zeta), with {m=0,n=1...n_max},{m=1...m_max,n=-n_max...n_max}"
    hdr+= "\n                           if lasym=1, full fourier series is taken for X,Y"
    hdr+= "\n  * 'boundary/ntheta'    : number of points in theta (>=2*m_max+1)"
    hdr+= "\n  * 'boundary/nzeta'     : number of points in zeta  (>=2*n_max+1), can be different to 'axis/nzeta' !"
    hdr+= "\n  * 'boundary/theta(:)'  : theta positions, 1D array of size 'boundary/ntheta',  theta[i]=theta[1] + (i-1)/ntheta*(2pi), starting value arbitrary"
    hdr+= "\n  * 'boundary/zeta(:)'   : zeta positions, 1D array of size 'boundary/nzeta', for one field period! zeta[i]=zeta[1] + (i-1)/nzeta*(2pi/nfp). starting value arbitrary"
    hdr+= "\n  * 'boundary/X(::)',"
    hdr+= "\n    'boundary/Y(::)'     : boundary position X,Y in the N-B plane of the axis frame, in one field period, 2D array of size(ntheta, nzeta),  with"
    hdr+= "\n                              X[i, j]=X(theta[i],zeta[j])"
    hdr+= "\n                              Y[i, j]=Y(theta[i],zeta[j]), i=0...ntheta-1,j=0...nzeta-1"

    ncfile.header = hdr
    ncfile.close()

    print(('NETCDF FILE "%s" WRITTEN! '%(ncout_file)))

if __name__ == "__main__":
    convert_pyQSC_to_GVEC_ncfile(r_bound=0.2)