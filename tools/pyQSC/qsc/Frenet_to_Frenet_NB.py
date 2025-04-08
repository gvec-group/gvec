"""
This module contains the routines to compute
a given flux surface shape, defined in by the Frenet TNB frame (tangential, normal, binormal):

 rtilde(theta,phi) = r0(phi) + X(theta,phi)*N(phi) +Y(theta,phi)*B(phi) + Z(theta,phi)*T(phi)  (1)

and transform to the N-B plane of the same frenet frame (here zeta is only used to distinguish from phi, its the same angular parameter), such that the flux surface is only represented by two variables

 rhat(thet,zeta) = r0(zeta) + Xhat(theta,zeta)*N(zeta) + Yhat(theta,zeta)*B(zeta)   (2)

We will compute Xhat,Yhat at interpolation points theta_j,zeta_k. We have to find phi in (1) along a theta_j such that it lies in the N(zeta_k)-B(zeta_k) plane

   (rtilde(theta_j,phi) - r0(zeta_k)) . T(zeta_k) =0     => root finding for phi_jk(theta_j,zeta_k)

to finally compute

   Xhat(theta_j,zeta_k) = (r(theta_j,phi_jk)-r0(zeta_k)) . N(zeta_k)
   Yhat(theta_j,zeta_k) = (r(theta_j,phi_jk)-r0(zeta_k)) . B(zeta_k)

"""

import numpy as np
from scipy.optimize import root_scalar

def Frenet_to_Frenet_NB_residual_func(phi0, phi_target, qsc):
    """
    This function takes a point defined by r(phi0) and computes the distance to the N-B-plane given by phi_target
    (r(phi0) - r0(phi_target)). T(phi_target)

    since the T is the normal vector of the N-B-plane
    Args:
        phi0 (float): position along r0 defining r(phi0)
        phi_target (float) :  target position along r0
    """
    total_x,total_y,total_z = Frenet_to_xyz_1_point(phi0, qsc)

    sinphi = np.sin(phi_target)
    cosphi = np.cos(phi_target)

    R0          = qsc.R0_func(phi_target)
    Z0          = qsc.Z0_func(phi_target)
    tangent_R   = qsc.tangent_R_spline(phi_target)
    tangent_phi = qsc.tangent_phi_spline(phi_target)
    tangent_z   = qsc.tangent_z_spline(phi_target)

    tangent_x = tangent_R * cosphi - tangent_phi * sinphi
    tangent_y = tangent_R * sinphi + tangent_phi * cosphi


    distance_to_NB =( (total_x - R0*cosphi)*tangent_x
                     +(total_y - R0*sinphi)*tangent_y
                     +(total_z - Z0       )*tangent_z  )

    return distance_to_NB

def get_XYhat(phi0, phi_target, qsc):
    """
    This function takes a point defined by r(phi0) that lies already in the N-B plane (distance residual =0),
    and computes the new variables Xhat,Yhat in the N-B-plane given by phi_target
    Xhat=(r(phi0) - r0(phi_target)). N(phi_target)
    Yhat=(r(phi0) - r0(phi_target)). B(phi_target)

    Args:
        phi0 (float): position along r0 defining r(phi0)
        phi_target (float) :  target position along r0
    """

    total_x,total_y,total_z = Frenet_to_xyz_1_point(phi0, qsc)

    sinphi = np.sin(phi_target)
    cosphi = np.cos(phi_target)

    R0          = qsc.R0_func(phi_target)
    Z0          = qsc.Z0_func(phi_target)
    normal_R   = qsc.normal_R_spline(phi_target)
    normal_phi = qsc.normal_phi_spline(phi_target)
    normal_z   = qsc.normal_z_spline(phi_target)

    normal_x = normal_R * cosphi - normal_phi * sinphi
    normal_y = normal_R * sinphi + normal_phi * cosphi


    Xhat  =(  (total_x - R0*cosphi)*normal_x
             +(total_y - R0*sinphi)*normal_y
             +(total_z - Z0       )*normal_z  )

    binormal_R   = qsc.binormal_R_spline(phi_target)
    binormal_phi = qsc.binormal_phi_spline(phi_target)
    binormal_z   = qsc.binormal_z_spline(phi_target)

    binormal_x = binormal_R * cosphi - binormal_phi * sinphi
    binormal_y = binormal_R * sinphi + binormal_phi * cosphi

    Yhat  =(  (total_x - R0*cosphi)*binormal_x
             +(total_y - R0*sinphi)*binormal_y
             +(total_z - Z0       )*binormal_z  )

    return Xhat,Yhat

def Frenet_to_xyz_1_point(phi0,qsc):
    """
    This function takes a point on the magnetic axis with a given
    toroidal angle phi0 and computes the cartesian x,y,z coordinates of a point on the surface,
    using previously defined functions qsc.X_spline(phi0),qsc.Y_spline(phi0),qsc.Z_spline(phi0)

    rtilde(phi) = r0(phi) + X(phi)*N(phi) +Y(phi)*B(phi) + Z(phi)*T(phi)

    Args:
        phi0: toroidal angle on the axis
    """
    X_at_phi0    = qsc.X_spline(phi0)
    Y_at_phi0    = qsc.Y_spline(phi0)

    sinphi0 = np.sin(phi0)
    cosphi0 = np.cos(phi0)
    R0_at_phi0   = qsc.R0_func(phi0)
    z0_at_phi0   = qsc.Z0_func(phi0)
    normal_R     = qsc.normal_R_spline(phi0)
    normal_phi   = qsc.normal_phi_spline(phi0)
    normal_z     = qsc.normal_z_spline(phi0)
    binormal_R   = qsc.binormal_R_spline(phi0)
    binormal_phi = qsc.binormal_phi_spline(phi0)
    binormal_z   = qsc.binormal_z_spline(phi0)

    normal_x   =   normal_R * cosphi0 -   normal_phi * sinphi0
    normal_y   =   normal_R * sinphi0 +   normal_phi * cosphi0
    binormal_x = binormal_R * cosphi0 - binormal_phi * sinphi0
    binormal_y = binormal_R * sinphi0 + binormal_phi * cosphi0

    total_x = R0_at_phi0 * cosphi0 + X_at_phi0 * normal_x + Y_at_phi0 * binormal_x
    total_y = R0_at_phi0 * sinphi0 + X_at_phi0 * normal_y + Y_at_phi0 * binormal_y

    total_z = z0_at_phi0           + X_at_phi0 * normal_z + Y_at_phi0 * binormal_z

    if qsc.order != 'r1':
        Z_at_phi0    = qsc.Z_spline(phi0)
        tangent_R   = qsc.tangent_R_spline(phi0)
        tangent_phi = qsc.tangent_phi_spline(phi0)
        tangent_z   = qsc.tangent_z_spline(phi0)

        tangent_x = tangent_R * cosphi0 - tangent_phi * sinphi0
        tangent_y = tangent_R * sinphi0 + tangent_phi * cosphi0

        total_x += Z_at_phi0 * tangent_x
        total_y += Z_at_phi0 * tangent_y
        total_z += Z_at_phi0 * tangent_z


    return total_x, total_y, total_z

def Frenet_to_Frenet_NB(self, r, ntheta=None, theta_in = np.linspace(0,2*np.pi,20,endpoint=False)):
    r"""
    For a given minor radius coordinate :math:`r`, compute the
    shape of the flux surface in standard cylindrical coordinates
    :math:`(R, \phi, Z)`.  This function returns :math:`R` and
    :math:`Z` as 2D arrays corresponding to dimensions ``(theta,
    phi)``, where ``theta`` is the Boozer poloidal angle and ``phi``
    is the standard toroidal angle.  Also returned is ``phi0(theta,
    phi)``, defined as follows: for given ``(theta, phi)``, move to
    the magnetic axis while holding the Boozer poloidal and toroidal
    angles fixed; the standard toroidal angle at that resulting point
    on the axis is ``phi0``.

    Args:
        r: near-axis radius r of the desired boundary surface
        ntheta: resolution in the poloidal angle theta

    Returns: 2 element tuple containing ``(Xhat, Yhat, phi0 )``. Each entry has shape ``(ntheta, nphi)``.

    """
    nphi_conversion = self.nphi
    if(ntheta):
        theta = np.linspace(0,2*np.pi,ntheta,endpoint=False)
    else:
        ntheta = theta_in.size
        theta  = theta_in
    phi_conversion = np.linspace(0,2*np.pi/self.nfp,nphi_conversion,endpoint=False)
    Xhat_2D = np.zeros((ntheta,nphi_conversion))
    Yhat_2D = np.zeros((ntheta,nphi_conversion))
    phi0_2D = np.zeros((ntheta,nphi_conversion))
    for j_theta in range(ntheta):
        costheta = np.cos(theta[j_theta])
        sintheta = np.sin(theta[j_theta])
        X_at_this_theta = r * (self.X1c_untwisted * costheta + self.X1s_untwisted * sintheta)
        Y_at_this_theta = r * (self.Y1c_untwisted * costheta + self.Y1s_untwisted * sintheta)
        Z_at_this_theta = 0 * X_at_this_theta
        if self.order != 'r1':
            # We need O(r^2) terms:
            cos2theta = np.cos(2 * theta[j_theta])
            sin2theta = np.sin(2 * theta[j_theta])
            X_at_this_theta += r * r * (self.X20_untwisted + self.X2c_untwisted * cos2theta + self.X2s_untwisted * sin2theta)
            Y_at_this_theta += r * r * (self.Y20_untwisted + self.Y2c_untwisted * cos2theta + self.Y2s_untwisted * sin2theta)
            Z_at_this_theta += r * r * (self.Z20_untwisted + self.Z2c_untwisted * cos2theta + self.Z2s_untwisted * sin2theta)
            if self.order == 'r3':
                # We need O(r^3) terms:
                costheta  = np.cos(theta[j_theta])
                sintheta  = np.sin(theta[j_theta])
                cos3theta = np.cos(3 * theta[j_theta])
                sin3theta = np.sin(3 * theta[j_theta])
                r3 = r * r * r
                X_at_this_theta += r3 * (self.X3c1_untwisted * costheta + self.X3s1_untwisted * sintheta
                                         + self.X3c3_untwisted * cos3theta + self.X3s3_untwisted * sin3theta)
                Y_at_this_theta += r3 * (self.Y3c1_untwisted * costheta + self.Y3s1_untwisted * sintheta
                                         + self.Y3c3_untwisted * cos3theta + self.Y3s3_untwisted * sin3theta)
                Z_at_this_theta += r3 * (self.Z3c1_untwisted * costheta + self.Z3s1_untwisted * sintheta
                                         + self.Z3c3_untwisted * cos3theta + self.Z3s3_untwisted * sin3theta)
        self.X_spline = self.convert_to_spline(X_at_this_theta)
        self.Y_spline = self.convert_to_spline(Y_at_this_theta)
        self.Z_spline = self.convert_to_spline(Z_at_this_theta)
        for j_phi in range(nphi_conversion):
            # Solve for the phi0
            phi_target = phi_conversion[j_phi]
            phi0_rootSolve_min = phi_target - 1.0 / self.nfp
            phi0_rootSolve_max = phi_target + 1.0 / self.nfp
            res = root_scalar(Frenet_to_Frenet_NB_residual_func, xtol=1e-15, rtol=1e-15, maxiter=1000,\
                              args=(phi_target, self), bracket=[phi0_rootSolve_min, phi0_rootSolve_max], x0=phi_target)
            phi0_solution = res.root
            Xhat_2D[j_theta,j_phi] , Yhat_2D[j_theta,j_phi] = get_XYhat(phi0_solution,phi_target,self)
            phi0_2D[j_theta,j_phi] = phi0_solution

    return Xhat_2D, Yhat_2D, phi0_2D

def to_xyz(self,points):
    """
    Function to convert a set of points in (r,theta,phi0) coordinates
    where r=sqrt(2*psi/B0) is the near-axis radius, theta is the
    Boozer poloidal angle and phi0 is the cylindrical angle phi
    on the axis to cartesian x,y,z coordinates

    Args:
        points: an array of floats with dimension Nx3 with N the
        number of points to evaluate with each points having
        the (r,theta,phi0) values to evaluate
    """
    xcart = []
    ycart = []
    zcart = []
    for point in points:
        r      = point[0]
        theta  = point[1]
        phi0   = point[2]
        costheta = np.cos(theta)
        sintheta = np.sin(theta)
        cos2theta = np.cos(2*theta)
        sin2theta = np.sin(2*theta)
        cos3theta = np.cos(3*theta)
        sin3theta = np.sin(3*theta)
        X_at_this_theta = r * (self.X1c_untwisted * costheta + self.X1s_untwisted * sintheta)
        Y_at_this_theta = r * (self.Y1c_untwisted * costheta + self.Y1s_untwisted * sintheta)
        Z_at_this_theta = 0 * X_at_this_theta
        if self.order != 'r1':
            # We need O(r^2) terms:
            X_at_this_theta += r * r * (self.X20_untwisted + self.X2c_untwisted * cos2theta + self.X2s_untwisted * sin2theta)
            Y_at_this_theta += r * r * (self.Y20_untwisted + self.Y2c_untwisted * cos2theta + self.Y2s_untwisted * sin2theta)
            Z_at_this_theta += r * r * (self.Z20_untwisted + self.Z2c_untwisted * cos2theta + self.Z2s_untwisted * sin2theta)
            if self.order == 'r3':
                # We need O(r^3) terms:
                r3 = r * r * r
                X_at_this_theta += r3 * (self.X3c1_untwisted * costheta + self.X3s1_untwisted * sintheta
                                         + self.X3c3_untwisted * cos3theta + self.X3s3_untwisted * sin3theta)
                Y_at_this_theta += r3 * (self.Y3c1_untwisted * costheta + self.Y3s1_untwisted * sintheta
                                         + self.Y3c3_untwisted * cos3theta + self.Y3s3_untwisted * sin3theta)
                Z_at_this_theta += r3 * (self.Z3c1_untwisted * costheta + self.Z3s1_untwisted * sintheta
                                         + self.Z3c3_untwisted * cos3theta + self.Z3s3_untwisted * sin3theta)
        self.X_spline = self.convert_to_spline(X_at_this_theta)
        self.Y_spline = self.convert_to_spline(Y_at_this_theta)
        self.Z_spline = self.convert_to_spline(Z_at_this_theta)
        final_x, final_y, final_z = Frenet_to_xyz_1_point(phi0, self)
        xcart.append(final_x)
        ycart.append(final_y)
        zcart.append(final_z)

    return xcart, ycart, zcart
