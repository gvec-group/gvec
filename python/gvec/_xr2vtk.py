# OUTPUT to VTK (.vts file) for visualization in paraview
import xarray as xr
import numpy as np
from numpy.typing import ArrayLike
import warnings
from pathlib import Path

def extend_period(ndarr: np.ndarray,extendperiod: ArrayLike = None):
    """
    This function takes a multidimensional array and extends any dimension with 1 element depending on the entries of extendperiod:
    - 'off': No extension is performed on the dimension.
    - 'copy_first': Copies the first element of the dimension to extend it.
    - 'extrapolate': Uses linear extrapolation of the last two slices to extend the dimension.

    Args:
    ndarr (np.ndarray): Input array to be extended.
    extendperiod (ArrayLike): An iterable of strings specifying the extension method for each dimension.
                         Options are 'off', 'copy_first', and 'extrapolate'. The iterable length
                         must match the number of dimensions in ndarr.

    Returns:
    np.ndarray: An array with extended dimensions according to the specified methods.

    Raises:
    AssertionError: If extendperiod length does not match ndarr dimensions, or if an invalid 
                    method is specified, or if a dimension with 'copy_first' or 'extrapolate' 
                    method has a length less than or equal to 2.
    """

    ndims=ndarr.ndim
    allowed_methods=['off','copy_first','extrapolate']
    if(len(extendperiod)==0) or (extendperiod is None):
        return ndarr
    else:
        assert len(extendperiod)==ndims , f"extentperiod must be an 1d-array or iterable of strings of length={ndims}, with one entry per dimension of the input ndarray"
        for idim,extper in enumerate(extendperiod):
            assert any(extper == s for s in allowed_methods), f"extendperiod only allows: {allowed_methods}, but was set to '{extper}'"
            if(extper !="off"):
                assert (ndarr.shape[idim] > 2) , f"cannot extendperiod in dimension {idim}, as it must have a length >2"
        

    ext_dims=np.asarray([(1 if s!="off" else 0) for s in extendperiod])

    new_shape=ndarr.shape+ext_dims
    ndarr_out=np.zeros(new_shape)+np.nan

    start=np.asarray([0]*ndims)
    end=np.asarray(ndarr.shape)
    incr=np.asarray([1]*ndims)
    ndslice=tuple([slice(a,b,c) for a,b,c in zip(start,end,incr)])
    ndarr_out[ndslice]=ndarr

    for idim,extper in enumerate(extendperiod):
        if(extper =="off"):
            continue
        start[idim]=-1; incr[idim]=-1;end[idim]=start[idim]-1
        ndslice_last=tuple([slice(a,b,c) for a,b,c in zip(start,end,incr)])
        if(extper=="copy_first"):
            start[idim]=0; incr[idim]=1;end[idim]=1
            ndslice_first=tuple([slice(a,b,c) for a,b,c in zip(start,end,incr)])
            ndarr_out[ndslice_last]=ndarr_out[ndslice_first]
        if(extper=="extrapolate"):
            start[idim]=-2; incr[idim]=-1;end[idim]=start[idim]-1
            ndslice_seclast=tuple([slice(a,b,c) for a,b,c in zip(start,end,incr)])
            start[idim]=-3; incr[idim]=-1;end[idim]=start[idim]-1
            ndslice_trdlast=tuple([slice(a,b,c) for a,b,c in zip(start,end,incr)])
            ndarr_out[ndslice_last]=2*ndarr_out[ndslice_seclast]-ndarr_out[ndslice_trdlast]
        # set dimension to new size
        start[idim]=0
        incr[idim]=1
        end[idim]=ndarr_out.shape[idim] 
    try:
        assert (not np.isnan(np.sum(ndarr_out)))
    except AssertionError:
        warnings.warn("NaN values in extend_period", UserWarning)
    return ndarr_out



def write_to_vtk(
    filename: Path,
    xyz_pos: np.ndarray,
    dim_names = ArrayLike | str,
    scals: dict = {},
    vecs: dict = {},
    extend_pol: bool = False,
    extend_tor: bool = False,
    coord_arrays: dict = {},
    coord_varnames = {"rad":"rho","pol":"theta","tor":"zeta"}
    ):

    """
    Writes grid and vector data to a VTK file for visualization.

    This function takes grid positions, scalar fields, and vector fields, and
    writes them into a VTK file format suitable for visualization in tools like
    ParaView. It supports extending the grid for periodic boundary conditions
    in specified dimensions.

    Args:
        filename (Path): The name of the output VTK file, without the extension.
        xyz_pos (np.ndarray): N-dimensional array representing the grid positions.
        dim_names (ArrayLike or str, optional): Iterable of strings representing the names of each dimension.
        scals (dict, optional): Dictionary of scalar fields to include, keyed by field name.
        vecs (dict, optional): Dictionary of vector fields to include, keyed by field name.
        extend_pol (bool, optional): Whether to extend the grid in the poloidal direction.
        extend_tor (bool, optional): Whether to extend the grid in the toroidal direction.
        coord_arrays (dict, optional): Dictionary of coordinate arrays, keyed by dimension name.
        coord_varnames (dict, optional): Dictionary mapping coordinate directions to variable names.

    Returns:
        str: A message indicating the VTK file has been written, or an error message if the VTK library is not available.

    Raises:
        AssertionError: If input arrays do not match expected dimensions or conditions.

    Notes:
        - The function requires the `pyevtk` library to be installed.
        - Supports up to 3-dimensional grids.
        - Only variables with dimensions that are subsets of the grid dimensions are processed.
    """

    try:
        from pyevtk.hl import gridToVTK
    except ImportError as e:
        warnings.warn(f"failed to import pyevtk.hl.gridToVTK: {e}",UserWarning)
        return "no vtk file written"
    
    ndims_vec  = xyz_pos.ndim   # number of dimensions of vector fields
    ndims_grid = ndims_vec-1 # number of dimensions of scalar fields

    assert len(dim_names)==ndims_vec , f"dim_names must be an 1d-array or strings of length={ndims_vec}, with one entry per dimension of the input ndarray xyz_pos"
    
    allowed_dim_names=["xyz","rad","pol","tor"]
    dims={"xyz":{"idim":None,"pos":None},
          "rad":{"idim":None,"pos":None},
          "pol":{"idim":None,"pos":None},
          "tor":{"idim":None,"pos":None}}
    
    for idim,dim_name in enumerate(dim_names):
        assert any(dim_name == s for s in allowed_dim_names), f"dim_names only allows: '{allowed_dim_names}', but was set to '{dim_name}'"
        dims[dim_name]["idim"]=idim
    
    for dim_name,pos in coord_arrays.items():
        assert any(dim_name == s for s in allowed_dim_names), f"key in coord_arrays only allows: '{allowed_dim_names}', but was set to '{dim_name}'"
        assert dims[dim_name]["idim"] is not None , "key in coord_arrays must appear also in dim_names"
        dims[dim_name]["pos"]=pos 
        assert any(dim_name == s for s in coord_varnames.keys()) , "key in coord_arrays must also be a key in coord_varnames"
    
    
    start=np.asarray([0]*ndims_vec)
    end  =np.asarray([None]*ndims_vec) # until the end

    idim=dims["xyz"]["idim"]

    start[idim]=0
    end[idim]=1
    ndslice_x =tuple([slice(a,b) for a,b in zip(start,end)])

    start[idim]=1
    end[idim]=2
    ndslice_y =tuple([slice(a,b) for a,b in zip(start,end)])

    start[idim]=2
    end[idim]=3
    ndslice_z =tuple([slice(a,b) for a,b in zip(start,end)])
    # get the shape of the mesh, if not 3d, arrays need to be reshaped with one or two extra dimensions
    
    shape_grid_in=np.asarray([1]*ndims_grid)
    shape_grid_out=np.asarray([1]*ndims_grid)

    copy_first_vec=["off"]*ndims_vec
    copy_first_grid=["off"]*ndims_grid

    extrapolate_vec=["off"]*ndims_vec
    extrapolate_grid=["off"]*ndims_grid

    to3d_grid=np.asarray([1]*3)

    i=0
    for idim,dim_name in enumerate(dim_names):
        if dim_name == "xyz":
            continue
        dims[dim_name]["i"]=i
        shape_grid_in[i]=xyz_pos.shape[idim]
        shape_grid_out[i]=xyz_pos.shape[idim]
        if (dim_name == "pol" and extend_pol) or (dim_name == "tor" and extend_tor): 
            shape_grid_out[i]+=1
            copy_first_vec[idim]="copy_first"
            copy_first_grid[i]="copy_first"
            extrapolate_vec[idim]="extrapolate"
            extrapolate_grid[i]="extrapolate"
        i+=1
    
    to3d_grid[0:ndims_grid]=shape_grid_out

    coords=extend_period(xyz_pos,extendperiod=copy_first_vec)
    xcoord=coords[ndslice_x].reshape(to3d_grid)
    ycoord=coords[ndslice_y].reshape(to3d_grid)
    zcoord=coords[ndslice_z].reshape(to3d_grid)
    
    ptdata_out={}
    scal_in=np.zeros(shape_grid_in)
    # add coord_arrays to output, and extend period in poloidal and torodial direction by extrapolation
    for dim_name,pos in coord_arrays.items():
        if((pos.shape == shape_grid_in).all()):
            ptdata_out[coord_varnames[dim_name]]=extend_period(
                pos,extendperiod=extrapolate_grid).reshape(to3d_grid)
        elif(pos.ndim == 1):
            i=dims[dim_name]['i']
            assert (len(pos)==scal_in.shape[i]), f"1d position array of dimension {dim_name} has size {len(pos)} but should match the grid size (={scal_in.shape[i]})"
            if(i==0):
                pos_out=np.expand_dims(pos,axis=-1)#behind
                if(ndims_grid==3):
                    pos_out=np.expand_dims(pos_out,axis=-1)#bedind
            elif(i==1):
                pos_out=np.expand_dims(pos,axis=0)#before
                if(ndims_grid==3):
                    pos_out=np.expand_dims(pos_out,axis=-1)#bedind
            elif(i==2):
                pos_out=np.expand_dims(np.expand_dims(pos,axis=0),axis=0) #can only be 3d
            ptdata_out[coord_varnames[dim_name]]=extend_period(
                pos_out+scal_in,extendperiod=extrapolate_grid).reshape(to3d_grid)
    
            
    for key,val in scals.items():
        assert (val.shape == shape_grid_in).all() , f"scalar field {key} must be of the same shape as the one of the xyz_pos components"
        ptdata_out[key]=extend_period(val,extendperiod=copy_first_grid).reshape(to3d_grid)
    for key,val in vecs.items():
        assert (val.shape == xyz_pos.shape) , f"vector field {key} must be of the same shape as  xyz_pos"
        vecval=extend_period(val,extendperiod=copy_first_vec)
        ptdata_out[key]=(vecval[ndslice_x].reshape(to3d_grid),
                         vecval[ndslice_y].reshape(to3d_grid),
                         vecval[ndslice_z].reshape(to3d_grid))
    fn=gridToVTK(
         filename,
         xcoord,ycoord,zcoord,
         pointData=ptdata_out,
         )
    return (fn + " written.")


def xr2vtk(
    filename: Path,
    xrds:xr.Dataset,
    dim_vec: str = "xyz",
    grid_pos: str = "pos",
    extend_dims: ArrayLike = None,
    coord_names: dict = {"rad":"rho","pol":"theta","tor":"zeta"}
    ):   
    """
    Converts an xarray.Dataset to a VTK file for visualization.

    This function takes an xarray.Dataset, extracts grid and vector information, and writes it to a VTK file in a structured grid format.

    Args:
        filename (Path): The name of the output VTK file, without the file extension.
        xrds (xarray.Dataset): The input dataset containing grid and variable information.
        dim_vec (str, optional): The name of the vector dimension in the dataset. Defaults to "xyz".
        grid_pos (str, optional): The variable name in the dataset that contains grid positions. Defaults to "pos".
        extend_dims (ArrayLike, optional): Iterable of dimensions to extend by one point for periodicity. Defaults to None.
        coord_names (dict, optional): A dictionary mapping coordinate directions to their variable names in the dataset. Defaults to {"rad": "rho", "pol": "theta", "tor": "zeta"}.

    Returns:
        str: A message indicating the VTK file has been written, or an error message if the VTK library is not available.

    Raises:
        AssertionError: If the input dataset or its variables do not meet the required conditions for VTK output.
        UserWarning: If the VTK library is not available, a warning is issued and no file is written.

    Notes:
        The function requires the `pyevtk` library to be installed for writing the VTK file.
        Supports up to 3-dimensional grids. Only variables with dimensions that are subsets of the grid dimensions are processed.
    """

    try:
        from pyevtk.hl import gridToVTK
    except ImportError as e:
        warnings.warn(f"failed to import pyevtk.hl.gridToVTK: {e}",UserWarning)
        return "no vtk file written"
    
    assert (grid_pos in xrds) , f"grid_pos variable ({grid_pos}) not found in dataset"
    dim_names=xrds[grid_pos].dims
    assert (dim_vec in dim_names) , f"{grid_pos} variable does not have the vector dimension {dim_vec}"
    dims_grid=[s for s in dim_names if s != dim_vec]

    ndims_grid=len(dims_grid)
    assert ndims_grid <=3 , f"Expected 3 dimensional grids for vtk output but found {ndims_grid}."

    if extend_dims is not None:
        for s in extend_dims:
            assert (s in dims_grid) ,f"extend_dims has a dimension name {s} that is not part of the grid_pos dimensions {dims_grid}"
    else:
        extend_dims=[]


    #take only variables with dimension names being a subset of dims_grid, ignoring dim_vec
    outvars=[var for var in xrds.data_vars if np.all([dim in dims_grid for dim in [s for s in xrds[var].dims if s != "xyz"] ])]
    scalar_vars=[var for var in outvars if (dim_vec not in xrds[var].dims)]
    vector_vars=[var for var in outvars if (dim_vec     in xrds[var].dims)]
        

    dims_vec=np.concatenate(([dim_vec],dims_grid))

    xcoord,ycoord,zcoord = xrds[grid_pos].transpose(*dims_vec).values
    
    shape_scal=np.asarray([1]*ndims_grid)
    extrapolate_grid=["off"]*ndims_grid
    copy_first_grid=["off"]*ndims_grid

    for idim,dim in enumerate(dims_grid):   
        shape_scal[idim]=xrds.sizes[dim]
        if dim in extend_dims:
            extrapolate_grid[idim]="extrapolate"
            copy_first_grid[idim]="copy_first"
    

    to3d_grid=np.asarray([1]*3)
    to3d_grid[0:ndims_grid]=shape_scal

    ptdata_out={}
    #coordinates need to be broadcasted, and extrapolated
    bc_scal=0*xrds.coords[coord_names[dims_grid[0]]]
    for dim in dims_grid[1:]:
        bc_scal=bc_scal+0*xrds.coords[coord_names[dim]]
    for dim in dims_grid:
        cname=coord_names[dim]
        coord=(xrds.coords[cname]+bc_scal).transpose(*dims_grid).copy()
        ptdata_out[cname]=np.ascontiguousarray(extend_period(coord,extendperiod=extrapolate_grid))
  

    ds_scal=xr.broadcast(xrds[scalar_vars])[0].transpose(*dims_grid)
    for var in scalar_vars:
        scal=ds_scal[var].values.copy()
        ptdata_out[var]=np.ascontiguousarray(extend_period(scal,extendperiod=copy_first_grid))
    
    ds_vec=xr.broadcast(xrds[vector_vars])[0].transpose(*dims_vec) # put vector to first dimension
    for var in vector_vars:
        if(var == grid_pos): 
            continue
        vx,vy,vz=ds_vec[var].values.copy()
        ptdata_out[var]=(np.ascontiguousarray(extend_period(vx,extendperiod=copy_first_grid)),
                         np.ascontiguousarray(extend_period(vy,extendperiod=copy_first_grid)),
                         np.ascontiguousarray(extend_period(vz,extendperiod=copy_first_grid)))
    fn=gridToVTK(
         filename,
         np.ascontiguousarray(extend_period(xcoord,extendperiod=copy_first_grid)),
         np.ascontiguousarray(extend_period(ycoord,extendperiod=copy_first_grid)),
         np.ascontiguousarray(extend_period(zcoord,extendperiod=copy_first_grid)),
         pointData=ptdata_out,
         )
    return (fn + " written.")
