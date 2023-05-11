import os
import numpy as np
import paraview.simple as pvs


TEMP_CSV_LOC = "temp.csv"


def get_data_from_csv(csv_file, remove_csv=True):
    """Pulls in data from a csv file.
    
    Parameters
    ----------
    csv_file : str
        File to pull the data in from.

    remove_csv : bool, optional
        Whether to delete the file the data are pulled from. Defaults to True.

    Returns
    -------
    column_headers : ndarray
        Array of column headers.

    data : ndarray
        Data array. Columns correspond to column headers.
    """

    # Read into arrays
    with open(csv_file, 'r') as data_file:
        column_headers = [x.strip().replace('"', '') for x in data_file.readline().split(',')]
    cell_data = np.genfromtxt(csv_file, delimiter=',', skip_header=1)

    # Remove csv
    if remove_csv:
        os.remove(csv_file)

    return column_headers, cell_data


def _save(proxy, which_data):
    # Saves the data from the proxy

    pvs.SaveData(TEMP_CSV_LOC, proxy=proxy, Precision=12, FieldAssociation="{0} Data".format(which_data.title()))


def extract_all_data(data_file, which_data='point'):
    """Gets all data in the given data file.
    
    Parameters
    ----------
    data_file : str
        Data filename. Must be .vtk.

    which_data : str, optional
        'point' or 'cell'. Defaults to 'point'.

    Returns
    -------
    column_headers : ndarray
        Array of column headers.

    data : ndarray
        Data array. Columns correspond to column headers.
    """

    # Read into ParaView
    thing = pvs.LegacyVTKReader(registrationName="name", FileNames=data_file)

    _save(thing, which_data)

    return get_data_from_csv(TEMP_CSV_LOC)


def get_data_column_from_array(headers, data, col_des):
    """Returns the data vector from data with the header col_des.
    
    Parameters
    ----------
    headers : list of str
        Lsit of column headers.

    data : ndarray
        Data array.

    col_des : str
        Header of column to pull data from.

    Returns
    -------
    vector
        Data in desired column.

    Raises
    ------
    ValueError
        If col_des is not found in headers.
    """

    # Get column index
    try:
        ind = headers.index(col_des)
    except ValueError:
        raise ValueError("Desired column header not found. Headers available: {0}".format(headers))

    return data[:,ind].flatten()


def extract_plane_slice(data_file, normal_vector, plane_origin, filter=False, which_data='point'):
    """Extracts a slice of the data on the plane defined by the given origin and normal vector.
    
    Parameters
    ----------
    data_filename : str
        Data filename. Must be .vtk.

    normal_vector : ndarray
        Normal vector of plane on which to slice.

    plane_origin : ndarray
        Origin of plane.

    filter : bool, optional
        Whether to apply a cell-data-to-point-data filter before slicing. Defaults to False.

    which_data : str, optional
        'point' or 'cell'. Defaults to 'point'.

    Returns
    -------
    column_headers : ndarray
        Array of column headers.

    slice_data : ndarray
        Data on the slice. Columns correspond to column headers.
    """

    # Read into ParaView
    thing = pvs.LegacyVTKReader(registrationName='name', FileNames=data_file)

    # Filter
    if filter:
        thing = pvs.CellDatatoPointData(registrationName='filtered', Input=thing)

    # Create slice
    slice = pvs.Slice(registrationName='slice', Input=thing)
    slice.SliceType.Origin = plane_origin
    slice.SliceType.Normal = normal_vector

    _save(slice, which_data)

    return get_data_from_csv(TEMP_CSV_LOC)