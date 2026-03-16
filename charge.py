import math
import numpy as np
import re
import copy
from openpyxl import Workbook
import atoms
import nucleus


"""Please forgive me for not being very proficient in python. 
You don't need to specifically understand the content of my code. 
I'll mark the functions of the code and its input and output for you"""

def chargeRead(pathofCHGCAR, NGF, nucl={'Fe': 8, 'Cr': 12, 'V': 13}):
    """
    Extract the charge density on the FFT grid from a CHGCAR file.

    The function reads the CHGCAR file, determines the total number of electrons
    from the atomic species and their counts, calculates the box volume, locates
    the FFT grid dimensions, and reads the charge density values. It then normalizes
    the charge density so that the integrated density matches the total number of
    electrons.

    Parameters:
        pathofCHGCAR (str): Absolute or relative path to the CHGCAR file.
        NGF (list of int): FFT grid dimensions [NGXF, NGYF, NGZF].
        nucl (dict): Dictionary mapping element symbols to their number of electrons.
                     Default: {'Fe': 8, 'Cr': 12, 'V': 13}.

    Returns:
        dict: A dictionary where keys are reduced grid indices (as returned by
              location_reduction()) and values are the normalized charge densities
              at those grid points. The normalization ensures that the integrated
              charge density over the volume equals the total number of electrons.
    """
    with open(file=pathofCHGCAR, mode='rt', encoding='UTF-8') as f:
        chgcarlines = f.readlines()

    # 1. Calculate total number of electrons from atomic species
    types_line = chgcarlines[5]
    types_num_line = chgcarlines[6]
    types_list = re.split(r'\s+', types_line.strip())
    types_num_list = re.split(r'\s+', types_num_line.strip())
    num_electrons = 0
    for i in np.arange(len(types_num_list)):
        num_electrons += nucl[types_list[i]] * int(types_num_list[i])

    # 2. Compute box volume
    vx_line = chgcarlines[2]
    vy_line = chgcarlines[3]
    vz_line = chgcarlines[4]
    vx = np.array(np.float64(re.split(r'\s+', vx_line.strip())))
    vy = np.array(np.float64(re.split(r'\s+', vy_line.strip())))
    vz = np.array(np.float64(re.split(r'\s+', vz_line.strip())))
    vol_box = np.linalg.norm(vx) * np.linalg.norm(vy) * np.linalg.norm(vz)

    # 3. Locate the line containing the FFT grid dimensions
    NG_F_lin = '  {}  {}  {}\n'.format(NGF[0], NGF[1], NGF[2])
    NG_F_lin_index = chgcarlines.index(NG_F_lin)

    # 4. Read all subsequent lines and split into numbers
    contents = []
    for lin in chgcarlines[NG_F_lin_index + 1:]:
        lin_content = re.split(r'\s+', lin.strip())
        contents += lin_content

    # 5. Extract charge density values (first NGF[0]*NGF[1]*NGF[2] numbers)
    charge_dot = np.array(contents[0:NGF[0] * NGF[1] * NGF[2]], dtype=np.float64)

    # 6. First pass: compute sum of raw charge densities for normalization
    location_density_dic = {}
    total_sum_density = 0
    for dot_index in np.arange(NGF[0] * NGF[1] * NGF[2]):
        location_density_dic[location_reduction(dot_index, NGF)] = charge_dot[dot_index]
        total_sum_density += charge_dot[dot_index]

    # 7. Calculate normalization factor so that integrated density equals total electrons
    reduced_coefficient = num_electrons / (total_sum_density * vol_box) * (NGF[0] * NGF[1] * NGF[2])

    # 8. Second pass: apply normalization to obtain corrected charge densities
    location_density_dic = {}
    for dot_index in np.arange(NGF[0] * NGF[1] * NGF[2]):
        location_density_dic[location_reduction(dot_index, NGF)] = charge_dot[dot_index] * reduced_coefficient
        # Note: The normalization factor appears to produce incorrect results.
        # The user comment indicates: "Why does this yield incorrect results?"

    return location_density_dic


def chg_redistribution_calc(location_density_dic_1, location_density_dic_2, vBox):
    """
    Calculate the charge redistribution between two systems.

    This function computes the integrated absolute difference in charge density
    between two configurations (e.g., perfect system vs. system with a vacancy),
    normalized by the volume and number of grid points.

    Parameters:
        location_density_dic_1 (dict): Charge density dictionary for the reference system
                                       (e.g., fully relaxed system without vacancy).
        location_density_dic_2 (dict): Charge density dictionary for the system with a vacancy
                                       (fully relaxed).
        vBox (list of np.ndarray): Simulation box vectors as returned by boxGET().
                                   Assumes orthorhombic box where box lengths are the
                                   diagonal elements vBox[i][i].

    Returns:
        float: The charge redistribution value Δn, calculated as the sum over all
               grid points of |ρ₂ - ρ₁| * (volume / number of grid points).
    """
    # Calculate volume assuming orthorhombic cell (box vectors aligned with axes)
    vol = vBox[0][0] * vBox[1][1] * vBox[2][2]
    list_keys = location_density_dic_1.keys()
    delta_n = 0
    for key in list_keys:
        # Sum absolute difference multiplied by volume element (vol / N)
        delta_n += np.abs(location_density_dic_2[key] - location_density_dic_1[key]) * vol / len(list_keys)
    return delta_n

"""The remnant three functions are just attachment of the first one. """
def location_reduction(lin_index, NGF: list[int]):
    nNGXF = (lin_index + 1) % NGF[0]
    lin_index_1 = (lin_index + 1) // NGF[0]
    nNGYF = lin_index_1 % NGF[1]
    lin_index_2 = lin_index_1 // NGF[1]
    nNGZF = lin_index_2 % NGF[2]
    return tuple([nNGXF, nNGYF, nNGZF])


def coord_to_location(coord, vBOX, NGF):
    location = [np.rint(coord[0] / vBOX[0][0] * NGF[0]),
        np.rint(coord[1] / vBOX[1][1] * NGF[1]),
        np.rint(coord[2] / vBOX[2][2] * NGF[2])]
    for i in np.arange(3):
        if location[i] >= NGF[i]:
            location[i] -= NGF[i]
        elif location[i] < 0:
            location[i] += NGF[i]
    return tuple(location)

def location_to_coord(location, vBox, NGF):
    coord = np.array([location[0] / NGF[0] * vBox[0][0], location[1] / NGF[1] * vBox[1][1], location[2] / NGF[2] * vBox[2][2]])
    return coord