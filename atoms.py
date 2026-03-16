import copy
import time

import numpy as np
import re

"""Please forgive me for not being very proficient in python. 
You don't need to specifically understand the content of my code. 
I'll mark the functions of the code and its input and output for you"""

class atom:
    """
    A class representing an atom with its type and Cartesian coordinates.

    Attributes:
        typ: Atomic type/species (can be string, integer, etc.).
        coord (np.ndarray): 3D coordinates as a numpy array of float64.
    """

    def __init__(self, typ, x: np.float64, y: np.float64, z: np.float64):
        """
        Initialize an atom with its type and coordinates.

        Parameters:
            typ: Atomic type/species.
            x (np.float64): x-coordinate.
            y (np.float64): y-coordinate.
            z (np.float64): z-coordinate.
        """
        self.typ = typ
        self.coord = np.array([x, y, z], dtype=np.float64)


def free_energy_obtain(pathofOUTCAR: str):
    """
    Extract the final free energy from an OUTCAR file.

    Parameters:
        pathofOUTCAR (str): Absolute or relative path to the OUTCAR file.

    Returns:
        float: The final free energy value (in eV) extracted from the line
               starting with 'free'. Returns 0.00 if the energy is not found.
    """
    with open(file=pathofOUTCAR, mode='rt', encoding='UTF-8') as f:
        lines_list_outcar = f.readlines()

    energy_final = 0.00
    # Search backwards from the end of the file to find the last occurrence of the 'free' line
    for lin in lines_list_outcar[::-1]:
        lin_list = re.split(r'\s+', lin.strip())
        try:
            if lin_list[0] == 'free':
                energy_final = np.float64(lin_list[4])
                # Debug print (commented out)
                # print(energy_final)
                break
        except:
            # Ignore lines that do not contain the expected format
            pass

    return energy_final

def boxGET(pathofPOSCAR: str):
    """
    Extract the simulation box vectors from a POSCAR or CONTCAR file.

    Parameters:
        pathofPOSCAR (str): Absolute or relative path to the POSCAR/CONTCAR file.

    Returns:
        list: A list of three numpy arrays, each representing a box vector (v1, v2, v3).
              The vectors are converted to float64.
    """
    with open(file=pathofPOSCAR, mode='rt', encoding='UTF-8') as f:
        listoflinesR = f.readlines()

    # Remove empty lines and strip trailing 'T' characters if present (selective dynamics)
    listoflines = []
    for lin in listoflinesR:
        if lin.strip() != '':
            # Remove selective dynamics flags if present (e.g., 'T   T   T')
            listoflines.append(lin.strip().strip('   T   T   T'))

    # The box vectors are typically on lines 3, 4, 5 (0-based index: 2, 3, 4)
    vx = re.split(r'\s+', listoflines[2])
    vy = re.split(r'\s+', listoflines[3])
    vz = re.split(r'\s+', listoflines[4])

    # Convert each component to float64 and store as numpy arrays
    vBox = [np.array(np.float64(vx)), np.array(np.float64(vy)), np.array(np.float64(vz))]
    return vBox

def atomsGET(pathofPOSCAR, vBox):
    """
    Extract all atom positions from a POSCAR or CONTCAR file and create atom objects.

    The function reads the POSCAR/CONTCAR file, parses the element types and counts,
    converts coordinates from direct (fractional) to Cartesian if necessary,
    and returns a list of atom instances.

    Parameters:
        pathofPOSCAR (str): Absolute or relative path to the POSCAR/CONTCAR file.
        vBox (list of np.ndarray): Simulation box vectors as returned by boxGET().
                                   Each vector is a 3-element float64 array.

    Returns:
        list: A list of atom objects, each containing type and Cartesian coordinates.
              The order corresponds to the order in the POSCAR file.
    """
    with open(file=pathofPOSCAR, mode='rt', encoding='UTF-8') as f:
        listoflines = f.readlines()

    # Strip newline and possible selective dynamics flags from each line
    lines = []
    for i in np.arange(len(listoflines)):
        lines.append(listoflines[i].strip().strip('   T   T   T'))

    # Line 5: element names (0-indexed), line 6: number of each element
    elementline = lines[5]
    numElementsline = lines[6]

    elements = re.split(r'\s+', elementline)
    numElementslist = np.array(re.split(r'\s+', numElementsline), dtype=np.int32)

    # Expand element list: repeat each element according to its count
    elementlist = []
    for i in np.arange(len(elements)):
        for j in np.arange(numElementslist[i]):
            elementlist.append(elements[i])

    # Remove the 'Selective dynamics' line if present (it may appear after the lattice vectors)
    try:
        lines.remove(r'Selective dynamics')
    except:
        # No selective dynamics tag found, continue
        pass

    # Extract atomic positions (fractional or Cartesian depending on line 7)
    positionlist = []
    positions = []
    # Atomic coordinates start at line 8 (0-indexed) and go for total number of atoms
    for lin in lines[8:8 + sum(numElementslist)]:
        positionlist.append(np.array(re.split(r'\s+', lin), dtype=np.float64))

    # Check coordinate type: 'Direct' indicates fractional coordinates
    if lines[7] == 'Direct':
        # Convert fractional to Cartesian: pos = frac[0]*v1 + frac[1]*v2 + frac[2]*v3
        for position in positionlist:
            positions.append(position[0] * vBox[0] + position[1] * vBox[1] + position[2] * vBox[2])
    else:
        # Already Cartesian coordinates
        positions = positionlist

    # Create atom objects for each atom
    atomslist = []
    for i in np.arange(sum(numElementslist)):
        atomslist.append(atom(typ=elementlist[i], x=positions[i][0], y=positions[i][1], z=positions[i][2]))

    return atomslist


def poscarWrite(atomslist: list[atom], vBox, poscarname):
    """
    Write a POSCAR file based on the given atom list and simulation box vectors.

    The function generates a standard POSCAR format file with Cartesian coordinates.
    It automatically determines the unique element types and their counts from the atom list.

    Parameters:
        atomslist (list[atom]): List of atom objects, each containing type and Cartesian coordinates.
        vBox (list of np.ndarray): Simulation box vectors as a list of three 3-element float64 arrays.
        poscarname (str): Absolute or relative path for the output POSCAR/CONTCAR file.

    Returns:
        None. The function writes the POSCAR content to the specified file.
    """
    # Extract element types from each atom
    elementslist = []
    for thisatom in atomslist:
        elementslist.append(thisatom.typ)

    # Get unique element types in the order they first appear
    elements = []
    for element in elementslist:
        if element not in elements:
            elements.append(element)

    # Count number of atoms for each unique element type
    numElementslist = []
    for element in elements:
        numElementslist.append(elementslist.count(element))

    # Build POSCAR header lines
    firstlin = 'POSCAR\n'
    secondlin = '   1.00000000000000\n'  # Universal scaling factor

    # Format box vectors with 3x3 components
    vBoxlines = ('     {}   {}   {}\n'.format(vBox[0][0], vBox[0][1], vBox[0][2]) +
                 '     {}   {}   {}\n'.format(vBox[1][0], vBox[1][1], vBox[1][2]) +
                 '     {}   {}   {}\n'.format(vBox[2][0], vBox[2][1], vBox[2][2]))

    # Element names line
    elementlin = ''
    for element in elements:
        elementlin += '  ' + element + '  '
    elementlin += '\n'

    # Number of atoms per element line
    numElementslin = ''
    for numElement in numElementslist:
        numElementslin += '  ' + str(numElement) + '  '
    numElementslin += '\n'

    # Coordinate type (Cartesian in this implementation)
    middlelines = 'Cartesian\n'

    # Atomic coordinates (Cartesian)
    positionlines = ''
    for i in np.arange(len(atomslist)):
        positionlines += '  {:.6f}  {:.6f}  {:.6f}\n'.format(
            atomslist[i].coord[0],
            atomslist[i].coord[1],
            atomslist[i].coord[2]
        )

    # Combine all parts into final content
    content = (firstlin + secondlin + vBoxlines + elementlin +
               numElementslin + middlelines + positionlines)

    # Write to file
    with open(file=poscarname, mode='wt', encoding='UTF-8') as f_obj:
        f_obj.write(content)



def distenceCalc(coord2, coord1, vBox):
    """
    Calculate the distance between two atoms considering periodic boundary conditions.

    This function computes the minimum image distance between two points in a
    periodic simulation box. The box is assumed to be orthorhombic (aligned with
    Cartesian axes) and the vectors are extracted from the boxGET function.

    Parameters:
        coord2 (array-like): Position vector of atom 2 (x, y, z).
        coord1 (array-like): Position vector of atom 1 (x, y, z).
        vBox (list of np.ndarray): Simulation box vectors as returned by boxGET().
                                   Each vector is a 3-element float64 array.
                                   The diagonal elements vBox[i][i] represent the
                                   box lengths along each axis.

    Returns:
        float: The minimum image distance between coord1 and coord2.
    """
    # Calculate differences with minimum image convention for each Cartesian component

    # x-component
    if coord2[0] - coord1[0] > vBox[0][0] / 2:
        delta_coord_0 = coord2[0] - coord1[0] - vBox[0][0]
    elif coord1[0] - coord2[0] > vBox[0][0] / 2:
        delta_coord_0 = coord2[0] - coord1[0] + vBox[0][0]
    else:
        delta_coord_0 = coord2[0] - coord1[0]

    # y-component
    if coord2[1] - coord1[1] > vBox[1][1] / 2:
        delta_coord_1 = coord2[1] - coord1[1] - vBox[1][1]
    elif coord1[1] - coord2[1] > vBox[1][1] / 2:
        delta_coord_1 = coord2[1] - coord1[1] + vBox[1][1]
    else:
        delta_coord_1 = coord2[1] - coord1[1]

    # z-component
    if coord2[2] - coord1[2] > vBox[2][2] / 2:
        delta_coord_2 = coord2[2] - coord1[2] - vBox[2][2]
    elif coord1[2] - coord2[2] > vBox[2][2] / 2:
        delta_coord_2 = coord2[2] - coord1[2] + vBox[2][2]
    else:
        delta_coord_2 = coord2[2] - coord1[2]

    delta_coord = np.array([delta_coord_0, delta_coord_1, delta_coord_2])
    return np.linalg.norm(delta_coord)


import numpy as np

def direct2cart(vBox, coorddirect: list[np.float64]):
    """
    Transform atomic coordinates from direct (fractional) to Cartesian format.

    Parameters:
        vBox (list of np.ndarray): Simulation box vectors as a list of three
                                   3-element float64 arrays. Assumes vectors are
                                   arranged as rows (v1, v2, v3).
        coorddirect (list[np.float64]): Fractional coordinates of an atom as a
                                        list or array of three floats.

    Returns:
        np.ndarray: Cartesian coordinates as a 3-element float64 array.
    """
    coordcart = np.dot(coorddirect, vBox)
    return coordcart

def cart2direct(vBox, coordcart):
    """
    Transform atomic coordinates from Cartesian to direct (fractional) format.

    Parameters:
        vBox (list of np.ndarray): Simulation box vectors as a list of three
                                   3-element float64 arrays. Assumes vectors are
                                   arranged as rows (v1, v2, v3).
        coordcart (array-like): Cartesian coordinates of an atom (x, y, z).

    Returns:
        np.ndarray: Fractional coordinates as a 3-element float64 array.
    """
    coorddirect = np.dot(coordcart, np.linalg.inv(vBox))
    return coorddirect

# Calculate the parameter "d" (atomic displacements) as defined in the article.
def d_calc(atomslist_1: list[atom], atomslist_2: list[atom], vBox):
    """
    Calculate the atomic displacement parameter 'd' between two atomic configurations.

    This function computes the root mean square displacement of atoms between
    the unrelaxed and fully relaxed vacancy systems, as defined in the referenced article.

    Parameters:
        atomslist_1 (list[atom]): List of atom objects for the vacancy system without relaxation.
        atomslist_2 (list[atom]): List of atom objects for the same vacancy system after full relaxation.
        vBox (list of np.ndarray): Simulation box vectors as returned by boxGET().

    Returns:
        float: The displacement parameter d, calculated as the square root of the sum
               of squared distances between corresponding atoms in the two configurations.
    """
    d_2 = 0
    for i in np.arange(len(atomslist_1)):
        d_2 += distenceCalc(atomslist_1[i].coord, atomslist_2[i].coord, vBox) ** 2
    d = np.sqrt(d_2)
    return d

def extract_forces_before_iteration(outcar_path):
    """
    Extract the forces on the first atom from all "TOTAL-FORCE (eV/Angst)" blocks
    that appear before the first occurrence of "Iteration      2(   1)" in the OUTCAR file.

    Parameters:
        outcar_path (str): Path to the OUTCAR file.

    Returns:
        np.ndarray: An array of shape (N, 3) where N is the number of force blocks
                    meeting the condition. Each row contains the force components
                    (Fx, Fy, Fz) of the first atom in that block.
                    Returns an empty array np.array([]) if no force blocks are found.
    """
    forces_first_atom = []  # Store forces of the first atom from each block

    try:
        with open(outcar_path, 'r') as f:
            lines = f.readlines()
    except FileNotFoundError:
        raise FileNotFoundError(f"File not found: {outcar_path}")
    except Exception as e:
        raise IOError(f"Error reading file: {e}")

    stop_flag = False
    i = 0
    n_lines = len(lines)

    while i < n_lines:
        line = lines[i]

        # Check for the stop marker
        if 'Iteration      2(   1)' in line:
            stop_flag = True
            i += 1
            continue

        if stop_flag:
            i += 1
            continue

        # Locate the start of a force block
        if 'TOTAL-FORCE (eV/Angst)' in line:
            # Skip the header line
            i += 1
            if i >= n_lines:
                break

            # Skip possible separator lines (e.g., "----------------------------------------------------")
            while i < n_lines and lines[i].strip().startswith('---'):
                i += 1
            if i >= n_lines:
                break

            # Read data lines
            data_lines = []
            while i < n_lines:
                data_line = lines[i].strip()
                # An empty line or the next header ends the current block
                if not data_line or 'TOTAL-FORCE' in data_line or 'Iteration' in data_line:
                    # Step back one line because this line belongs to the next block or marker,
                    # to be handled in the outer loop.
                    i -= 1
                    break
                data_lines.append(data_line)
                i += 1

            # Parse data lines, extract forces for the first atom (last three columns)
            for j, data_line in enumerate(data_lines):
                parts = data_line.split()
                # Expect at least 6 numbers (3 coordinates + 3 forces), possibly more (e.g., index)
                if len(parts) >= 6:
                    try:
                        fx, fy, fz = float(parts[-3]), float(parts[-2]), float(parts[-1])
                    except ValueError:
                        continue  # Skip lines that cannot be converted
                    if j == 0:  # Take only the first atom
                        forces_first_atom.append([fx, fy, fz])
                        break  # First atom done, move to next block
        i += 1

    return np.array(forces_first_atom)

