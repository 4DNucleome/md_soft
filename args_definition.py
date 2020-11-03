import datetime
import re
from dataclasses import dataclass
from math import pi
from typing import Union

import simtk
import simtk.openmm as mm
from simtk.unit import Quantity


@dataclass
class Arg(object):
    name: str
    help: str
    type: type
    default: Union[str, float, int, bool, Quantity, None]
    val: Union[str, float, int, bool, Quantity, None]


class ListOfArgs(list):
    quantity_regexp = re.compile(r'(?P<value>[-+]?\d+(?:\.\d+)?) ?(?P<unit>\w+)')

    def get_arg(self, name: str) -> Arg:
        """Stupid arg search in list of args"""
        name = name.upper()
        for i in self:
            if i.name == name:
                return i
        raise ValueError(f"No such arg: {name}")

    def __getattr__(self, item):
        return self.get_arg(item).val

    def parse_quantity(self, val: str) -> Union[Quantity, None]:
        if val == '':
            return None
        match_obj = self.quantity_regexp.match(val)
        value, unit = match_obj.groups()
        try:
            unit = getattr(simtk.unit, unit)
        except AttributeError:
            raise ValueError(f"I Can't recognise unit {unit} in expresion {val}. Example of valid quantity: 12.3 femtosecond.")
        return Quantity(value=float(value), unit=unit)

    def to_python(self):
        """Casts string args to ints, floats, bool..."""
        for i in self:
            if i.val == '':
                i.val = None
            elif i.name == "HR_K_PARAM":  # Workaround for complex unit
                i.val = Quantity(float(i.val), simtk.unit.kilojoule_per_mole / simtk.unit.nanometer ** 2)
            elif i.type == str:
                continue
            elif i.type == int:
                i.val = int(i.val)
            elif i.type == float:
                i.val = float(i.val)
            elif i.type == bool:
                if i.val.lower() in ['true', '1', 'y', 'yes']:
                    i.val = True
                elif i.val.lower() in ['false', '0', 'n', 'no']:
                    i.val = False
                else:
                    raise ValueError(f"Can't convert {i.val} into bool type.")
            elif i.type == Quantity:
                try:
                    i.val = self.parse_quantity(i.val)
                except AttributeError:
                    raise ValueError(f"Can't parse: {i.name} = {i.val}")
            else:
                raise ValueError(f"Can't parse: {i.name} = {i.val}")

    def get_complete_config(self) -> str:
        w = "####################\n"
        w += "#   Spring Model   #\n"
        w += "####################\n\n"
        w += "# This is automatically generated config file.\n"
        w += f"# Generated at: {datetime.datetime.now().isoformat()}\n\n"
        w += "# Notes:\n"
        w += "# Some fields require units. Units are represented as objects from simtk.units module.\n"
        w += "# Simple units are parsed directly. For example: \n"
        w += "# HR_R0_PARAM = 0.2 nanometer\n"
        w += "# But more complex units does not have any more sophisticated parser written, and will fail.'\n"
        w += "# In such cases the unit is fixed (and noted in comment), so please convert complex units manually if needed.\n"
        w += "# <float> and <int> types does not require any unit. Quantity require unit.\n\n"
        w += "# Default values does not mean valid value. In many places it's only a empty field that need to be filled.\n\n"

        w += '[Main]'
        for i in self:
            w += f'; {i.help}, type: {i.type.__name__}, default: {i.default}\n'
            if i.val is None:
                w += f'{i.name} = \n\n'
            else:
                if i.type == Quantity:
                    # noinspection PyProtectedMember
                    w += f'{i.name} = {i.val._value} {i.val.unit.get_name()}\n\n'
                else:
                    w += f'{i.name} = {i.val}\n\n'
        w = w[:-2]
        return w

    def write_config_file(self):
        auto_config_filename = 'config_auto.ini'
        with open(auto_config_filename, 'w') as f:
            f.write(self.get_complete_config())
        print(f"Automatically generated config file saved in {auto_config_filename}")


available_platforms = [mm.Platform.getPlatform(i).getName() for i in range(mm.Platform.getNumPlatforms())]

# Every single arguments must be listed here.
# Not all of them must be provided by user.
# Invalid arguments should rise ValueError.
# Default args ar overwritten by config.ini, and then they are overwritten by command line.
# Defaults value must be strings. They will be converted to python object later when ListOfArgs.to_python() will be called
args = ListOfArgs([
    # Platform settings
    Arg('PLATFORM', help=f"name of the platform. Available choices: {' '.join(available_platforms)}", type=str, default='', val=''),
    Arg('DEVICE', help="device index for CUDA or OpenCL (count from 0)", type=str, default='', val=''),

    # Input data
    Arg('INITIAL_STRUCTURE_PATH', help="Path to PDB file.", type=str, default='', val=''),
    Arg('FORCEFIELD_PATH', help="Path to XML file with forcefield.", type=str, default='', val=''),

    # Basic Polymer Bond
    Arg('POL_USE_HARMONIC_BOND', help="Use harmonic bond interaction.", type=bool, default='True', val='True'),
    Arg('POL_HARMONIC_BOND_R0', help="harmonic bond distance equilibrium constant", type=Quantity, default='0.1 nanometer', val='0.1 nanometer'),
    Arg('POL_HARMONIC_BOND_K', help="harmonic bond force constant (fixed unit: kJ/mol/nm^2)", type=float, default='300000.0', val='300000.0'),
    Arg('POL_USE_CONSTRAINTS', help="Use fixed bond length instead of harmonic interaction", type=bool, default='False', val='False'),
    Arg('POL_CONSTRAINT_DISTANCE', help="Fixed constraints length", type=Quantity, default='0.1 nanometer', val='0.1 nanometer'),

    # Basic polymer stiffness
    Arg('POL_USE_HARMONIC_ANGLE', help="Use harmonic angle interaction.", type=bool, default='True', val='True'),
    Arg('POL_HARMONIC_ANGLE_R0', help="harmonic angle distance equilibrium constant", type=float, default=str(pi), val=str(pi)),
    Arg('POL_HARMONIC_ANGLE_K', help="harmonic angle force constant (fixed unit: kJ/mol/radian^2)", type=float, default='10.0', val='10.0'),

    # Excluded Volume
    Arg('EV_USE_EXCLUDED_VOLUME', help="Use excluded volume.", type=bool, default='False', val='False'),
    Arg('EV_EPSILON', help="Epsilon parameter.", type=float, default='2.86', val='2.86'),
    Arg('EV_SIGMA', help="Sigma parameter.", type=Quantity, default='0.05 nanometer', val='0.05 nanometer'),

    # Harmonic restraints
    Arg('HR_USE_HARMONIC_RESTRAINTS', help="Use long range interactions or not", type=bool, default='False', val='False'),
    Arg('HR_USE_FLAT_BOTTOM_FORCE', help="Use flat bottom force instead of standard harmonic force.", type=bool, default='False', val='False'),
    Arg('HR_RESTRAINTS_PATH', help='Path to .rst file with indices', type=str, default='', val=''),
    Arg('HR_R0_PARAM', help='distance constant, this value will be used only if it is missing in rst file', type=Quantity, default='', val=''),
    Arg('HR_K_PARAM', help='force constant, this value will be used only if it is missing in rst file (Fixed unit: kilojoule_per_mole/nanometer**2 - only float number needed.)', type=float, default='', val=''),
    Arg('HR_K_SCALE', help='force constant multiplicator applied if you are using values from rst file.', type=float, default='1', val='1'),

    # Spherical container
    Arg('SC_USE_SPHERICAL_CONTAINER', help='Use Spherical container', type=bool, default='False', val='False'),
    Arg('SC_CENTER_X', help='Spherical container location x, fixed unit: nanometers', type=float, default='', val=''),
    Arg('SC_CENTER_Y', help='Spherical container location y, fixed unit: nanometers', type=float, default='', val=''),
    Arg('SC_CENTER_Z', help='Spherical container location z, fixed unit: nanometers', type=float, default='', val=''),
    Arg('SC_RADIUS', help='Spherical container radius, fixed unit: nanometers', type=float, default='30', val='30'),
    Arg('SC_SCALE', help='Spherical container scaling factor', type=float, default='1000', val='1000'),

    # Energy minimization
    Arg('MINIMIZE', help='should initial structure be minimized? - This is spring model main functionality.', type=bool, default='True', val='True'),
    Arg('MINIMIZED_FILE', help='If left empty result file will have name based on initial structure file name with _min.pdb ending.', type=str, default='', val=''),

    # Simulation parameters
    Arg('SIM_RUN_SIMULATION', help='Do you want to run MD simulation?', type=bool, default='False', val='False'),
    Arg('SIM_INTEGRATOR_TYPE', help='Alternative: langevin, verlet', type=str, default='verlet', val='verlet'),
    Arg('SIM_FRICTION_COEFF', help='Friction coefficient (Used only with langevin integrator)', type=float, default='', val=''),
    Arg('SIM_N_STEPS', help='Number of steps in MD simulation', type=int, default='', val=''),
    Arg('SIM_TIME_STEP', help='Time step (use time unit from simtk.unit module)', type=Quantity, default='', val=''),
    Arg('SIM_TEMP', help='Temperature (use temperature unit from simtk.unit module)', type=Quantity, default='', val=''),
    Arg('SIM_RANDOM_SEED', help='Random seed. Set to 0 for random seed.', type=int, default='0', val='0'),
    Arg('SIM_SET_INITIAL_VELOCITIES', help='Sets initial velocities based on Boltzmann distribution', type=bool, default='False', val='False'),

    # Trajectory settings
    Arg('TRJ_FRAMES', help='Number of trajectory frames to save.', type=int, default='2000', val='2000'),
    Arg('TRJ_FILENAME_DCD', help='Write trajectory in DCD file format, leave empty if you do not want to save.', type=str, default='', val=''),
    Arg('TRJ_FILENAME_PDB', help='Write trajectory in PDB file format, leave empty if you do not want to save.', type=str, default='', val=''),
    Arg('TRJ_LAST_FRAME_PDB', help='Write last frame of trajectory in PDB file format, leave empty if you do not want to save.', type=str, default='', val=''),

    # State reporting
    Arg('REP_STATE_N_SCREEN', help='Number of states reported on screen', type=int, default='20', val='20'),
    Arg('REP_STATE_N_FILE', help='Number of states reported to file screen', type=int, default='1000', val='1000'),
    Arg('REP_STATE_FILE_PATH', help='Filepath to save state. Leave empty if not needed.', type=str, default='state.csv', val='state.csv'),
    Arg('REP_STATE_FILE_H5_PATH', help="H5 file to save velocities. Leave empty if not needed.", type=str, default="state.h5", val='state.h5'),
    Arg('REP_PLOT_FILE_NAME', help='Filepath to save energy plot. Leave empty if not needed.', type=str, default='energy.pdf', val='energy.pdf'),

    # External Field parameters
    Arg('EF_USE_EXTERNAL_FIELD', help='External force', type=bool, default='False', val='False'),
    Arg('EF_PATH', help='npy file, that defines regular 3D grid with external field values', type=str, default='', val=''),
    Arg('EF_VOXEL_SIZE_X', help='External Field Voxel size X', type=Quantity, default='', val=''),
    Arg('EF_VOXEL_SIZE_Y', help='External Field Voxel size Y', type=Quantity, default='', val=''),
    Arg('EF_VOXEL_SIZE_Z', help='External Field Voxel size Z', type=Quantity, default='', val=''),
    Arg('EF_NORMALIZE', help='Should the field be normalized to [0;1]?', type=bool, default='False', val='False'),
    Arg('EF_SCALING_FACTOR', help='External field scaling factor', type=float, default='1.0', val='1.0'),
])
