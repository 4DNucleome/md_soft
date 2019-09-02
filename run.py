#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import configparser
import datetime
import os
import re
import sys
from dataclasses import dataclass
from typing import List, Tuple, Union

import numpy as np
import simtk
import simtk.openmm as mm
from scipy import ndimage
from simtk.openmm.app import PDBFile, ForceField, Simulation, PDBReporter, DCDReporter, StateDataReporter
from simtk.unit import Quantity

from md_utils import sizeof_fmt, plot_data


@dataclass
class Arg(object):
    name: str
    help: str
    type: type
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
            raise ValueError(f"I Can't recognise unit {unit}. Example of valid quantity: 12.3 femtosecond.")
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
                i.val = self.parse_quantity(i.val)
            else:
                raise ValueError(f"Can't parse value: {i.val}")

    def conf_file(self) -> str:
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
        w += "# In such cases the unit is fixed (and noted in comment), so please convert complex units manually if needed.\n\n"
        w += '[Main]'
        for i in self:
            w += f'; {i.help}\n'
            if i.val is None:
                w += f'{i.name} = \n\n'
            else:
                if i.type == Quantity:
                    w += f'{i.name} = {i.val._value} {i.val.unit.get_name()}\n\n'
                else:
                    w += f'{i.name} = {i.val}\n\n'
        w = w[:-2]
        return w


def add_funnel(img, mask_n):
    funnel = ndimage.distance_transform_edt(mask_n)
    funnel = (funnel - funnel.min(initial=None)) / (funnel.max(initial=None) - funnel.min(initial=None)) * 1
    return img + funnel


def standardize_image(img):
    """returns image (0, -1)"""
    return - (img - img.min()) / (img.max() - img.min())


def my_config_parser(config_parser: configparser.ConfigParser) -> List[Tuple[str, str]]:
    """Helper function that makes flat list arg name, and it's value from ConfigParser object."""
    sections = config_parser.sections()
    all_nested_fields = [dict(config_parser[s]) for s in sections]
    args_cp = []
    for section_fields in all_nested_fields:
        for name, value in section_fields.items():
            args_cp.append((name, value))
    return args_cp


def main():
    # Every single arguments must be listed here.
    # Arguments may be missing.
    # Invalid arguments should rise ValueError.
    # Default args ar overwritten by config.ini, and then they are overwritten by command line.
    # Defaults value must be strings. They will be converted to python object later when ListOfArgs.to_python() will be called
    args = ListOfArgs([
        Arg('INITIAL_STRUCTURE_PATH', help="Path to PDB file.", type=str, val=''),
        Arg('FORCEFIELD_PATH', help="Path to XML file with forcefield.", type=str, val=''),

        # Harmonic restraints
        Arg('HR_USE_HARMONIC_RESTRAINTS', help="Use long range interactions or not (True/False)", type=bool, val='False'),
        Arg('HR_USE_FLAT_BOTTOM_FORCE', help="Use flat bottom force instead of standard harmonic force. (True/False)", type=bool, val='False'),
        Arg('HR_RESTRAINTS_PATH', help='Path to .rst file with indices', type=str, val=''),
        Arg('HR_R0_PARAM', help='distance constant, this value will be used only if it is missing in rst file', type=Quantity, val=''),
        Arg('HR_K_PARAM', help='force constant, this value will be used only if it is missing in rst file (Fixed unit: kilojoule_per_mole/nanometer**2 - only float number needed.)', type=float, val=''),

        # Spherical container
        Arg('SC_USE_SPHERICAL_CONTAINER', help='Use Spherical container (True/False)', type=bool, val='False'),
        Arg('SC_CENTER_X', help='Spherical container location x', type=Quantity, val=''),
        Arg('SC_CENTER_Y', help='Spherical container location y', type=Quantity, val=''),
        Arg('SC_CENTER_Z', help='Spherical container location z', type=Quantity, val=''),
        Arg('SC_RADIUS', help='Spherical container radius', type=Quantity, val=''),
        Arg('SC_SCALE', help='Spherical container scaling factor', type=Quantity, val=''),

        # Energy minimization
        Arg('MINIMIZE', help='should initial structure be minimized? (True/False) - This is spring model main functionality.', type=bool, val='True'),
        Arg('MINIMIZED_FILE', help='If left empty result file will have name based on initial structure file name with _min.pdb ending.', type=str, val=''),

        # Simulation parameters
        Arg('SIM_RUN_SIMULATION', help='Do you want to run MD simulation? (True/False)', type=bool, val='False'),
        Arg('SIM_INTEGRATOR_TYPE', help='Alternative: langevin, verlet', type=str, val='verlet'),
        Arg('SIM_FRICTION_COEFF', help='Friction coefficient (Used only with langevin integrator)', type=float, val=''),
        Arg('SIM_N_STEPS', help='Number of steps in MD simulation', type=int, val=''),
        Arg('SIM_TIME_STEP', help='Time step (use time unit from simtk.unit module)', type=Quantity, val=''),
        Arg('SIM_TEMP', help='Temperature (use temperature unit from simtk.unit module)', type=Quantity, val=''),
        Arg('SIM_RANDOM_SEED', help='Random seed. Set to 0 for random seed.', type=int, val='0'),
        Arg('SIM_SET_INITIAL_VELOCITIES', help='Sets initial velocities based on Boltzmann distribution (True/False)', type=bool, val='True'),

        # Trajectory settings
        Arg('TRJ_FRAMES', help='Number of trajectory frames to save.', type=int, val='2000'),
        Arg('TRJ_FILENAME_DCD', help='Write trajectory in DCD file format, leave empty if you do not want to save.', type=str, val=''),
        Arg('TRJ_FILENAME_PDB', help='Write trajectory in PDB file format, leave empty if you do not want to save.', type=str, val=''),

        # State reporting
        Arg('REP_STATE_N_SCREEN', help='Number of states reported on screen', type=int, val='20'),
        Arg('REP_STATE_N_FILE', help='Number of states reported to file screen', type=int, val='1000'),
        Arg('REP_STATE_FILE_PATH', help='Filepath to save state. Leave empty if not needed.', type=str, val='state.csv'),
        Arg('REP_PLOT_FILE_NAME', help='Filepath to save energy plot. Leave empty if not needed.', type=str, val='energy.pdf'),

        # External Field parameters
        Arg('EF_USE_EXTERNAL_FIELD', help='External force', type=bool, val='False'),
        Arg('EF_PATH', help='npy file, that defines regular 3D grid with external field values', type=str, val=''),
        Arg('EF_VOXEL_SIZE_X', help='External Field Voxel size X', type=Quantity, val=''),
        Arg('EF_VOXEL_SIZE_Y', help='External Field Voxel size Y', type=Quantity, val=''),
        Arg('EF_VOXEL_SIZE_Z', help='External Field Voxel size Z', type=Quantity, val=''),
        Arg('EF_NORMALIZE', help='Should the field be normalized to [0;1]? (True/False)', type=bool, val='False'),
        Arg('EF_SCALING_FACTOR', help='External field scaling factor', type=float, val='1.0'),
    ])

    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument('-c', '--config_file', help="Specify config file (ini format)", metavar="FILE")
    for arg in args:
        arg_parser.add_argument(f"--{arg.name.lower()}", help=arg.help)
    args_ap = arg_parser.parse_args()  # args from argparse

    config_parser = configparser.ConfigParser()
    config_parser.read(args_ap.config_file)
    args_cp = my_config_parser(config_parser)

    # Override defaults args with values from config file
    for cp_arg in args_cp:
        name, value = cp_arg
        arg = args.get_arg(name)
        arg.val = value

    # Now again override args with values from command line.
    for ap_arg in args_ap.__dict__:
        if ap_arg != 'config_file':
            name, value = ap_arg, getattr(args_ap, ap_arg)
            if value is not None:
                arg = args.get_arg(name)
                arg.val = value

    args.to_python()

    print("   OpenMM version:                  {}".format(mm.__version__))
    auto_config_filename = 'config_auto.ini'
    with open(auto_config_filename, 'w') as f:
        f.write(args.conf_file())
    print(f"Automatically generated config file saved in {auto_config_filename}")

    print("Initial setup...")
    if args.SIM_RANDOM_SEED == 0:
        random_seed = np.random.randint(2147483647)
    else:
        random_seed = args.SIM_RANDOM_SEED

    print("Loading initial structure:\n  {}".format(args.INITIAL_STRUCTURE_PATH))
    pdb = PDBFile(args.INITIAL_STRUCTURE_PATH)
    print("Loading forcefield file:\n  {}".format(args.FORCEFIELD_PATH))
    forcefield = ForceField(args.FORCEFIELD_PATH)
    print("Building system...")
    system = forcefield.createSystem(pdb.topology)

    if args.HR_USE_HARMONIC_RESTRAINTS:
        print('Loading restraints...')
        if args.HR_USE_FLAT_BOTTOM_FORCE:
            contact_force = mm.CustomBondForce('step(r-r0) * (k/2) * (r-r0)^2')
            contact_force.addPerBondParameter('r0')
            contact_force.addPerBondParameter('k')
        else:
            contact_force = mm.HarmonicBondForce()
        system.addForce(contact_force)

        with open(args.HR_RESTRAINTS_PATH) as input_file:
            counter = 0
            for line in input_file:
                columns = line.split()
                atom_index_i = int(columns[0][1:]) - 1
                atom_index_j = int(columns[1][1:]) - 1
                try:
                    r0 = float(columns[3])
                    k = float(columns[4])
                except IndexError:
                    r0 = args.HR_R0_PARAM
                    k = args.HR_K_PARAM
                if args.HR_USE_FLAT_BOTTOM_FORCE:
                    contact_force.addBond(atom_index_i, atom_index_j, [r0, k])
                else:
                    contact_force.addBond(atom_index_i, atom_index_j, r0, k)
                counter += 1
        print("  {} restraints added.".format(counter))

    if args.SC_USE_SPHERICAL_CONTAINER:
        container_force = mm.CustomExternalForce(
            '{}*max(0, r-{})^2; r=sqrt((x-{})^2+(y-{})^2+(z-{})^2)'.format(args.SC_SCALE,
                                                                           args.SC_RADIUS,
                                                                           args.SC_CENTER_X,
                                                                           args.SC_CENTER_Y,
                                                                           args.SC_CENTER_Z,
                                                                           ))
        system.addForce(container_force)
        for i in range(system.getNumParticles()):
            container_force.addParticle(i, [])
        print("  Spherical container with radius {} nm and force {} applied.".format(args.SC_SCALE,
                                                                                     args.SC_RADIUS))
    else:
        print("  Spherical container will NOT be used.")

    if args.EF_USE_EXTERNAL_FIELD:
        print('Loading external field...')
        size = os.stat(args.EF_PATH).st_size
        print("   Reading {} file ({})...".format(args.EF_PATH, sizeof_fmt(size)))
        img = np.load(args.EF_PATH)
        print("   Array of shape {} loaded.".format(img.shape))
        print("   Number of values: {}".format(img.size))
        print("   Min: {}".format(np.min(img)))
        print("   Max: {}".format(np.max(img)))
        if args.EF_NORMALIZE:
            print('   [INFO] Field will be normalized to [0, -1]')
            img = standardize_image(img)
        print(f'   [INFO] IMG min = {np.min(img)}, max = {np.max(img)}')
        print(f'   [INFO] Adding funnel like border to image')
        mask_p = (img < -0.1)
        mask_n = np.logical_not(mask_p)
        img = add_funnel(img, mask_n)
        print("  Creating a force based on density...")
        voxel_size = np.array((args.EF_VOXEL_SIZE_X, args.EF_VOXEL_SIZE_Y, args.EF_VOXEL_SIZE_Z))
        real_size = img.shape * voxel_size
        density_fun_args = dict(
            xsize=img.shape[2],
            ysize=img.shape[1],
            zsize=img.shape[0],
            values=img.flatten().astype(np.float64),
            xmin=0 * simtk.unit.angstrom - 0.5 * voxel_size[0],
            ymin=0 * simtk.unit.angstrom - 0.5 * voxel_size[1],
            zmin=0 * simtk.unit.angstrom - 0.5 * voxel_size[2],
            xmax=(img.shape[0] - 1) * voxel_size[0] + 0.5 * voxel_size[0],
            ymax=(img.shape[1] - 1) * voxel_size[1] + 0.5 * voxel_size[1],
            zmax=(img.shape[2] - 1) * voxel_size[2] + 0.5 * voxel_size[2])

        print(f'   [INFO] Voxel size: ({args.EF_VOXEL_SIZE_X}, {args.EF_VOXEL_SIZE_Y}, {args.EF_VOXEL_SIZE_Z})')
        print(f'   [INFO] Real size (Shape * voxel size): ({real_size[0]}, {real_size[1]}, {real_size[2]})')
        print(
            f"   [INFO] begin coords: ({density_fun_args['xmin']}, {density_fun_args['ymin']}, {density_fun_args['zmin']})")
        print(
            f"   [INFO] end coords:   ({density_fun_args['xmax']}, {density_fun_args['ymax']}, {density_fun_args['zmax']})")
        center_x = (density_fun_args['xmax'] - density_fun_args['xmin']) / 2 + density_fun_args['xmin']
        center_y = (density_fun_args['ymax'] - density_fun_args['ymin']) / 2 + density_fun_args['ymin']
        center_z = (density_fun_args['zmax'] - density_fun_args['zmin']) / 2 + density_fun_args['zmin']
        print(f"   [INFO] Image central point: ({center_x}, {center_y}, {center_z}) ")
        field_function = mm.Continuous3DFunction(**density_fun_args)
        field_force = mm.CustomCompoundBondForce(1, 'ksi*fi(x1,y1,z1)')
        field_force.addTabulatedFunction('fi', field_function)
        field_force.addGlobalParameter('ksi', args.EF_SCALING_FACTOR)
        print("  Adding force to the system...")
        for i in range(system.getNumParticles()):
            field_force.addBond([i], [])
        system.addForce(field_force)
    else:
        print("External force (density) will NOT be used.")

    print("Integrator initialization...")
    integrator = mm.VerletIntegrator(10 * simtk.unit.femtosecond)  # default integrator
    if args.SIM_RUN_SIMULATION:
        if args.SIM_INTEGRATOR_TYPE == "langevin":
            integrator = mm.LangevinIntegrator(args.SIM_TEMP, args.SIM_FRICTION_COEFF, args.SIM_TIME_STEP)
            integrator.setRandomNumberSeed(random_seed)
        elif args.SIM_INTEGRATOR_TYPE == "verlet":
            integrator = mm.VerletIntegrator(args.SIM_TIME_STEP)

    print("Setting up simulation...")
    print('dupa')
    simulation = Simulation(pdb.topology, system, integrator)
    simulation.context.setPositions(pdb.positions)
    if args.MINIMIZE:
        print('  Energy minimizing...')
        simulation.minimizeEnergy(tolerance=0.01 * simtk.unit.kilojoules_per_mole)
        if args.MINIMIZED_FILE:
            if args.MINIMIZED_FILE is None:
                base, _ = os.path.splitext(args.INITIAL_STRUCTURE_PATH)
                args.MINIMIZED_FILE = f'{base}_min.pdb'
            print('  Saving minimized structure in {}'.format(args.MINIMIZED_FILE))
            state = simulation.context.getState(getPositions=True)
            PDBFile.writeFile(pdb.topology, state.getPositions(), open(args.MINIMIZED_FILE, 'w'))
    if args.SIM_SET_INITIAL_VELOCITIES:
        print(f"Setting up initial velocities at temperature {args.SIM_TEMP}")
        simulation.context.setVelocitiesToTemperature(args.SIM_TEMP)

    print('Setting up reporters...')

    if args.SIM_RUN_SIMULATION:
        reporting_to_screen_freq = max(1, int(round(args.SIM_N_STEPS / args.REP_STATE_N_SCREEN)))
        reporting_to_file_freq = max(1, int(round(args.SIM_N_STEPS / args.REP_STATE_N_FILE)))
        trajectory_freq = max(1, int(round(args.SIM_N_STEPS / args.TRJ_FRAMES)))

        total_time = args.SIM_N_STEPS * args.SIM_TIME_STEP
        print("   Number of steps:                 {} steps".format(args.SIM_N_STEPS))
        print("   Time step:                       {}".format(args.SIM_TIME_STEP))
        print("   Temperature:                     {}".format(args.TEMPERATURE))
        print("   Total simulation time:           {}".format(total_time.in_units_of(simtk.unit.nanoseconds)))
        print("   Number of state reads:           {} reads".format(args.REP_STATE_N_SCREEN))
        print("   State reporting to screen every: {} step".format(reporting_to_screen_freq))
        print("   State reporting to file every:   {} step".format(reporting_to_file_freq))
        print("   Number of trajectory frames:     {} frames".format(args.TRJ_FRAMES))
        print("   Trajectory frame every:          {} step".format(trajectory_freq))
        print('   Random seed:', random_seed)
        print()
        if args.TRJ_FILENAME_PDB:
            simulation.reporters.append(PDBReporter(args.TRJ_FILENAME_PDB, trajectory_freq))
        if args.TRAJECTORY_FILENAME_DCD:
            simulation.reporters.append(DCDReporter(args.TRJ_FILENAME_DCD, trajectory_freq))
        simulation.reporters.append(StateDataReporter(sys.stdout, reporting_to_screen_freq,
                                                      step=True, potentialEnergy=True, kineticEnergy=False,
                                                      totalEnergy=False,
                                                      temperature=False))

        simulation.reporters.append(StateDataReporter(args.REP_STATE_FILE_PATH, reporting_to_file_freq,
                                                      step=True, potentialEnergy=True, kineticEnergy=False,
                                                      totalEnergy=False,
                                                      temperature=False))

        print('Running simulation...')
        simulation.step(args.SIM_N_STEPS)
        if args.PLOT_DATA:
            plot_data(args.REP_STATE_FILE_PATH, args.REP_PLOT_FILE_NAME)
    print()
    print("Everything is done")


if __name__ == '__main__':
    main()
