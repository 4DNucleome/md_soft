#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2019 Michał Kadlof <m.kadlof@cent.uw.edu.pl>
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
# IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
# DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
# OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE
# OR OTHER DEALINGS IN THE SOFTWARE.

import argparse
import configparser
import os
import sys
from typing import List, Tuple

import numpy as np
import simtk
import simtk.openmm as mm
from mdtraj.reporters import HDF5Reporter
from scipy import ndimage
from simtk.openmm.app import PDBFile, ForceField, Simulation, PDBReporter, DCDReporter, StateDataReporter
from simtk.unit import Quantity

from args_definition import ListOfArgs
from md_utils import sizeof_fmt, plot_data


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


def add_forces_to_system(system: mm.System, args: ListOfArgs):
    """Helper function, that add forces to the system."""
    print("   Adding forces...")
    if args.POL_USE_HARMONIC_BOND:
        add_harmonic_bond(system, args)
    elif args.POL_USE_CONSTRAINTS:
        add_constraints(system, args)
    if args.POL_USE_HARMONIC_ANGLE:
        add_harmonic_angle(system, args)
    if args.EV_USE_EXCLUDED_VOLUME:
        add_excluded_volume(system, args)
    if args.HR_USE_HARMONIC_RESTRAINTS:
        add_harmonic_restraints(system, args)
    if args.SC_USE_SPHERICAL_CONTAINER:
        add_spherical_container(system, args)
    if args.EF_USE_EXTERNAL_FIELD:
        add_external_field(system, args)


def add_harmonic_bond(system: mm.System, args: ListOfArgs):
    print("      Adding harmonic bonds...")
    print(f"         r0 = {args.POL_HARMONIC_BOND_R0}")
    print(f"         k = {args.POL_HARMONIC_BOND_K} kJ/mol/nm^2")
    bond_force = mm.HarmonicBondForce()
    system.addForce(bond_force)
    counter = 0
    for i in range(system.getNumParticles() - 1):
        bond_force.addBond(i, i + 1, args.POL_HARMONIC_BOND_R0, args.POL_HARMONIC_BOND_K)
        counter += 1
    print(f"         {counter} harmonic bonds added.")


def add_constraints(system: mm.System, args: ListOfArgs):
    print("      Adding constraints...")
    print(f"         r = {args.POL_CONSTRAINT_DISTANCE}")
    counter = 0
    for i in range(system.getNumParticles() - 1):
        system.addConstraint(i, i + 1, args.POL_CONSTRAINT_DISTANCE)
        counter += 1
    print(f"         {counter} constraints added.")


def add_harmonic_angle(system: mm.System, args: ListOfArgs):
    print("      Adding harmonic angles...")
    print(f"         r0 = {args.POL_HARMONIC_ANGLE_R0}")
    print(f"         k = {args.POL_HARMONIC_ANGLE_K} kJ/mol/radian^2")
    bond_force = mm.HarmonicAngleForce()
    system.addForce(bond_force)
    counter = 0
    for i in range(system.getNumParticles() - 2):
        bond_force.addAngle(i, i + 1, i + 2, args.POL_HARMONIC_ANGLE_R0, args.POL_HARMONIC_ANGLE_K)
        counter += 1
    print(f"         {counter} harmonic bonds added.")


def add_excluded_volume(system: mm.System, args: ListOfArgs):
    print("      Adding excluded volume...")
    print(f"         epsilon = {args.EV_EPSILON}")
    print(f"         sigma = {args.EV_SIGMA}")
    ev_force = mm.CustomNonbondedForce('epsilon*((sigma1+sigma2)/r)^12')
    ev_force.addGlobalParameter('epsilon', defaultValue=args.EV_EPSILON)
    ev_force.addPerParticleParameter('sigma')
    system.addForce(ev_force)
    counter = 0
    for i in range(system.getNumParticles()):
        ev_force.addParticle([args.EV_SIGMA])
        counter += 1
    print(f"         {counter} ev interactions added.")


def add_harmonic_restraints(system: mm.System, args: ListOfArgs):
    """Restraints format is different here than in SM webservice.
    Example record: :10 :151\n"""
    print("      Adding harmonic restraints")
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
    print(f"         {counter} restraints added.")


def add_spherical_container(system: mm.System, args: ListOfArgs):
    print("      Adding spherical container...")
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
    print(f"         Spherical container added.")
    print(f"            radius: {args.SC_RADIUS} nm")
    print(f"            scale:  {args.SC_SCALE} ")
    print(f"            center: ({args.SC_CENTER_X}, {args.SC_CENTER_Y}, {args.SC_CENTER_Z})")


def add_external_field(system: mm.System, args: ListOfArgs):
    """Add external forcefield for image-driven modelling purposes."""
    print('      Adding external forcefield.')
    size = os.stat(args.EF_PATH).st_size
    print(f"   Reading {args.EF_PATH} file ({sizeof_fmt(size)})...")
    img = np.load(args.EF_PATH)
    print(f"   Array of shape {img.shape} loaded.")
    print(f"   Number of values: {img.size}")
    print(f"   Min: {np.min(img)}")
    print(f"   Max: {np.max(img)}")
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


def get_integrator(random_seed: int, args: ListOfArgs) -> mm.Integrator:
    """Helper function that returns requested integrator."""
    print("   Integrator initialization...")
    integrator = mm.VerletIntegrator(10 * simtk.unit.femtosecond)  # default integrator
    if args.SIM_RUN_SIMULATION:
        if args.SIM_INTEGRATOR_TYPE == "langevin":
            integrator = mm.LangevinIntegrator(args.SIM_TEMP, args.SIM_FRICTION_COEFF, args.SIM_TIME_STEP)
            integrator.setRandomNumberSeed(random_seed)
        elif args.SIM_INTEGRATOR_TYPE == "verlet":
            integrator = mm.VerletIntegrator(args.SIM_TIME_STEP)
    return integrator


def get_config() -> ListOfArgs:
    """This function prepares the list of arguments.
    At first List of args with defaults is read.
    Then it's overwritten by args from config file (ini file).
    In the end config is overwritten by argparse options."""

    print(f"Reading config...")
    from args_definition import args
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
    args.write_config_file()
    return args


def setup(args: ListOfArgs) -> Tuple[PDBFile, int, Simulation]:
    print("Initialization...")
    if args.SIM_RANDOM_SEED == 0:
        random_seed = np.random.randint(2147483647)
    else:
        random_seed = args.SIM_RANDOM_SEED
    print(f"   Loading initial structure: {args.INITIAL_STRUCTURE_PATH}")
    pdb = PDBFile(args.INITIAL_STRUCTURE_PATH)
    print(f"   Loading forcefield file:  {args.FORCEFIELD_PATH}")
    forcefield = ForceField(args.FORCEFIELD_PATH)
    print("   Building system...")
    system = forcefield.createSystem(pdb.topology)
    add_forces_to_system(system, args)
    integrator = get_integrator(random_seed, args)
    print("   Setting up simulation...")
    simulation = Simulation(pdb.topology, system, integrator)
    simulation.context.setPositions(pdb.positions)
    return pdb, random_seed, simulation


def minimize_energy(pdb: PDBFile, simulation: Simulation, args: ListOfArgs):
    if args.MINIMIZE:
        print('Energy minimizing...')
        simulation.minimizeEnergy(tolerance=0.01 * simtk.unit.kilojoules_per_mole)
        if not args.MINIMIZED_FILE:
            base, _ = os.path.splitext(args.INITIAL_STRUCTURE_PATH)
            minimized_file_name = f'{base}_min.pdb'
        else:
            minimized_file_name = args.MINIMIZED_FILE.val
        print(f'  Saving minimized structure in {minimized_file_name}')
        state = simulation.context.getState(getPositions=True)
        PDBFile.writeFile(pdb.topology, state.getPositions(), open(minimized_file_name, 'w'))


def run_md_simulation(random_seed, simulation, args):
    if args.SIM_RUN_SIMULATION:
        print("Running simulation...")
        if args.SIM_SET_INITIAL_VELOCITIES:
            print(f"   Setting up initial velocities at temperature {args.SIM_TEMP}")
            simulation.context.setVelocitiesToTemperature(args.SIM_TEMP, random_seed)
        reporting_to_screen_freq = max(1, int(round(args.SIM_N_STEPS / args.REP_STATE_N_SCREEN)))
        reporting_to_file_freq = max(1, int(round(args.SIM_N_STEPS / args.REP_STATE_N_FILE)))
        trajectory_freq = max(1, int(round(args.SIM_N_STEPS / args.TRJ_FRAMES)))

        total_time = args.SIM_N_STEPS * args.SIM_TIME_STEP
        print("   Number of steps:                 {} steps".format(args.SIM_N_STEPS))
        print("   Time step:                       {}".format(args.SIM_TIME_STEP))
        print("   Temperature:                     {}".format(args.SIM_TEMP))
        print("   Total simulation time:           {}".format(total_time.in_units_of(simtk.unit.nanoseconds)))
        print("   Number of state reads:           {} reads".format(args.REP_STATE_N_SCREEN))
        print("   State reporting to screen every: {} step".format(reporting_to_screen_freq))
        print("   State reporting to file every:   {} step".format(reporting_to_file_freq))
        print("   Number of trajectory frames:     {} frames".format(args.TRJ_FRAMES))
        print("   Trajectory frame every:          {} step".format(trajectory_freq))
        print("   Trajectory frame every:          {}".format(trajectory_freq * args.SIM_TIME_STEP))
        print('   Random seed:', random_seed)
        print()
        if args.TRJ_FILENAME_PDB:
            simulation.reporters.append(PDBReporter(args.TRJ_FILENAME_PDB, trajectory_freq))
        if args.TRJ_FILENAME_DCD:
            simulation.reporters.append(DCDReporter(args.TRJ_FILENAME_DCD, trajectory_freq))
        simulation.reporters.append(StateDataReporter(sys.stdout, reporting_to_screen_freq,
                                                      step=True, progress=True, potentialEnergy=True,
                                                      totalSteps=args.SIM_N_STEPS))
        if args.REP_STATE_FILE_PATH:
            simulation.reporters.append(StateDataReporter(args.REP_STATE_FILE_PATH, reporting_to_file_freq,
                                                          step=True, potentialEnergy=True))
        if args.REP_STATE_FILE_H5_PATH:
            simulation.reporters.append(HDF5Reporter(args.REP_STATE_FILE_H5_PATH, reporting_to_file_freq, velocities=True))

        print('Running simulation...')
        simulation.step(args.SIM_N_STEPS)
        if args.REP_PLOT_FILE_NAME:
            plot_data(args.REP_STATE_FILE_PATH, args.REP_PLOT_FILE_NAME)


def main():
    print(f"Spring model - in house version.")
    print(f"OpenMM version: {mm.__version__}")
    args = get_config()
    pdb, random_seed, simulation = setup(args)
    minimize_energy(pdb, simulation, args)
    run_md_simulation(random_seed, simulation, args)
    print()
    print("Everything is done")


if __name__ == '__main__':
    main()
