#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import importlib
import os
import sys

import numpy as np
import simtk.openmm as mm
import simtk.unit as u
from md_utils import sizeof_fmt, plot_data
from scipy import ndimage
from simtk.openmm.app import PDBFile, ForceField, Simulation, PDBReporter, DCDReporter, StateDataReporter


def add_funnel(img, mask_n):
    brzegi = ndimage.distance_transform_edt(mask_n)
    brzegi = (brzegi - brzegi.min()) / (brzegi.max() - brzegi.min()) * 1
    return img + brzegi


def standardize_image(img):
    """returns image (0, -1)"""
    return - (img - img.min()) / (img.max() - img.min())


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('config', help="config file (python module)")
    args = parser.parse_args()

    settings_file = args.config
    if settings_file.endswith('.py'):
        settings_file = settings_file[:-3]
    if '.' in settings_file:
        sys.exit('settings file must not contain dot in its file name (except .py extension)')

    sys.path.append(os.getcwd())
    cfg = importlib.import_module(settings_file)
    sys.path.pop()

    print("Initial setup...")
    total_time = cfg.N_STEPS * cfg.TIME_STEP
    if cfg.RANDOM_SEED == 0:
        random_seed = np.random.randint(2147483647)
    else:
        random_seed = cfg.RANDOM_SEED
    reporting_to_screen_freq = max(1, int(round(cfg.N_STEPS / cfg.N_OF_STATE_READS_REPORTED_TO_SCREEN)))
    reporting_to_file_freq = max(1, int(round(cfg.N_STEPS / cfg.N_OF_STATE_READS_REPORTED_TO_FILE)))
    trajectory_freq = max(1, int(round(cfg.N_STEPS / cfg.TRAJECTORY_FRAMES)))

    print("   OpenMM version:                  {}".format(mm.__version__))
    print("   Number of steps:                 {} steps".format(cfg.N_STEPS))
    print("   Time step:                       {}".format(cfg.TIME_STEP))
    print("   Temperature:                     {}".format(cfg.TEMPERATURE))
    print("   Total simulation time:           {}".format(total_time.in_units_of(u.nanoseconds)))
    print("   Number of state reads:           {} reads".format(cfg.N_OF_STATE_READS_REPORTED_TO_SCREEN))
    print("   State reporting to screen every: {} step".format(reporting_to_screen_freq))
    print("   State reporting to file every:   {} step".format(reporting_to_file_freq))
    print("   Number of trajectory frames:     {} frames".format(cfg.TRAJECTORY_FRAMES))
    print("   Trajectory frame every:          {} step".format(trajectory_freq))
    print('   Random seed:', random_seed)
    print()

    print("Loading initial structure:\n  {}".format(cfg.INITIAL_STRUCTURE_FILENAME))
    pdb = PDBFile(cfg.INITIAL_STRUCTURE_FILENAME)
    print("Loading forcefield file:\n  {}".format(cfg.FORCEFIELD_FILE))
    forcefield = ForceField(cfg.FORCEFIELD_FILE)
    unmatched_residues = forcefield.getUnmatchedResidues(pdb.topology)
    print("Building system...")
    system = forcefield.createSystem(pdb.topology)

    if cfg.USE_RESTRAINTS:
        print('Loading restraints...')
        flat_bottom_force = mm.CustomBondForce('step(r-r0) * (k/2) * (r-r0)^2')
        flat_bottom_force.addPerBondParameter('r0')
        flat_bottom_force.addPerBondParameter('k')
        system.addForce(flat_bottom_force)

        with open(cfg.RESTRAINTS_FILE) as input_file:
            counter = 0
            for line in input_file:
                columns = line.split()
                atom_index_i = int(columns[0][1:]) - 1
                atom_index_j = int(columns[1][1:]) - 1
                try:
                    r0 = float(columns[3])
                    k = float(columns[4])
                except IndexError:
                    r0 = cfg.R0_PARAM
                    k = cfg.K_PARAM
                flat_bottom_force.addBond(atom_index_i, atom_index_j, [r0, k])
                counter += 1
        print("  {} restraints added.".format(counter))

    if cfg.USE_SPHERICAL_CONTAINER:
        container_force = mm.CustomExternalForce(
            '{}*max(0, r-{})^2; r=sqrt((x-{})^2+(y-{})^2+(z-{})^2)'.format(cfg.SPHERICAL_CONTAINER_SCALE,
                                                                           cfg.SPHERICAL_CONTAINER_RADIUS,
                                                                           cfg.SPHERICAL_CONTAINER_X,
                                                                           cfg.SPHERICAL_CONTAINER_Y,
                                                                           cfg.SPHERICAL_CONTAINER_Z,
                                                                           ))
        system.addForce(container_force)
        for i in range(system.getNumParticles()):
            container_force.addParticle(i, [])
        print("  Spherical container with radius {} nm and force {} applied.".format(cfg.SPHERICAL_CONTAINER_SCALE,
                                                                                     cfg.SPHERICAL_CONTAINER_RADIUS))
    else:
        print("  Spherical container will NOT be used.")

    if cfg.USE_EXTERNAL_FIELD:
        print('Loading external field...')
        size = os.stat(cfg.EXTERNAL_FIELD_FILE).st_size
        print("   Reading {} file ({})...".format(cfg.EXTERNAL_FIELD_FILE, sizeof_fmt(size)))
        img = np.load(cfg.EXTERNAL_FIELD_FILE)
        print("   Array of shape {} loaded.".format(img.shape))
        print("   Number of values: {}".format(img.size))
        print("   Min: {}".format(np.min(img)))
        print("   Max: {}".format(np.max(img)))
        if cfg.F_NORMALIZE:
            print('   [INFO] Field will be normalized to [0, -1]')
            img = standardize_image(img)
        print(f'   [INFO] IMG min = {np.min(img)}, max = {np.max(img)}')
        print(f'   [INFO] Adding funnel like border to image')
        maska_p = (img < -0.1)
        maska_n = np.logical_not(maska_p)
        img = add_funnel(img, maska_n)
        print("  Creating a force based on density...")
        voxel_size = np.array(cfg.F_VOXEL_SIZE)
        real_size = img.shape * voxel_size
        density_fun_args = {
            'xsize': img.shape[2],
            'ysize': img.shape[1],
            'zsize': img.shape[0],
            'values': img.flatten().astype(np.float64),
            'xmin': 0 * u.angstrom - 0.5 * voxel_size[0],
            'ymin': 0 * u.angstrom - 0.5 * voxel_size[1],
            'zmin': 0 * u.angstrom - 0.5 * voxel_size[2],
            'xmax': (img.shape[0] - 1) * voxel_size[0] + 0.5 * voxel_size[0],
            'ymax': (img.shape[1] - 1) * voxel_size[1] + 0.5 * voxel_size[1],
            'zmax': (img.shape[2] - 1) * voxel_size[2] + 0.5 * voxel_size[2],
        }

        print(f'   [INFO] Voxel size: ({cfg.F_VOXEL_SIZE[0]}, {cfg.F_VOXEL_SIZE[1]}, {cfg.F_VOXEL_SIZE[2]})')
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
        field_force.addGlobalParameter('ksi', cfg.F_SCALING_FACTOR)
        print("  Adding force to the system...")
        for i in range(system.getNumParticles()):
            field_force.addBond([i], [])
        system.addForce(field_force)
    else:
        print("External force (density) will NOT be used.")

    print("Integrator initialization...")
    if cfg.INTEGRATOR_TYPE == "langevin":
        integrator = mm.LangevinIntegrator(cfg.TEMPERATURE, cfg.FRICTION_COEFF, cfg.TIME_STEP)
        integrator.setRandomNumberSeed(random_seed)
    elif cfg.INTEGRATOR_TYPE == "verlet":
        integrator = mm.VerletIntegrator(cfg.TIME_STEP)
    else:
        sys.exit("Integrator initialization error!")

    print("Setting up simulation...")
    simulation = Simulation(pdb.topology, system, integrator)
    simulation.context.setPositions(pdb.positions)
    if cfg.MINIMIZE:
        print('  Energy minimizing...')
        simulation.minimizeEnergy(tolerance=0.01 * u.kilojoules_per_mole)
        if cfg.MINIMIZED_FILE:
            if cfg.MINIMIZED_FILE == ' ':
                base, _ = os.path.splitext(cfg.INITIAL_STRUCTURE_FILENAME)
                cfg.MINIMIZED_FILE = f'{base}_min.pdb'
            print('  Saving minimized structure in {}'.format(cfg.MINIMIZED_FILE))
            state = simulation.context.getState(getPositions=True)
            PDBFile.writeFile(pdb.topology, state.getPositions(), open(cfg.MINIMIZED_FILE, 'w'))
    if cfg.SET_INITIAL_VELOCITIES:
        print(f"Setting up initial velocities at temperature {cfg.TEMPERATURE}")
        simulation.context.setVelocitiesToTemperature(cfg.TEMPERATURE)

    print('Setting up reporters...')

    if cfg.TRAJECTORY_FILENAME_PDB:
        if cfg.TRAJECTORY_FILENAME_PDB == ' ':
            base, _ = os.path.splitext(cfg.INITIAL_STRUCTURE_FILENAME)
            cfg.TRAJECTORY_FILENAME_PDB = f'{base}_trj.pdb'
        simulation.reporters.append(PDBReporter(cfg.TRAJECTORY_FILENAME_PDB, trajectory_freq))
    if cfg.TRAJECTORY_FILENAME_DCD:
        if cfg.TRAJECTORY_FILENAME_DCD == ' ':
            base, _ = os.path.splitext(cfg.INITIAL_STRUCTURE_FILENAME)
        simulation.reporters.append(DCDReporter(cfg.TRAJECTORY_FILENAME_DCD, trajectory_freq))
    simulation.reporters.append(StateDataReporter(sys.stdout, reporting_to_screen_freq,
                                                  step=True, potentialEnergy=True, kineticEnergy=False,
                                                  totalEnergy=False,
                                                  temperature=False))

    simulation.reporters.append(StateDataReporter(cfg.STATE_FILE_NAME, reporting_to_file_freq,
                                                  step=True, potentialEnergy=True, kineticEnergy=False,
                                                  totalEnergy=False,
                                                  temperature=False))
    if cfg.RUN_SIMULATION:
        print('Running simulation...')
        simulation.step(cfg.N_STEPS)
        if cfg.PLOT_DATA:
            plot_data(cfg.STATE_FILE_NAME, cfg.PLOT_FILE_NAME)
    print()
    print("Everything is done")


if __name__ == '__main__':
    main()
