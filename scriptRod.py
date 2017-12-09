from __future__ import division
import hoomd
import hoomd.md
import numpy as np


hoomd.context.initialize("");

# Place the type R central particles
uc = hoomd.lattice.unitcell(N = 1,
                            a1 = [11, 0,   0],
                            a2 = [0,    11, 0],
                            a3 = [0,    0,   1],
                            dimensions = 3,
                            position = [[0,0,0]],
                            type_name = ['R'],
                            mass = [1.0],
                            moment_inertia = [[0,
                                               0,
                                               1/12*1.0*8**2]],
                            orientation = [[1, 0, 0, 0]]);
system = hoomd.init.create_lattice(unitcell=uc, n=[1,1,10]);

# Add consituent particles of type A and create the rods
system.particles.types.add('A');
rigid = hoomd.md.constrain.rigid();
rigid.set_param('R',
                types=['A']*8,
                positions=[(-1,0,0),(-2,0,0),(-3,0,0),(-4,0,0),
                           (1,0,0),(2,0,0),(3,0,0),(4,0,0)]);

rigid.create_bodies()

rigidTest = hoomd.group.rigid_center();
hoomd.md.integrate.langevin(group=rigidTest, kT=1.0, seed=42);

hoomd.dump.gsd("trajectory.gsd", period=20, group=rigidTest, overwrite=True);

hoomd.run(5e4);

