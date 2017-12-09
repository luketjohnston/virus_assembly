from __future__ import division
import hoomd
import hoomd.md
import numpy as np


hoomd.context.initialize("");

positions = [[5*x, 0,0] for x in range(-10,10)]
num_particles = 20
types = ['R','A','B','C','D']
bond_types = []
box_width = 200

# Place the type R central particles
# 2d box, 10x10
box = hoomd.data.boxdim(L=box_width,dimensions=2) # dimensions=2 forces Lz=1
snapshot = hoomd.data.make_snapshot(N=num_particles,
                                    box=box,
                                    particle_types=types,
                                    bond_types=bond_types);


snapshot.particles.position[:] = positions
snapshot.particles.typeid[:] = [0 for x in range(-10,10)]

# set gaussian random velocity for all proteins
snapshot.particles.velocity[:] = np.random.normal(0.0,
  np.sqrt(0.8 / 1.0), [snapshot.particles.N, 3]);

# initialize hoomd with this snapshot
hoomd.init.read_snapshot(snapshot)

rigid = hoomd.md.constrain.rigid();
rigid.set_param('R',
                types=['A','B','C','D'],
                positions=[(0,0,0),(1.5,9,0),
                           (2.5,9,0),(4,0,0)]);
rigid.create_bodies()


nl = hoomd.md.nlist.cell();
# use gaussian interaction because it is mathematicaaly simple
gauss = hoomd.md.pair.gauss(r_cut=10, nlist=nl)
# same type particles repel
gauss.pair_coeff.set('A', 'A',epsilon=-100, sigma=1.0);
gauss.pair_coeff.set('B', 'B',epsilon=-100, sigma=1.0);
gauss.pair_coeff.set('C', 'C',epsilon=-100, sigma=1.0);
gauss.pair_coeff.set('D', 'D',epsilon=-100, sigma=1.0);
# these need to attract so proteins as0le
gauss.pair_coeff.set('A', 'D',epsilon=100, sigma=1.0);
gauss.pair_coeff.set('B', 'C',epsilon=100, sigma=1.0);
# everything else repels
gauss.pair_coeff.set('A', 'B',epsilon=-100, sigma=1.0);
gauss.pair_coeff.set('A', 'C',epsilon=-100, sigma=1.0);
gauss.pair_coeff.set('B', 'D',epsilon=-100, sigma=1.0);
gauss.pair_coeff.set('C', 'D',epsilon=-100, sigma=1.0);

# center of mass particles don't do anything
gauss.pair_coeff.set('R', 'R',epsilon=0, sigma=1.0);
gauss.pair_coeff.set('R', 'A',epsilon=0, sigma=1.0);
gauss.pair_coeff.set('R', 'B',epsilon=0, sigma=1.0);
gauss.pair_coeff.set('R', 'C',epsilon=0, sigma=1.0);
gauss.pair_coeff.set('R', 'D',epsilon=0, sigma=1.0);


rigidTest = hoomd.group.rigid_center();
hoomd.md.integrate.mode_standard(dt=0.001);
hoomd.md.integrate.nve(group=rigidTest);


hoomd.dump.gsd("trajectory.gsd", period=20, group=hoomd.group.all(), overwrite=True);

hoomd.run(5e4);

