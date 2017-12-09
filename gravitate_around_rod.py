from __future__ import division
import hoomd
import hoomd.md
import numpy as np


hoomd.context.initialize("");
box_width = 200
positions = [[5*x, 0,0] for x in range(-10,10)]
for i in range(box_width):
  positions.append([i - box_width / 2,0,0])
num_particles = 20+box_width
types = ['R','A','B','C','D','E']
bond_types = []

# Place the type R central particles
# 2d box, 10x10
box = hoomd.data.boxdim(L=box_width,dimensions=3) # dimensions=2 forces Lz=1
snapshot = hoomd.data.make_snapshot(N=num_particles,
                                    box=box,
                                    particle_types=types,
                                    bond_types=bond_types);

print len(positions)
snapshot.particles.position[:] = positions
snapshot.particles.typeid[:] = [0 for x in positions]
for i in range(box_width):
  snapshot.particles.typeid[i+num_particles-box_width] = 5

snapshot.particles.moment_inertia[:] = [[0,0,10] for x in range(num_particles)]

# set gaussian random velocity for all proteins
snapshot.particles.velocity[:] = np.random.normal(0.0,
  np.sqrt(2 / 1.0), [snapshot.particles.N, 3]);

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
gauss.pair_coeff.set('A', 'A',epsilon=20, sigma=1.0);
gauss.pair_coeff.set('B', 'B',epsilon=20, sigma=1.0);
gauss.pair_coeff.set('C', 'C',epsilon=100, sigma=1.0);
gauss.pair_coeff.set('D', 'D',epsilon=100, sigma=1.0);
# these need to attract so proteins 
gauss.pair_coeff.set('A', 'D',epsilon=0, sigma=1.0);
gauss.pair_coeff.set('B', 'C',epsilon=0, sigma=1.0);
# everything else is neutral
gauss.pair_coeff.set('A', 'B',epsilon=0, sigma=1.0);
gauss.pair_coeff.set('A', 'C',epsilon=0, sigma=1.0);
gauss.pair_coeff.set('B', 'D',epsilon=0, sigma=1.0);
gauss.pair_coeff.set('C', 'D',epsilon=0, sigma=1.0);

# center of mass particles don't do anything
gauss.pair_coeff.set('R', 'R',epsilon=20, sigma=1.0);
gauss.pair_coeff.set('R', 'A',epsilon=0, sigma=1.0);
gauss.pair_coeff.set('R', 'B',epsilon=0, sigma=1.0);
gauss.pair_coeff.set('R', 'C',epsilon=0, sigma=1.0);
gauss.pair_coeff.set('R', 'D',epsilon=0, sigma=1.0);

gauss.pair_coeff.set('E', 'R',epsilon=0, sigma=1.0);
gauss.pair_coeff.set('E', 'A',epsilon=-200, sigma=1.0);
gauss.pair_coeff.set('E', 'B',epsilon=-200, sigma=1.0);
gauss.pair_coeff.set('E', 'C',epsilon=-100, sigma=1.0);
gauss.pair_coeff.set('E', 'D',epsilon=-100, sigma=1.0);
gauss.pair_coeff.set('E', 'E',epsilon=000, sigma=1.0);


rigidTest = hoomd.group.rigid_center();
hoomd.md.integrate.mode_standard(dt=0.001);
hoomd.md.integrate.nve(group=rigidTest);


hoomd.dump.gsd("trajectory.gsd", period=20, group=hoomd.group.all(), overwrite=True);

hoomd.run(5e4);

