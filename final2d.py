from __future__ import division
import hoomd
import hoomd.md
import numpy as np


hoomd.context.initialize("");

positions = [[5*x, 0,0] for x in range(-10,10)]
num_particles = 20
types = ['R','A','B','C','D','E']
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

snapshot.particles.moment_inertia[:] = [[0,0,10] for x in range(-10,10)]

# set gaussian random velocity for all proteins
snapshot.particles.velocity[:] = np.random.normal(0.0,
  np.sqrt(0.8 / 1.0), [snapshot.particles.N, 3]);

# initialize hoomd with this snapshot
hoomd.init.read_snapshot(snapshot)

rigid = hoomd.md.constrain.rigid();

theta = np.arccos(1.5/9)
scale = 1 # wanted to be able to scale proteins up in size but couldn't get COM coods right
start_x = -2
start_y = -3
protein_edge_particles = []
edge_particle_types = []
# particles along first edge
for i in range(1 + 9*scale):
  x = np.cos(theta)*i + (start_x*scale)
  y = np.sin(theta)*i + (start_y*scale)
  protein_edge_particles.append([x,y,0])
  if i <=3:
    edge_particle_types.append('A')
  else:
    edge_particle_types.append('B')
# particles along opposite edge
for i in range(1 + 9*scale):
  x = -np.cos(theta)*i + (start_x + 4)*scale
  y = np.sin(theta)*i + (start_y)*scale
  protein_edge_particles.append([x,y,0])
  if i <= 3:
    edge_particle_types.append('D')
  else:
    edge_particle_types.append('C')
# particles along top edge
for i in range(1,1*scale):
  x = (start_x + 1.5)*scale + i
  y = (start_y + 9)*scale
  protein_edge_particles.append([x,y,0])
  edge_particle_types.append('E')
# particles along bottom edge
for i in range(1,4*scale):
  x = (start_x)*scale + i
  y = (start_y)*scale
  protein_edge_particles.append([x,y,0])
  edge_particle_types.append('E')
  



# trapezoid particles, adjusted so R is approximately COM (just eyeballing)
rigid.set_param('R',
                types=edge_particle_types,
                positions=protein_edge_particles);
rigid.create_bodies()


nl = hoomd.md.nlist.cell();
# use gaussian interaction because it is mathematicaaly simple
gauss = hoomd.md.pair.gauss(r_cut=2, nlist=nl)
# same type particles repel
gauss.pair_coeff.set('A', 'A',epsilon=10, sigma=0.5);
gauss.pair_coeff.set('B', 'B',epsilon=10, sigma=0.5);
gauss.pair_coeff.set('C', 'C',epsilon=10, sigma=0.5);
gauss.pair_coeff.set('D', 'D',epsilon=10, sigma=0.5);
# these need to attract so proteins as0le
gauss.pair_coeff.set('A', 'D',epsilon=-10, sigma=0.5);
gauss.pair_coeff.set('B', 'C',epsilon=-10, sigma=0.5);
# everything else repels
gauss.pair_coeff.set('A', 'B',epsilon=10, sigma=0.5);
gauss.pair_coeff.set('A', 'C',epsilon=10, sigma=0.5);
gauss.pair_coeff.set('B', 'D',epsilon=10, sigma=0.5);
gauss.pair_coeff.set('C', 'D',epsilon=10, sigma=0.5);

# R does nothing
gauss.pair_coeff.set('R', 'R',epsilon=0, sigma=2.0);
gauss.pair_coeff.set('R', 'A',epsilon=0, sigma=2.0);
gauss.pair_coeff.set('R', 'B',epsilon=0, sigma=2.0);
gauss.pair_coeff.set('R', 'C',epsilon=0, sigma=2.0);
gauss.pair_coeff.set('R', 'D',epsilon=0, sigma=2.0);

gauss.pair_coeff.set('E', 'E',epsilon=10, sigma=2.0);
gauss.pair_coeff.set('E', 'R',epsilon=10, sigma=2.0);
gauss.pair_coeff.set('E', 'A',epsilon=10, sigma=2.0);
gauss.pair_coeff.set('E', 'B',epsilon=10, sigma=2.0);
gauss.pair_coeff.set('E', 'C',epsilon=10, sigma=2.0);
gauss.pair_coeff.set('E', 'D',epsilon=10, sigma=2.0);


rigidTest = hoomd.group.rigid_center();
hoomd.md.integrate.mode_standard(dt=0.001);
hoomd.md.integrate.nve(group=rigidTest);


hoomd.dump.gsd("trajectory.gsd", period=20, group=hoomd.group.all(), overwrite=True);

hoomd.run(5e4);

