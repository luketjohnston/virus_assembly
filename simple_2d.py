import hoomd
import hoomd.md
import numpy as np



hoomd.context.initialize("");


# will start with 4 particles
# box is 2d 
# two types of particles

class Protein:
  def __init__(self, positions, types, bonds):
    self.positions = np.copy(positions)
    self.origin_adjust_position() # adjust position
    self.bonds = bonds
    self.types = types
    self.num_particles = self.positions.shape[0]
    self.num_bonds = self.bonds.shape[0]
    # dimensions of particle
    self.particle_dims = np.amax(positions, axis=0)

  def shift_position(self, dx, dy, dz):
    for p in self.positions:
      p += np.array([dx, dy, dz])
  
  def origin_adjust_position(self):
    minx, miny, minz = positions[1]
    for p in positions:
      minx = min(p[0], minx)
      miny = min(p[0], miny)
      minz = min(p[0], minz)
    self.shift_position(-minx, -miny, -minz)



 



""" give a particle type and a number of particles, want to
- return a box that will contain all the particles
- get particle types and bond types and number of particles for snapshot
- get positions for snapshot
- get typid list mapping particles to types
- get number of bonds
- get list of bonds
"""

def get_n_proteins(protein, n, separation, box_width):
  ppp = protein.num_particles # "particles per protein"
  bpp = protein.num_bonds # "bonds per protein"
  pdims = protein.particle_dims
  positions = np.zeros((n * ppp, 3))
  bonds = np.zeros((n * bpp, 2))
  adjustment = np.zeros(3)
  # first, let's set all the positions
  for i in range(n):
    positions[ppp*i:ppp*(i+1),:] = (protein.positions + np.expand_dims(adjustment,0))
    adjustment[0] += pdims[0] + separation
    if adjustment[0] + pdims[0] > box_width:
      adjustment[0] = 0 # reset x adjustment, exceeded box width
      adjustment[1] += pdims[1] + separation # increment y adjustment instead
      if adjustment[1] + pdims[1] > box_width:
        adjustment[1] = 0 # reset y adjustment, exceeded box width
        adjustment[2] += pdims[2] + separation # increment z adjustment instead
  # now, lets set all the bonds
  for i in range(n):
    bonds[bpp*i:bpp*(i+1)] = protein.bonds + ppp*i
  types = protein.types * n
  return positions, bonds, types

types = ['A','B']
positions = np.array([[0,0,0],[0,1,0]])
typeids = [0,1]
bonds = np.array([[0,1]])
protein = Protein(positions, typeids, bonds)

box_width = 30
num_proteins = 4
positions, bonds, typeids = get_n_proteins(protein, num_proteins, separation=2, box_width=box_width)
num_bonds = bonds.shape[0]
num_particles = positions.shape[0]

bond_types = ['polymer']
    


# 2d box, 10x10
box = hoomd.data.boxdim(L=box_width,dimensions=2) # dimensions=2 forces Lz=1
snapshot = hoomd.data.make_snapshot(N=num_particles,
                                    box=box,
                                    particle_types=types,
                                    bond_types=bond_types);


# put particles at s
print(snapshot.particles.position.shape)
# have 3d coordinates because the way hoomd handles 2d simulations is
# by just operating in 3d with a small z-axis (Lz of the above
# box is equal to 1)
snapshot.particles.position[:] = positions

snapshot.particles.typeid[:] = typeids

snapshot.bonds.resize(num_bonds)
snapshot.bonds.group[:] = bonds


# set gaussian random velocity for all particles
snapshot.particles.velocity[:] = np.random.normal(0.0,
  np.sqrt(0.8 / 1.0), [snapshot.particles.N, 3]);

# initialize hoomd with this snapshot
hoomd.init.read_snapshot(snapshot)

# this tells hoomd what method to use to find the particles that interact
# with eachother. Each interaction type will have a cutoff length, and particles
# farther than that cutoff length will not interact. the "cell" method 
# apparently performs well when there is only one interaction type.
nl = hoomd.md.nlist.cell();

# use gaussian interaction because it is mathematicaaly simple
gauss = hoomd.md.pair.gauss(r_cut=10, nlist=nl)

# same type particles repel
gauss.pair_coeff.set('A', 'A',epsilon=-100.0, sigma=1.0);
gauss.pair_coeff.set('B', 'B',epsilon=-100.0, sigma=1.0);
# opposites attract
gauss.pair_coeff.set('A', 'B',epsilon=1.0, sigma=1.0);

# "apply harmonic bonds between the directly bonded particles"
harmonic = hoomd.md.bond.harmonic();
harmonic.bond_coeff.set('polymer', k=100.0, r0=0.5);


hoomd.md.integrate.mode_standard(dt=0.01);
group_all = hoomd.group.all();
hoomd.md.integrate.nve(group=group_all);

hoomd.analyze.log(filename="log-output.log",
                  quantities=['potential_energy', 'temperature'],
                  period=500,
                  overwrite=True);

hoomd.dump.gsd("trajectory.gsd", period=20, group=group_all, overwrite=True);

hoomd.run(5e4);






