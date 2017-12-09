from __future__ import division
import hoomd
import hoomd.md
import numpy as np



hoomd.context.initialize("");


# will start with 4 particles
# box is 2d 
# two types of particles

class Protein:
  def __init__(self, positions, types, bonds, bond_types):
    """ given positions of particles, their types (as an index, not 'A' or 'B' but 0 or 1), 
    and the bonds between particles (i.e. [[0,1]] is a bond between particle 0 and 1),
    constructs a protein that can be used in the method below this class"""
    self.positions = np.copy(positions)
    self.origin_adjust_position() # adjust position
    self.bonds = bonds
    self.types = types
    self.bond_types = bond_types
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



def get_n_proteins(protein, n, separation, box_width, two_d):
  """ given a protein made with the above class, a number n of proteins
  we want, a separation distance between each protein, and the width 
  of the simulation box we will be using, tiles the n proteins across the 
  box and returns 
  particles, a list of all the particles
  bonds, a list of all the bonds
  and typeids, a list of all the typids
  to be used with the hoomd simulation.

  2d is a bool for if we are doing 2d simulation (instead of 3d). I
  haven't tested 3d at all.
  """
  ppp = protein.num_particles # "particles per protein"
  bpp = protein.num_bonds # "bonds per protein"
  pdims = protein.particle_dims
  positions = np.zeros((n * ppp, 3))
  bonds = np.zeros((n * bpp, 2))
  adjustment = np.ones(3) * -box_width / 2.0
  if two_d:
    adjustment[2] = 0.
  # first, let's set all the positions
  for i in range(n):
    positions[ppp*i:ppp*(i+1),:] = (protein.positions + np.expand_dims(adjustment,0))
    adjustment[0] += pdims[0] + separation
    if adjustment[0] + pdims[0] > box_width / 2:
      adjustment[0] = -box_width/2 # reset x adjustment, exceeded box width
      adjustment[1] += pdims[1] + separation # increment y adjustment instead
      if adjustment[1] + pdims[1] > box_width / 2:
        adjustment[1] = -box_width/2 # reset y adjustment, exceeded box width
        adjustment[2] += pdims[2] + separation # increment z adjustment instead
  # now, lets set all the bonds
  all_bond_types = []
  for i in range(n):
    bonds[bpp*i:bpp*(i+1)] = protein.bonds + ppp*i
    all_bond_types += bond_types
  types = protein.types * n
  return positions, bonds, types, all_bond_types

"""
So now, this is the only area we need to change to experiment with different protein
structures. 
I imagine we will also want to experiment with different interaction types
and maybe bond types, but not sure if modularizing that would make it any easier.
"""
types = ['A','B','C','D']
positions = np.array([[0,0,0],[1.5,9,0],[2.5,9,0],[4,0,0]])
positions = np.array([[0,0,0],[0,1,0],[1,1,0],[1,0,0]])
typeids = [0,1,2,3]
bonds = np.array([[0,1],[1,2],[2,3],[3,0]])
bond_types = ['ab','bc','cd','da']
protein = Protein(positions, typeids, bonds, bond_types)

box_width = 30
num_proteins = 1
two_d = True
positions, bonds, typeids, bond_types = get_n_proteins(protein, num_proteins, separation=2, box_width=box_width, two_d=two_d)
num_bonds = bonds.shape[0]
num_particles = positions.shape[0]

bond_types = ['ab','bc','cd','da']
    


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
snapshot.bonds.types = bond_types


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
gauss.pair_coeff.set('A', 'A',epsilon=0, sigma=1.0);
gauss.pair_coeff.set('B', 'B',epsilon=0, sigma=1.0);
gauss.pair_coeff.set('C', 'C',epsilon=0, sigma=1.0);
gauss.pair_coeff.set('D', 'D',epsilon=0, sigma=1.0);
# these need to attract so proteins as0le
gauss.pair_coeff.set('A', 'D',epsilon=0, sigma=1.0);
gauss.pair_coeff.set('B', 'C',epsilon=0, sigma=1.0);
# everything else repels
gauss.pair_coeff.set('A', 'B',epsilon=0, sigma=1.0);
gauss.pair_coeff.set('A', 'C',epsilon=0, sigma=1.0);
gauss.pair_coeff.set('B', 'D',epsilon=0, sigma=1.0);
gauss.pair_coeff.set('C', 'D',epsilon=0, sigma=1.0);

#gauss.pair_coeff.set('A', 'A',epsilon=-1.0, sigma=1.0);
#gauss.pair_coeff.set('B', 'B',epsilon=-1.0, sigma=1.0);
#gauss.pair_coeff.set('C', 'C',epsilon=-1.0, sigma=1.0);
#gauss.pair_coeff.set('D', 'D',epsilon=-1.0, sigma=1.0);
## these need to attract so proteins assemble
#gauss.pair_coeff.set('A', 'D',epsilon=1.0, sigma=1.0);
#gauss.pair_coeff.set('B', 'C',epsilon=1.0, sigma=1.0);
## everything else repels
#gauss.pair_coeff.set('A', 'B',epsilon=-1.0, sigma=1.0);
#gauss.pair_coeff.set('A', 'C',epsilon=-1.0, sigma=1.0);
#gauss.pair_coeff.set('B', 'D',epsilon=-1.0, sigma=1.0);
#gauss.pair_coeff.set('C', 'D',epsilon=-1.0, sigma=1.0);


# "apply harmonic bonds between the directly bonded particles"
harmonic = hoomd.md.bond.harmonic();
#harmonic.bond_coeff.set('ab', k=100.0, r0=9.12);
#harmonic.bond_coeff.set('bc', k=100.0, r0=1);
#harmonic.bond_coeff.set('cd', k=100.0, r0=9.12);
#harmonic.bond_coeff.set('da', k=100.0, r0=4);

harmonic.bond_coeff.set('ab', k=100.0, r0=1);
harmonic.bond_coeff.set('bc', k=100.0, r0=1);
harmonic.bond_coeff.set('cd', k=100.0, r0=1);
harmonic.bond_coeff.set('da', k=100.0, r0=1);


hoomd.md.integrate.mode_standard(dt=0.001);
group_all = hoomd.group.all();
hoomd.md.integrate.nve(group=group_all);

hoomd.analyze.log(filename="log-output.log",
                  quantities=['potential_energy', 'temperature'],
                  period=500,
                  overwrite=True);

hoomd.dump.gsd("trajectory.gsd", period=20, group=group_all, overwrite=True);

hoomd.run(5e4);






