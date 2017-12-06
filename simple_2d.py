import hoomd
import hoomd.md
import numpy


hoomd.context.initialize("");


# will start with 4 particles
# box is 2d 
# two types of particles

# 2d box, 10x10
box = hoomd.data.boxdim(L=10,dimensions=2) # dimensions=2 forces Lz=1
snapshot = hoomd.data.make_snapshot(N=4,
                                    box=box,
                                    particle_types=['A', 'B'],
                                    bond_types=['polymer']);


# put particles at s
print(snapshot.particles.position.shape)
# have 3d coordinates because the way hoomd handles 2d simulations is
# by just operating in 3d with a small z-axis (Lz of the above
# box is equal to 1)
snapshot.particles.position[:] = [[0,0,0],[0,1,0],[5,0,0],[5,1,0]]

# particles of type 'A'
snapshot.particles.typeid[0] = 0 
snapshot.particles.typeid[2] = 0

# particles of type 'B'
snapshot.particles.typeid[1] = 1 
snapshot.particles.typeid[3] = 1

# we have two bonds
snapshot.bonds.resize(2)
snapshot.bonds.group[:] = [[0,1],[2,3]]


# set gaussian random velocity for all particles
snapshot.particles.velocity[:] = numpy.random.normal(0.0,
  numpy.sqrt(0.8 / 1.0), [snapshot.particles.N, 3]);

# initialize hoomd with this snapshot
hoomd.init.read_snapshot(snapshot)

# this tells hoomd what method to use to find the particles that interact
# with eachother. Each interaction type will have a cutoff length, and particles
# farther than that cutoff length will not interact. the "cell" method 
# apparently performs well when there is only one interaction type.
nl = hoomd.md.nlist.cell();

# for now just use this random interaction type
dpd = hoomd.md.pair.dpd(r_cut=1.0, nlist=nl, kT=0.8, seed=1);

dpd.pair_coeff.set('A', 'A', A=25.0, gamma = 1.0);
dpd.pair_coeff.set('A', 'B', A=100.0, gamma = 1.0);
dpd.pair_coeff.set('B', 'B', A=25.0, gamma = 1.0);
# dpd doesn't exclude bonded atoms, so reset exclusions to nothing (before this line they are the bonded atoms).
nl.reset_exclusions(exclusions = []);
# "apply harmonic bonds between the directly bonded particles"
harmonic = hoomd.md.bond.harmonic();
# dpd uses bond length zero to keep particles together (dpd pair repulsion will push them apart)
harmonic.bond_coeff.set('polymer', k=100.0, r0=0);

hoomd.md.integrate.mode_standard(dt=0.01);

group_all = hoomd.group.all();
hoomd.md.integrate.nve(group=group_all);

hoomd.analyze.log(filename="log-output.log",
                  quantities=['potential_energy', 'temperature'],
                  period=500,
                  overwrite=True);

hoomd.dump.gsd("trajectory.gsd", period=10e3, group=group_all, overwrite=True);

hoomd.run(5e4);






