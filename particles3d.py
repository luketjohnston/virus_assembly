import hoomd
import hoomd.md
import numpy
import random

hoomd.context.initialize("");


# will start with 4 particles
# box is 2d 
# two types of particles

# 2d box, 10x10
numParticles = 10
numAtomsPerParticle = 3
length = 40
box = hoomd.data.boxdim(L=length,dimensions=3) # dimensions=3 forces Lz=1
snapshot = hoomd.data.make_snapshot(N=numParticles*numAtomsPerParticle,
                                    box=box,
                                    particle_types=['A', 'B'],
                                    bond_types=['polymer']);


# put particles at s
print(snapshot.particles.position.shape)
# have 3d coordinates because the way hoomd handles 2d simulations is
# by just operating in 3d with a small z-axis (Lz of the above
# box is equal to 1)
initialPositions = []
for i in range(numParticles):
    curPosition = []
    for j in range(3):
        curPosition.append(random.random()*length / 2.0)
    initialPositions.append(curPosition)
fullPositions = []
for position in initialPositions:
    fullPositions.append(position)
    fullPositions.append([position[0]-1,position[1],position[2]])
    fullPositions.append([position[0],position[1]-1,position[2]])


snapshot.particles.position[:] = fullPositions
for i in range(numParticles*numAtomsPerParticle):
    if i % 3 == 0:
        snapshot.particles.typeid[i] = 0 
    else:
        snapshot.particles.typeid[i] = 1

snapshot.bonds.resize(numParticles*2)
bonds = []
for i in range(numParticles):
    bonds.append([i*3,i*3+1])
    bonds.append([i*3,i*3+2])
snapshot.bonds.group[:] = bonds


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

# use gaussian interaction because it is mathematicaaly simple
gauss = hoomd.md.pair.gauss(r_cut=10, nlist=nl)

# Points want to interact
gauss.pair_coeff.set('A', 'A',epsilon=400.0, sigma=1.0);
#Opposites want to repel slightlys
gauss.pair_coeff.set('B', 'B',epsilon=-40.0, sigma=1.0);
# opposites have weak repel
gauss.pair_coeff.set('A', 'B',epsilon=-40.0, sigma=1.0);

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

hoomd.run(10000);
