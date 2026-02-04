import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', required=True, help='Input LAMMPS Data File')
# parser.add_argument('-o', '--output', required=True, help='Output types.txt file')
args = parser.parse_args()

nparticle_types = 0
particle_types = []
nbond_types = 0
bond_types = []
nangle_types = 0
angle_types = []
ndihedral_types = 0
dihedral_types = []
nimproper_types = 0
improper_types = []

with open(args.input, 'r') as f:
    for line in f:
        line = line.split()
        if len(line) > 1 and line[1] == 'atom':
            nparticle_types = int(line[0])
        elif len(line) > 1 and line[1] == 'bond':
            nbond_types = int(line[0])
        elif len(line) > 1 and line[1] == 'angle':
            nangle_types = int(line[0])
        elif len(line) > 1 and line[1] == 'dihedral':
            ndihedral_types = int(line[0])
        elif len(line) > 1 and line[1] == 'improper':
            nimproper_types = int(line[0])
        elif len(line) > 0 and line[0] == 'Masses':
            next(f)
            for i in range(nparticle_types):
                line = next(f)
                particle_types.append(line.split('#')[1].strip())
        elif len(line) > 0 and line[0] == 'Bond':
            next(f)
            for i in range(nbond_types):
                line = next(f)
                bond_types.append(line.split('#')[1].strip())
        elif len(line) > 0 and line[0] == 'Angle':
            next(f)
            for i in range(nangle_types):
                line = next(f)
                angle_types.append(line.split('#')[1].strip())
        elif len(line) > 0 and line[0] == 'Dihedral':
            next(f)
            for i in range(ndihedral_types):
                line = next(f)
                dihedral_types.append(line.split('#')[1].strip())
        elif len(line) > 0 and line[0] == 'Improper':
            next(f)
            for i in range(nimproper_types):
                line = next(f)
                improper_types.append(line.split('#')[1].strip())


with open('types.txt', 'w') as fout:
    if nparticle_types > 0:
        fout.write('atom types\n')
        for i, pt in enumerate(particle_types):
            fout.write(' %4d \t%s \n' % (i+1, pt))
        fout.write('#\n')
    if nbond_types > 0:
        fout.write('bond types\n')
        for i, bt in enumerate(bond_types):
            fout.write(' %4d \t%s \n' % (i+1, bt))
        fout.write('#\n')
    if nangle_types > 0:
        fout.write('angle types\n')
        for i, at in enumerate(angle_types):
            fout.write(' %4d \t%s \n' % (i+1, at))
        fout.write('#\n')
    if ndihedral_types > 0:
        fout.write('dihedral types\n')
        for i, dt in enumerate(dihedral_types):
            fout.write(' %4d \t%s \n' % (i+1, dt))
        fout.write('#\n')
    if nimproper_types > 0:
        fout.write('improper types\n')
        for i, it in enumerate(improper_types):
            fout.write(' %4d \t%s \n' % (i+1, it))
        fout.write('#\n')

exit()

