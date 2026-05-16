import glob

# Get all .vtu files sorted by step number
files = sorted(glob.glob("solution/cahnhilliard_step*.vtu"))

with open("cahnhilliard.pvd", "w") as f:
    f.write('<?xml version="1.0"?>\n')
    f.write('<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">\n')
    f.write('  <Collection>\n')

    for i, filename in enumerate(files):
        f.write(f'    <DataSet timestep="{i}" group="" part="0" file="{filename}"/>\n')

    f.write('  </Collection>\n')
    f.write('</VTKFile>\n')

print(f"✅ Created cahnhilliard.pvd with {len(files)} files.")
