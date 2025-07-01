import MDAnalysis as mda
from MDAnalysis.analysis import dihedrals
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# Load your trajectory and topology
u = mda.Universe('topology', 'trajectory')

#Define the torsion indices (this is just an example, adjust according to your system)
# Ensure that each torsion contains exactly 4 atoms
torsion_indices = [
    [2813, 2814, 2829, 2830],  
    [2812, 2813, 2815, 2816]
# Add more torsions as needed, each with exactly 4 atoms
]

atomgroups = []
for indices in torsion_indices:
    try:
        ag = u.atoms[indices]
        if len(set(indices)) == 4:
            atomgroups.append(ag)
        else:
            print(f"Skipping torsion with duplicate atoms: {indices}")
    except IndexError:
        print(f"Skipping invalid torsion {indices}: index out of bounds")

if not atomgroups:
    raise ValueError("No valid torsions found for analysis.")

print(f"Number of valid torsion groups: {len(atomgroups)}")

dihedral_analysis = dihedrals.Dihedral(atomgroups).run()

angles_deg = dihedral_analysis.angles 

df = pd.DataFrame(angles_deg, columns=[f'Torsion_{i+1}' for i in range(angles_deg.shape[1])])
df.to_csv('torsion_angles.dat', sep='\t', index=False)
print("Torsion angles saved to torsion_angles.dat")

plt.figure(figsize=(12, 6))
for i in range(angles_deg.shape[1]):
    plt.plot(angles_deg[:, i], label=f'Torsion {i+1}')
plt.xlabel("Frame")
plt.ylabel("Torsion Angle (degrees)")
plt.title("Torsion Angles Over Time")
plt.legend()
plt.tight_layout()
plt.savefig("torsion_angles_plot.png")
plt.show()

