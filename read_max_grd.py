import pygmt
import numpy as np

# Directory containing the GRD files
directory = "2.5deg_grad_grds/"

azimuths = range(-170, 190, 10)

region = [-158, -142, 58, 69]

vmins = []
vmaxs = []

# Loop over azimuth values and extract min and max from each GRD file
for azimuth in azimuths:
    grdfile = f"{directory}grad_mask.{azimuth}.grd"
    info = pygmt.grdinfo(grdfile, per_column=True,R=region).strip().split()
    vmin = float(info[4])
    vmax = float(info[5])
    vmins.append(vmin)
    vmaxs.append(vmax)

# Calculate average min and max
avg_vmin = np.mean(vmins)
avg_vmax = np.mean(vmaxs)

# Output the results
print("Minimum values:", vmins,'\n')
print("Min of all values:", np.min(vmins))
print('#######')
print("Maximum values:", vmaxs)
print("Max of all values:", np.max(vmaxs))

print("Average minimum:", avg_vmin)
print("Average maximum:", avg_vmax)
