from branches.roots.athena_analysis import *

filename = file.data_loc + "Cyl_1/disk.out1.01000.athdf"

aa = Athena_Analysis(filename, grid_type = "Cylindrical")

aa.get_primaries(get_rho=True)
aa.native_grid(get_r=True, get_phi=True, get_z=True, get_dr=True, get_dphi=True, get_dz=True)

now = datetime.now()

rust_intV = aa.integrate(aa.rho, "All", RUSTY=True)

rust_timeV = datetime.now() - now

now = datetime.now()

python_intV = aa.integrate(aa.rho, "All", RUSTY=False)

python_timeV = datetime.now() - now

print("Python Time: ")
print(python_timeV)
print("Python Value: ")
print(python_intV)
print("Rust Time: ")
print(rust_timeV)
print("Rust Value: ")
print(rust_intV)

now = datetime.now()

rust_intA = aa.integrate(aa.rho, "Shell", RUSTY=True)

rust_timeA = datetime.now() - now

now = datetime.now()

python_intA = aa.integrate(aa.rho, "Shell", RUSTY=False)

python_timeA = datetime.now() - now

vert = 1
horz = 1
gs = gridspec.GridSpec(vert, horz)
fig = plt.figure(figsize=(horz*3, vert*3), dpi=300)

ax = fig.add_subplot(gs[0, 0])

ax.plot(aa.possible_r, python_intA, label="PyTime: " + str(python_timeA))
ax.plot(aa.possible_r, rust_intA, label="RustTime: " + str(rust_timeA), linestyle="--")

ax.legend()

plt.savefig("/home/morgan/Athena_Analysis/integration_bowl.png")

print(aa.integrate(aa.rho, "z", RUSTY=True)[0,0])
print(aa.integrate(aa.rho, "z", RUSTY=False)[0,0])
