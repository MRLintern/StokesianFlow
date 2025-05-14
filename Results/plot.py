import matplotlib.pyplot as plt

# import calculated velocity data
with open("uSolution.dat", "r") as f:
    lines = [line.strip() for line in f if line.strip()]

# Parse data (every 2 lines: position, then u-velocity)
x_vals, y_vals, u_vals, v_vals = [], [], [], []

i = 0
while i < len(lines) - 1:
    x, y = map(float, lines[i].split())

    # fluid velocity components
    u = float(lines[i + 1])
    v = 0.0

    # node coordinates  
    x_vals.append(x)
    y_vals.append(y)

    # velocity plotting components
    u_vals.append(u)
    v_vals.append(v)

    i += 2  # Move to the next node block

# quiver plot
plt.figure(figsize=(8, 6))
plt.quiver(x_vals, y_vals, u_vals, v_vals, angles='xy', scale_units='xy', scale=1, color='blue')
plt.xlabel('Node Position')
plt.ylabel('Velocity (m/s)')
plt.title('Velocity of Fluid Flowing through the Channel/River')
plt.grid(True)
plt.axis('equal')
plt.tight_layout()
plt.savefig('Velocity_Profile.png', dpi=300) # can also use: .jpg, pdf ...
plt.show()
