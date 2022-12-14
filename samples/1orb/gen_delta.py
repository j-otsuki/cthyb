import numpy as np

beta = 100.0
nw = 1024
vsq = 0.1

print(f"beta = {beta}")
print(f"nw = {nw}")
print(f"V^2 = {vsq}")

w = (2*np.arange(nw)+1) * np.pi / beta
delta_im = -np.arctan(1/w) * vsq

with open("delta_w.in", "w") as f:
    for d in delta_im:
        print(f"0.0 {d:.8e} 0.0 {d:.8e}", file=f)
print("'delta_w.in' generated")
