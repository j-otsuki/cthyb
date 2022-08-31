import numpy as np

beta = 10.0
nw = 512
vsq = 0.1
norb = 3

print(f"beta = {beta}")
print(f"nw = {nw}")
print(f"V^2 = {vsq}")
print(f"norb = {norb}")

w = (2*np.arange(nw)+1) * np.pi / beta
delta_im = -np.arctan(1/w) * vsq

with open("delta_w.in", "w") as f:
    for d in delta_im:
        for _ in range(norb*2):
            print(f" 0.0 {d:.8e}", file=f, end="")
        print("", file=f)
print("'delta_w.in' generated")
