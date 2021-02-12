import subprocess
import matplotlib.pyplot as plt
import numpy as np
import glob
import os
import random

algos = ['lagrange_naive', 'lagrange_fft', 'lagrange_div', 'vandermonde', 'fast']

print('rep: ')
rep = int(input())
print('n_max: ')
n_max = int(input())


def run(exe: str):
    print('Running', exe + '...')
    ret = []

    for n in range(1, n_max + 1):
        avg = .0
        for __ in range(rep):
            # Coefficients
            a = [random.uniform(-1, 1) for i in range(n)]

            # Evaluation
            def val(z):
                v = 0
                p = 1
                for c in a:
                    v += c * p
                    p *= z
                return v

            # Error
            def delta(_a):
                d = .0
                for e, _e, in zip(a, _a):
                    d += abs(_e - e) / max(1, e)
                return d

            # Sampling point
            x = [random.uniform(-1, 1) for i in range(n)]
            y = [val(z) for z in x]

            txt = str(n) + '\n' + ' '.join([str(s) for s in x]) + '\n' + ' '.join([str(s) for s in y]) + '\n'
            proc = subprocess.run('./' + exe, input=txt.encode(), stdout=subprocess.PIPE)
            _a = list(map(float, proc.stdout.split()))
            avg += delta(_a)
            # print(a, _a)

        avg /= rep
        avg /= n
        ret.append(avg)
    print('Done!')
    return ret

# Setup plot
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

ax.set_yscale('log')
ax.set_title('Precision of interpolation')
ax.set_xlabel('Degree')
ax.set_ylabel('Average of delta')

ax.axhline(y=1, color='gray', linestyle='dashed')

xtks=[0]
step=3
while xtks[-1]+step < n_max:
    xtks.append(xtks[-1]+step)
ax.set_xticks(xtks)


# Run each algorithm
for a in algos:
    with open(a + '.cpp', 'w') as src:
        src.write('#include "../src/all.h"\n')
        src.write('int main() { print(' + a + '(read())); }\n')
    print('Compiling', a + '.cpp ...')
    subprocess.call(["g++", a + ".cpp", "-std=c++17", "-o", a])
    ax.plot([i for i in range(n_max)], run(a), label=a)

ax.legend()

print('Save as:')
fig.savefig(input())
plt.close(fig)
