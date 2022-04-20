#!/usr/bin/env python3

import click
from subprocess import run
import os

@click.command()
@click.option('--iter', '-t', default=1, help='Number of iterations', type=int)
@click.option('--input', '-i', default='test', help='Name of the input file')
@click.option('--output', '-o', default='test', help='Name of the output file')
def main(iter, input, output):
    # input definition
    com_path = 'build'
    command = 'SpatialResol.exe'
    flags = ['-b', '-r']

    project_path = os.path.abspath(os.path.dirname(
                                   os.path.abspath(__file__)
                                   ) + '/../' + com_path)

    comm = [project_path + '/' + command]
    comm.extend(['-i', input])
    # comm.extend(['--param', project_path + '/../params/diag.ini'])

    for flag in flags:
        comm.append(flag)
    print(iter)
    for iterId in range(iter):
        exe = comm.copy()
        exe.extend(['-t', str(iterId)])
        exe.extend(['-o', f'{str(output)}_iter{iterId}.root'])
        print(exe)
        run(exe, check=True)

if __name__ == "__main__":
    main()