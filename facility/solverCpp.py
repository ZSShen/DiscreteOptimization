#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
from subprocess import Popen, PIPE

def solve_it(input_data, engine, time_limit):

    # Writes the inputData to a temporay file

    tmp_file_name = 'tmp.data'
    tmp_file = open(tmp_file_name, 'w')
    tmp_file.write(input_data)
    tmp_file.close()

    # Runs the command: java Solver -file=tmp.data

    # Please specify the environment variable DYLD_LIBRARY_PATH to indiciate
    # the path of google-or-tool library.
    env = os.environ
    #env["DYLD_LIBRARY_PATH"] = "PATH_TO_GOOGLE_OR_TOOLS"

    process = Popen(['./%s' % engine, tmp_file_name, time_limit], stdout=PIPE, env=env)
    (stdout, stderr) = process.communicate()

    # removes the temporay file
    os.remove(tmp_file_name)

    return stdout.strip()


import sys

if __name__ == '__main__':
    if len(sys.argv) > 1:
        file_location = sys.argv[1].strip()
        engine = sys.argv[2]
        time_limit = sys.argv[3].strip()

        with open(file_location, 'r') as input_data_file:
            input_data = input_data_file.read()

        print solve_it(input_data, engine, time_limit)
    else:
        print('This test requires an input file.  Please select one from the data directory. (i.e. python solver.py ./data/ks_4_0)')