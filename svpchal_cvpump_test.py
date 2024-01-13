#!/usr/bin/env python
# -*- coding: utf-8 -*-
####
#
#   Copyright (C) 2018-2021 Team G6K
#
#   This file is part of G6K. G6K is free software:
#   you can redistribute it and/or modify it under the terms of the
#   GNU General Public License as published by the Free Software Foundation,
#   either version 2 of the License, or (at your option) any later version.
#
#   G6K is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with G6K. If not, see <http://www.gnu.org/licenses/>.
#
####


"""
LWE Challenge Solving Command Line Client
"""

from __future__ import absolute_import
from __future__ import print_function
import copy
import re
import sys
import time

from collections import OrderedDict # noqa
from math import log

from fpylll import BKZ as fplll_bkz
from fpylll.algorithms.bkz2 import BKZReduction
from fpylll.tools.quality import basis_quality
from fpylll.util import gaussian_heuristic

from g6k.algorithms.bkz import pump_n_jump_bkz_tour
from g6k.algorithms.cvpump import cvpump
from g6k.algorithms.pump import pump
from g6k.siever import Siever
from g6k.utils.cli import parse_args, run_all, pop_prefixed_params
from g6k.utils.stats import SieveTreeTracer, dummy_tracer
from g6k.utils.util import load_lwe_challenge

from g6k.utils.lwe_estimation import gsa_params, lattice_basis
from six.moves import range
from copy import deepcopy
from fpylll import IntegerMatrix
from g6k.utils.util import load_svpchallenge_and_randomize
import numpy as np
import os
from g6k.algorithms.cworkout import cworkout
from g6k.algorithms.workout import workout
from g6k.utils.util import sanitize_params_names
import six
from g6k.algorithms.pumpcvp import pumpcvp

def load_svp_midmat(d):
    """
    Load svp challenge midmat from file.

    :param d: svp dimension.
    
    """

    start = "svpchallenge/"
    if not os.path.isdir(start):
        os.mkdir(start)
    
    end = "{d:03d}-midmat.txt".format(d=d)
    filename = os.path.join(start, end)
    try:
        data = open(filename, "r").readlines()
    except FileNotFoundError:
        return None
    B = eval(",".join([s_.replace(" ", ", ") for s_ in data[0:]]))
    B = IntegerMatrix.from_matrix(B)
 
    
    return B



def asvp_kernel(arg0, params=None, seed=None):
    """
    Run SVP instance.

    :param n: the dimension of the svp-challenge secret
    - verbose: print information throughout the lwe challenge attempt

    """

    # Pool.map only supports a single parameter
    if params is None and seed is None:
        n, params, seed = arg0
    else:
        n = arg0

    params = copy.copy(params)

    # misc
    verbose = params.pop("verbose")
    workout_params = pop_prefixed_params("workout", params)


    A = load_svp_midmat(n)
    if(A is None):
        A, _ = load_svpchallenge_and_randomize(n, s=0, seed=0)
    if verbose:
        print(("Loaded challenge dim %d" % n))

    
    g6k = Siever(A, params)
    gh = gaussian_heuristic([g6k.M.get_r(i, i) for i in range(n)])
    goal_r0 = (1.05 ** 2) * gh

    print("gh = %d, goal_r0 = %.3f" %(gh,goal_r0))
    
    
    print("GSO precision: ", g6k.M.float_type)

    g6k.lll(0, g6k.full_n)
    slope = basis_quality(g6k.M)["/"]
    print("Intial Slope = %.5f\n" % slope)
    
    
    # pump_params = {"down_sieve": True}
    # workout(
    #     g6k, dummy_tracer, 0, n, dim4free_min= 40, pump_params=pump_params, **workout_params, verbose=True)
    
    
    T0 = time.time()
    # pumpcvp(g6k,dummy_tracer, 0, n, 30, verbose=True,  goal_r0 = goal_r0, down_sieve = True, min_cvp_dim = 60)
    
    cworkout(
        g6k, dummy_tracer, 0, n, 70, goal_r0=goal_r0, pump_params={"down_sieve":True}, verbose=True, start_n= 73,**workout_params
    )
    
    g6k.lll(0,g6k.full_n)
    print("wall time: %.2f sec" %(time.time()-T0))
    
    ee = g6k.M.B[0]
    
    print(ee)
    
    norm_ee = sum([_**2 for _ in ee]) 
    print("norm(ee) = %d" %norm_ee)
    if  ee is not None and norm_ee < goal_r0:
        print("Finished!")
        print(ee)
        filename = 'svpchallenge/%03d-solution.txt' % (n)
        fn = open(filename, "w")
        fn.write(str(ee))
        fn.close()
        return

    raise ValueError("No solution found.")


def asvp():
    """
    Run a Workout until 1.05-approx-SVP on matrices with dimensions in ``range(lower_bound, upper_bound, step_size)``.
    """
    
    description = asvp.__doc__

    args, all_params = parse_args(
        description,
        # load_matrix=None,
        verbose=True,
        # challenge_seed=0,
        workout__dim4free_dec=3,
    )

    stats = run_all(
        asvp_kernel,
        list(all_params.values()),
        lower_bound=args.lower_bound,
        upper_bound=args.upper_bound,
        step_size=args.step_size,
        trials=args.trials,
        workers=args.workers,
        seed=args.seed,
    )

    inverse_all_params = OrderedDict([(v, k) for (k, v) in six.iteritems(all_params)])
    stats = sanitize_params_names(stats, inverse_all_params)
    

if __name__ == "__main__":
    asvp()
