# Author(s): Maurice Laveaux
# Copyright: see the accompanying file COPYING or copy at
# https://github.com/mCRL2org/mCRL2/blob/master/COPYING
#
# Distributed under the Boost Software License, Version 1.0.
# (See accompanying file LICENSE_1_0.txt or copy at
# http://www.boost.org/LICENSE_1_0.txt)
#

import os
import unittest

import LabelledTransitionSystem, LinearProcess from mcrl2py

dir_path = os.path.dirname(os.path.realpath(__file__))

class TestModule(unittest.TestCase):

    def test_exploration(self):
        with open(os.path.join(dir_path, "../../examples/academic/abp/abp.mcrl2"), "r", encoding="utf-8") as f:
            options = LinearisationOptions()
            options.method = Linearisation.stack

            # Read the mCRL2 specification and apply the linearisation options
            spec = LinearProcess(f.read(), options)

            with open(os.path.join(dir_path, "abp.lts"), "wb") as f:
                # Explore the state space and write it to disk
                lts = LabelledTransitionSystem(spec)
            

if __name__ == '__main__':
    unittest.main()
