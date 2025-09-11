Performance measurements
========================

To keep track of the run time performance of the most important tools of the toolset an automated process has been setup to perform benchmarks. These benchmarks also periodically performed on our build servers to measure how the performance changes over builds. If you have an account for our TeamCity servers you can view the results of these measurements `here <http://mcrl2build1.win.tue.nl:8111/viewType.html?buildTypeId=mcrl2_Nightly_Default_UbuntuAmd64_Benchmark>`_.

The details of the benchmarks are described below. To perform these benchmarks yourself, see the Setup section.

Benchmarks
-----------

The run time of the following tools are currently being measured:

  ========= ============= ============
  Tool	       Case 1        Case 2
  ========= ============= ============
  mcrl22lps ``[default]`` N/A
  lps2lts   ``[default]`` ``-rjittyc``
  pbessolve ``[default]`` ``-rjittyc``
  ========= ============= ============

Each tool is used on a number of mCRL2 examples. The process performed depends on the tool, i.e, linearisation, state space exploration and model checking of the `nodeadlock` property. The following examples present in our repository are used to perform the benchmarks on:

- academic/abp/abp.mcrl2
- academic/allow/allow.mcrl2
- academic/cabp/cabp.mcrl2
- academic/dining/dining8.mcrl2
- games/clobber/clobber.mcrl2
- games/domineering/domineering.mcrl2
- games/magic_square/magic_square.mcrl2
- games/othello/othello.mcrl2
- industrial/1394/1394-fin.mcrl2
- industrial/brp/brp.mcrl2
- industrial/chatbox/chatbox.mcrl2
- industrial/lift/lift3-final.mcrl2
- industrial/lift/lift3-init.mcrl2
 
Besides tool benchmarks there are also microbenchmarks that measure specific aspects of individual libraries.

Setup
-----------

For compiling the toolset itself and the build configuration, see :ref:`build-instructions`. If you are able to compile the toolset then enable the `MCRL2_ENABLE_BENCHMARKS`` CMake option to enable the generation of benchmark targets. Compile the ``benchmarks`` target, this will generate the necessary files to perform the benchmarks. It will also compile the tools and microbenchmarks if necessary.

After this step, performing the benchmarks can be done using ctest. For each tool, specification and option it generates a target named ``benchmark_<tool>_<specification>_<options>``. Each benchmark target has the ``benchmark`` label. For example to run all benchmarks pertaining the 1394-fin specification use ``ctest -L benchmark -R 1394-fin``, or to benchmark a specific tool use ``ctest -L benchmark -R lps2lts``. 