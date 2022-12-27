#!/bin/bash
python3 ~/OP2-Common-sycl/translator/c/python/op2.py volna.cpp volna_event.cpp volna_init.cpp volna_output.cpp volna_simulation.cpp
sed -i 's/SIMD_VEC 8/SIMD_VEC 16/g' sycl/volna_kernels.cpp
