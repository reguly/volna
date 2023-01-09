#!/bin/bash
if [[ -v OP2_OLD ]]; then
       python3 "${OP2_INSTALL_PATH}"/../translator/c/python/op2.py volna.cpp volna_event.cpp volna_init.cpp volna_output.cpp volna_simulation.cpp
else
       python3 "${OP2_INSTALL_PATH}"/../translator/c/op2.py volna.cpp volna_event.cpp volna_init.cpp volna_output.cpp volna_simulation.cpp
fi
sed -i 's/NumericalFluxes_veckernel.cpp/..\/seq\/NumericalFluxes_seqkernel.cpp/g' vec/volna_veckernels.cpp
sed -i 's/SIMD_VEC 4/SIMD_VEC 16/g' vec/volna_veckernels.cpp
sed -i 's/SIMD_VEC 8/SIMD_VEC 16/g' sycl/volna_kernels.cpp
