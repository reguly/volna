# Modern OP2 build for Volna.
#
# The current OP2 translator consumes the original application sources and
# writes backend-specific products under generated/.  Do not check generated
# sources into this tree.

include $(OP2_INSTALL_PATH)/../makefiles/common.mk

APP_NAME := volna
APP_SRC := \
  volna.cpp \
  volna_event.cpp \
  volna_init.cpp \
  volna_output.cpp \
  volna_simulation.cpp \
  volna_util.cpp \
  volna_writeVTK.cpp

# Volna uses the HDF5 C API directly, so expose its headers to both the
# translator's front end and normal sequential compilation.
APP_INC := -I. $(HDF5_SEQ_INC)
OP2_LIBS_WITH_HDF5 := true

# Besides the core HDF5 API used by OP2, Volna calls the HDF5 high-level
# routines (H5LT*).  Keep these flags local to this application.
override HDF5_SEQ_LIB += -lhdf5_hl
override HDF5_PAR_LIB += -lhdf5_hl

include $(OP2_INSTALL_PATH)/../makefiles/c_app.mk
