BUILD_DIR = build
CASCADELAKE_BUILD_DIR = build_CascadeLake
ICELAKE_BUILD_DIR = build_IceLake
DEBUG_BUILD_DIR = build_debug

RESIDUALS ?= OFF
FULL_DIAG ?= OFF

# Default target to build the project
all: $(BUILD_DIR)/Makefile
	@$(MAKE) -C $(BUILD_DIR)

$(BUILD_DIR)/Makefile: CMakeLists.txt
	@mkdir -p $(BUILD_DIR)
	@cd $(BUILD_DIR) && cmake \
		-DCMAKE_CXX_COMPILER=g++ \
		-DLATTICE_CUT_RESIDUALS=$(RESIDUALS) \
		-DLATTICE_CUT_FULL_DIAG=$(FULL_DIAG) \
		..

cascadelake: $(CASCADELAKE_BUILD_DIR)/Makefile
	@$(MAKE) -C $(CASCADELAKE_BUILD_DIR)

$(CASCADELAKE_BUILD_DIR)/Makefile: CMakeLists.txt
	@mkdir -p $(CASCADELAKE_BUILD_DIR)
	@cd $(CASCADELAKE_BUILD_DIR) && cmake \
	    -DCMAKE_CXX_COMPILER=g++ \
	    -DCLUSTER_BUILD=cascadelake \
	    -DLATTICE_CUT_RESIDUALS=$(RESIDUALS) \
	    -DLATTICE_CUT_FULL_DIAG=$(FULL_DIAG) \
	    ..

icelake: $(ICELAKE_BUILD_DIR)/Makefile
	@$(MAKE) -C $(ICELAKE_BUILD_DIR)

$(ICELAKE_BUILD_DIR)/Makefile: CMakeLists.txt
	@mkdir -p $(ICELAKE_BUILD_DIR)
	@cd $(ICELAKE_BUILD_DIR) && cmake \
	    -DCMAKE_CXX_COMPILER=g++ \
	    -DCLUSTER_BUILD=icelake \
	    -DLATTICE_CUT_RESIDUALS=$(RESIDUALS) \
	    -DLATTICE_CUT_FULL_DIAG=$(FULL_DIAG) \
	    ..

debug: $(DEBUG_BUILD_DIR)/Makefile
	@$(MAKE) -C $(DEBUG_BUILD_DIR)

$(DEBUG_BUILD_DIR)/Makefile: CMakeLists.txt
	@mkdir -p $(DEBUG_BUILD_DIR)
	@cd $(DEBUG_BUILD_DIR) && cmake \
		-DCMAKE_CXX_COMPILER=g++ \
		-DCMAKE_BUILD_TYPE=Debug \
		-DLATTICE_CUT_RESIDUALS=$(RESIDUALS) \
		-DLATTICE_CUT_FULL_DIAG=$(FULL_DIAG) \
		..

clean:
	@rm -rf $(BUILD_DIR) $(CASCADELAKE_BUILD_DIR) $(ICELAKE_BUILD_DIR) $(DEBUG_BUILD_DIR) build_header
	@rm -rf $(BUILD_DIR) $(CASCADELAKE_BUILD_DIR)_ed $(ICELAKE_BUILD_DIR)_ed $(DEBUG_BUILD_DIR)_ed
	@rm -rf auto_generated*

.PHONY: all clean icelake cascadelake debug