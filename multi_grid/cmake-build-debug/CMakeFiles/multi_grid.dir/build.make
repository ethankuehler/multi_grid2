# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.15

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /home/dcrush/.local/share/JetBrains/Toolbox/apps/CLion/ch-0/192.7142.39/bin/cmake/linux/bin/cmake

# The command to remove a file.
RM = /home/dcrush/.local/share/JetBrains/Toolbox/apps/CLion/ch-0/192.7142.39/bin/cmake/linux/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/dcrush/CLionProjects/multi_grid/multi_grid

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/dcrush/CLionProjects/multi_grid/multi_grid/cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/multi_grid.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/multi_grid.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/multi_grid.dir/flags.make

CMakeFiles/multi_grid.dir/multi_grid.c.o: CMakeFiles/multi_grid.dir/flags.make
CMakeFiles/multi_grid.dir/multi_grid.c.o: ../multi_grid.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/dcrush/CLionProjects/multi_grid/multi_grid/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object CMakeFiles/multi_grid.dir/multi_grid.c.o"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/multi_grid.dir/multi_grid.c.o   -c /home/dcrush/CLionProjects/multi_grid/multi_grid/multi_grid.c

CMakeFiles/multi_grid.dir/multi_grid.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/multi_grid.dir/multi_grid.c.i"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/dcrush/CLionProjects/multi_grid/multi_grid/multi_grid.c > CMakeFiles/multi_grid.dir/multi_grid.c.i

CMakeFiles/multi_grid.dir/multi_grid.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/multi_grid.dir/multi_grid.c.s"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/dcrush/CLionProjects/multi_grid/multi_grid/multi_grid.c -o CMakeFiles/multi_grid.dir/multi_grid.c.s

CMakeFiles/multi_grid.dir/operators.c.o: CMakeFiles/multi_grid.dir/flags.make
CMakeFiles/multi_grid.dir/operators.c.o: ../operators.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/dcrush/CLionProjects/multi_grid/multi_grid/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building C object CMakeFiles/multi_grid.dir/operators.c.o"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/multi_grid.dir/operators.c.o   -c /home/dcrush/CLionProjects/multi_grid/multi_grid/operators.c

CMakeFiles/multi_grid.dir/operators.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/multi_grid.dir/operators.c.i"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/dcrush/CLionProjects/multi_grid/multi_grid/operators.c > CMakeFiles/multi_grid.dir/operators.c.i

CMakeFiles/multi_grid.dir/operators.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/multi_grid.dir/operators.c.s"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/dcrush/CLionProjects/multi_grid/multi_grid/operators.c -o CMakeFiles/multi_grid.dir/operators.c.s

# Object files for target multi_grid
multi_grid_OBJECTS = \
"CMakeFiles/multi_grid.dir/multi_grid.c.o" \
"CMakeFiles/multi_grid.dir/operators.c.o"

# External object files for target multi_grid
multi_grid_EXTERNAL_OBJECTS =

libmulti_grid.a: CMakeFiles/multi_grid.dir/multi_grid.c.o
libmulti_grid.a: CMakeFiles/multi_grid.dir/operators.c.o
libmulti_grid.a: CMakeFiles/multi_grid.dir/build.make
libmulti_grid.a: CMakeFiles/multi_grid.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/dcrush/CLionProjects/multi_grid/multi_grid/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking C static library libmulti_grid.a"
	$(CMAKE_COMMAND) -P CMakeFiles/multi_grid.dir/cmake_clean_target.cmake
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/multi_grid.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/multi_grid.dir/build: libmulti_grid.a

.PHONY : CMakeFiles/multi_grid.dir/build

CMakeFiles/multi_grid.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/multi_grid.dir/cmake_clean.cmake
.PHONY : CMakeFiles/multi_grid.dir/clean

CMakeFiles/multi_grid.dir/depend:
	cd /home/dcrush/CLionProjects/multi_grid/multi_grid/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/dcrush/CLionProjects/multi_grid/multi_grid /home/dcrush/CLionProjects/multi_grid/multi_grid /home/dcrush/CLionProjects/multi_grid/multi_grid/cmake-build-debug /home/dcrush/CLionProjects/multi_grid/multi_grid/cmake-build-debug /home/dcrush/CLionProjects/multi_grid/multi_grid/cmake-build-debug/CMakeFiles/multi_grid.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/multi_grid.dir/depend

