# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.14

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
CMAKE_COMMAND = /Applications/CLion.app/Contents/bin/cmake/mac/bin/cmake

# The command to remove a file.
RM = /Applications/CLion.app/Contents/bin/cmake/mac/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/guilhermemorais/CLionProjects/icc/replacement

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/guilhermemorais/CLionProjects/icc/replacement/cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/gauss_elimination.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/gauss_elimination.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/gauss_elimination.dir/flags.make

CMakeFiles/gauss_elimination.dir/main.c.o: CMakeFiles/gauss_elimination.dir/flags.make
CMakeFiles/gauss_elimination.dir/main.c.o: ../main.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/guilhermemorais/CLionProjects/icc/replacement/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object CMakeFiles/gauss_elimination.dir/main.c.o"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/gauss_elimination.dir/main.c.o   -c /Users/guilhermemorais/CLionProjects/icc/replacement/main.c

CMakeFiles/gauss_elimination.dir/main.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/gauss_elimination.dir/main.c.i"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /Users/guilhermemorais/CLionProjects/icc/replacement/main.c > CMakeFiles/gauss_elimination.dir/main.c.i

CMakeFiles/gauss_elimination.dir/main.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/gauss_elimination.dir/main.c.s"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /Users/guilhermemorais/CLionProjects/icc/replacement/main.c -o CMakeFiles/gauss_elimination.dir/main.c.s

CMakeFiles/gauss_elimination.dir/libs/matrix.c.o: CMakeFiles/gauss_elimination.dir/flags.make
CMakeFiles/gauss_elimination.dir/libs/matrix.c.o: ../libs/matrix.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/guilhermemorais/CLionProjects/icc/replacement/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building C object CMakeFiles/gauss_elimination.dir/libs/matrix.c.o"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/gauss_elimination.dir/libs/matrix.c.o   -c /Users/guilhermemorais/CLionProjects/icc/replacement/libs/matrix.c

CMakeFiles/gauss_elimination.dir/libs/matrix.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/gauss_elimination.dir/libs/matrix.c.i"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /Users/guilhermemorais/CLionProjects/icc/replacement/libs/matrix.c > CMakeFiles/gauss_elimination.dir/libs/matrix.c.i

CMakeFiles/gauss_elimination.dir/libs/matrix.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/gauss_elimination.dir/libs/matrix.c.s"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /Users/guilhermemorais/CLionProjects/icc/replacement/libs/matrix.c -o CMakeFiles/gauss_elimination.dir/libs/matrix.c.s

# Object files for target gauss_elimination
gauss_elimination_OBJECTS = \
"CMakeFiles/gauss_elimination.dir/main.c.o" \
"CMakeFiles/gauss_elimination.dir/libs/matrix.c.o"

# External object files for target gauss_elimination
gauss_elimination_EXTERNAL_OBJECTS =

gauss_elimination: CMakeFiles/gauss_elimination.dir/main.c.o
gauss_elimination: CMakeFiles/gauss_elimination.dir/libs/matrix.c.o
gauss_elimination: CMakeFiles/gauss_elimination.dir/build.make
gauss_elimination: CMakeFiles/gauss_elimination.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/guilhermemorais/CLionProjects/icc/replacement/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking C executable gauss_elimination"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/gauss_elimination.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/gauss_elimination.dir/build: gauss_elimination

.PHONY : CMakeFiles/gauss_elimination.dir/build

CMakeFiles/gauss_elimination.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/gauss_elimination.dir/cmake_clean.cmake
.PHONY : CMakeFiles/gauss_elimination.dir/clean

CMakeFiles/gauss_elimination.dir/depend:
	cd /Users/guilhermemorais/CLionProjects/icc/replacement/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/guilhermemorais/CLionProjects/icc/replacement /Users/guilhermemorais/CLionProjects/icc/replacement /Users/guilhermemorais/CLionProjects/icc/replacement/cmake-build-debug /Users/guilhermemorais/CLionProjects/icc/replacement/cmake-build-debug /Users/guilhermemorais/CLionProjects/icc/replacement/cmake-build-debug/CMakeFiles/gauss_elimination.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/gauss_elimination.dir/depend

