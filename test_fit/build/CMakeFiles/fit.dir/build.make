# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.22

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/roushelfy/laser/test_fit

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/roushelfy/laser/test_fit/build

# Include any dependencies generated for this target.
include CMakeFiles/fit.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/fit.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/fit.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/fit.dir/flags.make

CMakeFiles/fit.dir/src/fit.cpp.o: CMakeFiles/fit.dir/flags.make
CMakeFiles/fit.dir/src/fit.cpp.o: ../src/fit.cpp
CMakeFiles/fit.dir/src/fit.cpp.o: CMakeFiles/fit.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/roushelfy/laser/test_fit/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/fit.dir/src/fit.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/fit.dir/src/fit.cpp.o -MF CMakeFiles/fit.dir/src/fit.cpp.o.d -o CMakeFiles/fit.dir/src/fit.cpp.o -c /home/roushelfy/laser/test_fit/src/fit.cpp

CMakeFiles/fit.dir/src/fit.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/fit.dir/src/fit.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/roushelfy/laser/test_fit/src/fit.cpp > CMakeFiles/fit.dir/src/fit.cpp.i

CMakeFiles/fit.dir/src/fit.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/fit.dir/src/fit.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/roushelfy/laser/test_fit/src/fit.cpp -o CMakeFiles/fit.dir/src/fit.cpp.s

# Object files for target fit
fit_OBJECTS = \
"CMakeFiles/fit.dir/src/fit.cpp.o"

# External object files for target fit
fit_EXTERNAL_OBJECTS =

fit: CMakeFiles/fit.dir/src/fit.cpp.o
fit: CMakeFiles/fit.dir/build.make
fit: CMakeFiles/fit.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/roushelfy/laser/test_fit/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable fit"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/fit.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/fit.dir/build: fit
.PHONY : CMakeFiles/fit.dir/build

CMakeFiles/fit.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/fit.dir/cmake_clean.cmake
.PHONY : CMakeFiles/fit.dir/clean

CMakeFiles/fit.dir/depend:
	cd /home/roushelfy/laser/test_fit/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/roushelfy/laser/test_fit /home/roushelfy/laser/test_fit /home/roushelfy/laser/test_fit/build /home/roushelfy/laser/test_fit/build /home/roushelfy/laser/test_fit/build/CMakeFiles/fit.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/fit.dir/depend

