rock_library(ls_pseudoinverse
    SOURCES Dummy.cpp
    HEADERS Dummy.hpp)

rock_executable(ls_pseudoinverse_bin Main.cpp
    DEPS ls_pseudoinverse)

