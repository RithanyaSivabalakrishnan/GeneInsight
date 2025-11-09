cmake_minimum_required(VERSION 3.5)
project(DNAQtTool)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTORCC ON)

find_package(Qt6 COMPONENTS Widgets REQUIRED)

qt_standard_project_setup()

qt_add_executable(DNAQtTool
    main.cpp
    mainwindow.cpp
    dna_analyzer.cpp
    mainwindow.h
    dna_analyzer.h
)

target_link_libraries(DNAQtTool PRIVATE Qt6::Widgets)
