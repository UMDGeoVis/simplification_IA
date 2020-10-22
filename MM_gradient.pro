TEMPLATE = app
CONFIG += console
CONFIG += c++11

# Directories
DESTDIR = dist/
OBJECTS_DIR = build/

#LIBS += -fopenmp


QMAKE_CXXFLAGS_RELEASE += -fpermissive
QMAKE_CXXFLAGS_DEBUG += -fpermissive
QMAKE_CXXFLAGS_RELEASE -= -O2
QMAKE_CXXFLAGS_RELEASE += -O3
#QMAKE_CXXFLAGS_RELEASE += -fopenmp

INCLUDEPATH += -L/usr/lib/


SOURCES += source/main.cpp \
    source/Lib_Tri/Vertex3D.cpp \
    source/Lib_Tri/Vertex2D.cpp \
    source/Lib_Tri/Triangle.cpp \
    source/Lib_Tri/Timer.cpp \
    source/Lib_Tri/Reader.cpp \
    source/Lib_Tri/formangradientvector.cpp \
    source/Lib_Tri/Edge.cpp \
    source/Lib_Tri/forman_feature.cpp \
    source/Lib_Tri/edge_collapse.cpp \
    source/Lib_Tri/IO.cpp \
    source/Lib_Tri/forman_refinement.cpp \
    source/Lib_Tri/forman_extraction_query.cpp \
    source/Lib_Tri/forman_simplify_mr.cpp \
    source/Lib_Tri/forman_simplify.cpp

HEADERS += \
    source/Lib_Tri/Vertex3D.h \
    source/Lib_Tri/Vertex2D.h \
    source/Lib_Tri/Triangle.h \
    source/Lib_Tri/Timer.h \
    source/Lib_Tri/Sorting.h \
    source/Lib_Tri/Reader.h \
    source/Lib_Tri/Mesh.h \
    source/Lib_Tri/formangradientvector.h \
    source/Lib_Tri/Edge.h \
    source/Lib_Tri/forman_arrow.h \
    source/Ig/ig.h \
    source/Dag/Dag.h \
    source/Lib_Tri/Matrix.h


