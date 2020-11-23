#
# Generated Makefile - do not edit!
#
# Edit the Makefile in the project folder instead (../Makefile). Each target
# has a -pre and a -post target defined where you can add customized code.
#
# This makefile implements configuration specific macros and targets.


# Environment
MKDIR=mkdir
CP=cp
GREP=grep
NM=nm
CCADMIN=CCadmin
RANLIB=ranlib
CC=gcc
CCC=g++
CXX=g++
FC=gfortran
AS=as

# Macros
CND_PLATFORM=GNU-Linux
CND_DLIB_EXT=so
CND_CONF=Debug
CND_DISTDIR=dist
CND_BUILDDIR=build

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/_ext/6cfc6540/Circuit.o \
	${OBJECTDIR}/_ext/6cfc6540/NewtonRaphson.o \
	${OBJECTDIR}/_ext/6cfc6540/PracticalSocket.o \
	${OBJECTDIR}/_ext/6cfc6540/Vcell.o \
	${OBJECTDIR}/_ext/6cfc6540/global.o \
	${OBJECTDIR}/_ext/6cfc6540/main.o \
	${OBJECTDIR}/_ext/6cfc6540/parser_helper.o \
	${OBJECTDIR}/_ext/6cfc6540/timer_interface.o


# C Compiler Flags
CFLAGS=

# CC Compiler Flags
CCFLAGS=-Wreturn-type -Wno-unused-result
CXXFLAGS=-Wreturn-type -Wno-unused-result

# Fortran Compiler Flags
FFLAGS=

# Assembler Flags
ASFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/ispd13

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/ispd13: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${LINK.cc} -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/ispd13 ${OBJECTFILES} ${LDLIBSOPTIONS}

${OBJECTDIR}/_ext/6cfc6540/Circuit.o: ../../src/ispd13/src/Circuit.cpp
	${MKDIR} -p ${OBJECTDIR}/_ext/6cfc6540
	${RM} "$@.d"
	$(COMPILE.cc) -g -D_GLIBCXX_DEBUG -I../../src -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/6cfc6540/Circuit.o ../../src/ispd13/src/Circuit.cpp

${OBJECTDIR}/_ext/6cfc6540/NewtonRaphson.o: ../../src/ispd13/src/NewtonRaphson.cpp
	${MKDIR} -p ${OBJECTDIR}/_ext/6cfc6540
	${RM} "$@.d"
	$(COMPILE.cc) -g -D_GLIBCXX_DEBUG -I../../src -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/6cfc6540/NewtonRaphson.o ../../src/ispd13/src/NewtonRaphson.cpp

${OBJECTDIR}/_ext/6cfc6540/PracticalSocket.o: ../../src/ispd13/src/PracticalSocket.cpp
	${MKDIR} -p ${OBJECTDIR}/_ext/6cfc6540
	${RM} "$@.d"
	$(COMPILE.cc) -g -D_GLIBCXX_DEBUG -I../../src -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/6cfc6540/PracticalSocket.o ../../src/ispd13/src/PracticalSocket.cpp

${OBJECTDIR}/_ext/6cfc6540/Vcell.o: ../../src/ispd13/src/Vcell.cpp
	${MKDIR} -p ${OBJECTDIR}/_ext/6cfc6540
	${RM} "$@.d"
	$(COMPILE.cc) -g -D_GLIBCXX_DEBUG -I../../src -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/6cfc6540/Vcell.o ../../src/ispd13/src/Vcell.cpp

${OBJECTDIR}/_ext/6cfc6540/global.o: ../../src/ispd13/src/global.cpp
	${MKDIR} -p ${OBJECTDIR}/_ext/6cfc6540
	${RM} "$@.d"
	$(COMPILE.cc) -g -D_GLIBCXX_DEBUG -I../../src -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/6cfc6540/global.o ../../src/ispd13/src/global.cpp

${OBJECTDIR}/_ext/6cfc6540/main.o: ../../src/ispd13/src/main.cpp
	${MKDIR} -p ${OBJECTDIR}/_ext/6cfc6540
	${RM} "$@.d"
	$(COMPILE.cc) -g -D_GLIBCXX_DEBUG -I../../src -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/6cfc6540/main.o ../../src/ispd13/src/main.cpp

${OBJECTDIR}/_ext/6cfc6540/parser_helper.o: ../../src/ispd13/src/parser_helper.cpp
	${MKDIR} -p ${OBJECTDIR}/_ext/6cfc6540
	${RM} "$@.d"
	$(COMPILE.cc) -g -D_GLIBCXX_DEBUG -I../../src -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/6cfc6540/parser_helper.o ../../src/ispd13/src/parser_helper.cpp

${OBJECTDIR}/_ext/6cfc6540/timer_interface.o: ../../src/ispd13/src/timer_interface.cpp
	${MKDIR} -p ${OBJECTDIR}/_ext/6cfc6540
	${RM} "$@.d"
	$(COMPILE.cc) -g -D_GLIBCXX_DEBUG -I../../src -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/6cfc6540/timer_interface.o ../../src/ispd13/src/timer_interface.cpp

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
