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
	${OBJECTDIR}/_ext/1854f76a/Circuit.o \
	${OBJECTDIR}/_ext/1854f76a/Vcell.o \
	${OBJECTDIR}/_ext/1854f76a/global.o \
	${OBJECTDIR}/_ext/1854f76a/main.o \
	${OBJECTDIR}/_ext/1854f76a/parser_helper.o \
	${OBJECTDIR}/_ext/1854f76a/timer_interface.o


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
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/ispd12

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/ispd12: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${LINK.cc} -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/ispd12 ${OBJECTFILES} ${LDLIBSOPTIONS}

${OBJECTDIR}/_ext/1854f76a/Circuit.o: ../../src/ispd12/Circuit.cpp
	${MKDIR} -p ${OBJECTDIR}/_ext/1854f76a
	${RM} "$@.d"
	$(COMPILE.cc) -g -D_GLIBCXX_DEBUG -I../../src -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1854f76a/Circuit.o ../../src/ispd12/Circuit.cpp

${OBJECTDIR}/_ext/1854f76a/Vcell.o: ../../src/ispd12/Vcell.cpp
	${MKDIR} -p ${OBJECTDIR}/_ext/1854f76a
	${RM} "$@.d"
	$(COMPILE.cc) -g -D_GLIBCXX_DEBUG -I../../src -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1854f76a/Vcell.o ../../src/ispd12/Vcell.cpp

${OBJECTDIR}/_ext/1854f76a/global.o: ../../src/ispd12/global.cpp
	${MKDIR} -p ${OBJECTDIR}/_ext/1854f76a
	${RM} "$@.d"
	$(COMPILE.cc) -g -D_GLIBCXX_DEBUG -I../../src -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1854f76a/global.o ../../src/ispd12/global.cpp

${OBJECTDIR}/_ext/1854f76a/main.o: ../../src/ispd12/main.cpp
	${MKDIR} -p ${OBJECTDIR}/_ext/1854f76a
	${RM} "$@.d"
	$(COMPILE.cc) -g -D_GLIBCXX_DEBUG -I../../src -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1854f76a/main.o ../../src/ispd12/main.cpp

${OBJECTDIR}/_ext/1854f76a/parser_helper.o: ../../src/ispd12/parser_helper.cpp
	${MKDIR} -p ${OBJECTDIR}/_ext/1854f76a
	${RM} "$@.d"
	$(COMPILE.cc) -g -D_GLIBCXX_DEBUG -I../../src -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1854f76a/parser_helper.o ../../src/ispd12/parser_helper.cpp

${OBJECTDIR}/_ext/1854f76a/timer_interface.o: ../../src/ispd12/timer_interface.cpp
	${MKDIR} -p ${OBJECTDIR}/_ext/1854f76a
	${RM} "$@.d"
	$(COMPILE.cc) -g -D_GLIBCXX_DEBUG -I../../src -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1854f76a/timer_interface.o ../../src/ispd12/timer_interface.cpp

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
