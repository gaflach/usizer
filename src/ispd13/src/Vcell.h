/*
 *  Vcell.h
 *  sizer
 *
 *  Created by Tiago Reimann on 11/01/12.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef _VCELL_H_
#define _VCELL_H_
#include <vector>
#include <set>
#include <string>
#include "parser_helper.h"

#define CLOUD_TYPE unsigned int

using namespace std;

enum EdgeTypeEnum { RISE = 0, FALL = 1 };

class EdgeType {
	
	friend ostream &operator<<(ostream &out, const EdgeType &edgeType ) {
		if ( edgeType == FALL )
			out << "FALL";
		else if ( edgeType == RISE )
			out << "RISE";
		else
			out << "INVALID_EDGE";
		
		return out;
	} // end operator
	
private:
	int clsValue;
public:
	
	EdgeType() {clsValue = -1;}
	EdgeType(const EdgeTypeEnum edgeTypeEnum) {clsValue = edgeTypeEnum; }
	EdgeType(const int edgeType) {clsValue = edgeType;}
	
	operator int() { return clsValue; }
	operator int() const { return clsValue; }
	
	void reverse() { clsValue = clsValue == RISE ? FALL : RISE; }
	EdgeType getReversed() const { return clsValue == RISE ? FALL : RISE; }
};

class Vcell {
	
public:
	vector< Vcell* > previousCells,nextCells;
	string instName;
	string instType;
	
	// The cloud id this cell belongs to.
	//int cloud;
	//int changes;
	//int trials;
	//set<CLOUD_TYPE> clouds;
    
    // The logical depth indicates the maximum depth over all cell inputs. The
	// depth of an input is the number of steps required to reach such input
	// from a primary input or a sequential element (e.g flip-flop). All flip-
	// flops and input drivers have zero logical depth.
	int logicalDepth;
    
	// Similar to logical depth, however the depths are computed walking from
	// path tails to path heads. Flip-flops and input driver have a non-zero
	// reverse logical depth and the maximum reverse logical depth cell is
	// always a flip-flop or input driver.
	int reverseLogicalDepth;
	
	// The delay this cell is expected to have in order to match the circuit
	// timing requirements. This delay is computed so that each cell contributes
	// proportionately to the path delay. For instance, a cell at a branch-less
	// path with depth d should have an expected delay of T/d. However as most
	// paths have branch, the expected delay is adjusted so that at end of each
	// deeper path the sum of expected delays is equal to T.
	double expectedDelay;
	double expectedArrivalTime;
    
	// Total cell load capacitance (sinks, wires, primary outputs).
	double actualLoad;
	
	// Load capacitance caused by wires.
	double wireLoad;
    
	// Load capacitance caused by sink cells.
	double outputLoad;
	
	// Load capacitance associated to a primary output of the circuit. Zero for
	// cells not driving a primary output.
	double portLoad;
	
    //Branching effort of the cell ((Con + Coff)/Con)
    double branchingEffort;
    
	// Output slew based solely on the cell look-up table.
	double outputRiseSlew;
	double outputFallSlew;
	
	// Output delay based solely on the cell look-up table.
	union {
		struct {
			double outputRiseDelay;
			double outputFallDelay;
		};
		double outputDelay[2];
	};
	
	// Slacks calculated at this cell. (clock period minus output arrival time)
	double actualRiseSlack;
	double actualFallSlack;
	
	// The worst path slack passing through this cell. These slacks are computed
	// at the end of paths and propagated back to each cell on the path.
	//double worstRiseSlack;
	//double worstFallSlack;
	
	//slew target, delta slew target and maximum slew limit for the output pin p (p is induced by the attached sinks)
	double fallSlewTarget;
	double riseSlewTarget;
	double deltaFallSlewTarget;
    double deltaRiseSlewTarget;

	//slew estimative of the final slew in p' (pin of a input cell drive)
	double slewFallEstimative;	
    double slewRiseEstimative;	
	
	//Local criticality, lc(cell) is the difference between the slacks of the input pins, and the minimum of the slacks of its fanins
	double localFallCriticality;	
    double localRiseCriticality;
	
	//Global criticality of the output pin of the cell (slack (p) = required time of the pin (rat(p) - at (p)) arrival time of the pin (p) 
	double globalCriticalityOutputPinFall;
    double globalCriticalityOutputPinRise;
	
	double mustafaCriticality;
	
	// The worst path slack passing through this cell. These slacks are computed
	// at the end of paths and propagated back to each cell according to the path.
	union {
		struct {
			double worstRisePathSlack;
			double worstFallPathSlack;
		};
		double worstPathSlack[2];
	};
	
	// Auxiliary: ???
	double worstRiseDelayToEndpoint;
	double worstFallDelayToEndpoint;
	
	// This cell index in the Circuit::icells vector.
	int vectorIndex;
	
	// This cell index in the Circuit::depthSortedIndex.
	int depthIndex;
	
	// Timing net index.
	// [WARNING] Assuming cell has only one output.
	// [TODO] Support multiple outputs.
	int sinkNetIndex;
	
	//clock period
	//double clkPeriod;
	
	//next cells counter
	int nextNum;
	
	//previous cells counter
	//int prevNum;
    
    //LR lambdas
    //double lambdaDelay, lambdaMaxCap, lambdaMaxSlew;
	
	LibParserCellInfo * actualInstType;
	int actualInstTypeIndex, footprintIndex, pinOk; //pinOk => number of pins with inputSlew ready
	bool dontTouch, changed;
	
	// Auxiliary. Indicate how many times this cell has been visited in a walk.
	// For instance, when this count reaches the number of inputs of a cell,
	// we know that all paths to such cell have been covered.
	// [TODO] Put this flag outside the Vcell structure to make it thread safe.
	int visitcount;
	
	vector< pair<string, string> > pinNetPairs;
	vector< pair<double,double> > delays;					//first = rise, second = fall
	//vector< pair<double,double> > arrivalTimes;				//first = rise, second = fall
	vector< pair<double,double> > inputSlews;				//first = rise, second = fall
	
	string returnNetConnectedToPin( const string &pinName ) const;
	int returnPinIndex( const string &pinName ) const;
		
	int getPathDepth() const { return logicalDepth + reverseLogicalDepth; }
	
	double getLeakagePower() const { return actualInstType->leakagePower; }
	
	//Vcell* cpyCellInfo();
	void cellTiming();
	//void cellTiming2();
	//void newCellTiming(LibParserCellInfo * cellInst);
	//void setNewCellSize();
	bool needTimer;
    //bool skip;
    int treeIndex;
	
	Vcell () :
    treeIndex(-1),
    //lambdaDelay(0),
    //lambdaMaxCap(0),
    //lambdaMaxSlew(0),
	needTimer (false),
    pinOk (0),
    //prevNum (0),
    nextNum (0),
    //worstRiseSlack (0.0),
    //worstFallSlack (0.0),
    outputRiseSlew (0.0),
    outputFallSlew (0.0),
    portLoad (0.0),
    actualRiseSlack (0.0),
    actualFallSlack (0.0),
    actualLoad (0.0),
    actualInstTypeIndex (0),
    footprintIndex (0),
    wireLoad (0.0),
    outputLoad (0.0),
    dontTouch (false),
    changed (false),
    actualInstType (NULL),
    logicalDepth(-1),
    sinkNetIndex(-1),
    vectorIndex(-1),
    depthIndex(-1),
    expectedDelay(0.0),
    expectedArrivalTime(0.0),
	mustafaCriticality(0)
    //cloud(-1),
    //changes(0),
    //skip(false),
    //trials(0)
	{}
	
};

#endif // _VCELL_H_