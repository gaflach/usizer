/*
 *  Vcell.cpp
 *  sizer
 *
 *  Created by Tiago Reimann on 11/01/12.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 */

#include <cassert>
#include <float.h>
#include "Vcell.h"
#include "Circuit.h"

string Vcell::returnNetConnectedToPin( const string &pinName ) const {
	for ( int i = 0; i < pinNetPairs.size(); i++ )
		if ( pinNetPairs[i].first == pinName )
			return pinNetPairs[i].second;
#ifndef NDEBUG
	cerr << "[BUG] @ Vcell::returnNetConnectedToPin() - Pin name '" << pinName << "' not found.\n";
#endif	
	return "";
} // end method

// -----------------------------------------------------------------------------

int Vcell::returnPinIndex( const string &pinName ) const {
	for ( int i = 0; i < actualInstType->pins.size(); i++ )
		if ( actualInstType->pins[i].name == pinName )
			return i;
#ifndef NDEBUG
	cerr << "[BUG] @ Vcell::returnPinIndex() - Pin name '" << pinName << "' not found.\n";
#endif	
	return -1;
} // end method

// -----------------------------------------------------------------------------

void Vcell::cellTiming(){
	//calc timing for tmpCell
	//rise and fall separately
	//keeps worst delay and worst slew (according to ISPD slides)
	
	//char nada;		//for degug
	double x = 0.0, y = 0.0, outRiseDelay = 0.0, outFallDelay = 0.0, outRiseSlew = 0.0, outFallSlew = 0.0;//, zeroLoadRiseDelay = 0.0, zeroLoadFallDelay = 0.0;
	double weightX = 0.0, weightY = 0.0;
	double xLower, xUpper, yLower, yUpper; //x = load, y = slew
	int xLowerIndex, xUpperIndex, yLowerIndex, yUpperIndex, xLimit, yLimit;
	LibParserCellInfo * cellInst;
	
	/*double prevRiseSlew, prevFallSlew;
	
	prevRiseSlew = this->outputRiseSlew;
	prevFallSlew = this->outputFallSlew;*/
	
	this->outputRiseSlew = 0.0;
	this->outputFallSlew = 0.0;
	
	this->outputRiseDelay = 0.0;
	this->outputFallDelay = 0.0;
	
	cellInst = this->actualInstType;
	//if (tmpCell->instName == "inputDriver") cout << endl << endl << tmpCell->instName << "\t(" << cellInst->name << ")" << endl;
	
	for (int i = 0; i < cellInst->timingArcs.size(); ++i) {
		
		int j = 0;
		
		for (j = 0; j < this->pinNetPairs.size(); ++j)
		{
			
			if (cellInst->timingArcs[i].fromPin == this->pinNetPairs[j].first)
				break;
		}		
		assert(j != this->pinNetPairs.size());
		
		if (cellInst->timingArcs[i].toPin != "o") continue;
		/*
		if (!cellInst->isSequential)
		{
			assert(this->inputSlews[j].first > 0.0);
			assert(this->inputSlews[j].second > 0.0);
		}
		*/
		//this->delays[j].first = 0.0;
		//this->delays[j].second = 0.0;
		
		//cout << "wire name: " << tmpCell->pinNetPairs[j].second << endl;
		//cout << "out wire name: " << tmpCell->pinNetPairs[tmpCell->pinNetPairs.size()-1].second << endl;
		//cout << tmpCell->pinNetPairs.size() << " " << i << " from pin (lib): " << cellInst->timingArcs[i].fromPin << " from pin (instance): " << tmpCell->pinNetPairs[j].first << endl;
		//calc rise delay time
		
		// [CHECK] For a sequential cell, should be y = fallDelay? Or is it
		// riseDelay since ffs are rising edge?
		
		x = this->actualLoad;
		if (cellInst->isSequential == true) y = cellInst->timingArcs[i].riseDelay.transitionIndices[0];
		else y = this->inputSlews[j].second;
			
		//cout << " loads: " << tmpCell->wireLoad << " " << tmpCell->outputLoad << " " << tmpCell->portLoad << endl;
		//cout << " load (x): " << x << " transitions (y): " << tmpCell->inputSlews[j].first << " " << tmpCell->inputSlews[j].second << " " << j << endl;
		xLowerIndex = xUpperIndex = yLowerIndex = yUpperIndex = 0;
		xLimit = cellInst->timingArcs[i].riseDelay.loadIndices.size()-2;
		yLimit = cellInst->timingArcs[i].riseDelay.transitionIndices.size()-2;
		
		//no loads viol. are accepted -> not anymore
		while ( (xLowerIndex < xLimit) && ( cellInst->timingArcs[i].riseDelay.loadIndices[xLowerIndex+1] <= x) ) ++xLowerIndex;
		while ( (xUpperIndex <= xLimit) && ( cellInst->timingArcs[i].riseDelay.loadIndices[xUpperIndex] <= x) ) ++xUpperIndex;
		while ( (yLowerIndex < yLimit) && ( cellInst->timingArcs[i].riseDelay.transitionIndices[yLowerIndex+1] <= y) ) ++yLowerIndex;
		while ( (yUpperIndex <= yLimit) && ( cellInst->timingArcs[i].riseDelay.transitionIndices[yUpperIndex] <= y) ) ++yUpperIndex;
		
		xLower = cellInst->timingArcs[i].riseDelay.loadIndices[xLowerIndex];
		xUpper = cellInst->timingArcs[i].riseDelay.loadIndices[xUpperIndex];
		yLower = cellInst->timingArcs[i].riseDelay.transitionIndices[yLowerIndex];
		yUpper = cellInst->timingArcs[i].riseDelay.transitionIndices[yUpperIndex];
		
		//if (tmpCell->instName == "inputDriver") cout << " indices: " << xLower << " " << xUpper << " " << yLower << " " << yUpper << endl;
		
		weightX = (x-xLower)/(xUpper-xLower);
		weightY = (y-yLower)/(yUpper-yLower);
		
		//if (tmpCell->instName == "inputDriver") cout << " pesos: " << weightX << " " << weightY << endl;
		
		outRiseDelay = (1.0-weightX)*(1.0-weightY)*(cellInst->timingArcs[i].riseDelay.tableVals[xLowerIndex][yLowerIndex]);
		outRiseDelay += (weightX)*(1.0-weightY)*(cellInst->timingArcs[i].riseDelay.tableVals[xUpperIndex][yLowerIndex]);
		outRiseDelay += (1.0-weightX)*(weightY)*(cellInst->timingArcs[i].riseDelay.tableVals[xLowerIndex][yUpperIndex]);
		outRiseDelay += (weightX)*(weightY)*(cellInst->timingArcs[i].riseDelay.tableVals[xUpperIndex][yUpperIndex]);
		
		//if (tmpCell->instName == "inputDriver") cout << " outRiseDelay: " << outRiseDelay << " " << cellInst->timingArcs[i].riseDelay.tableVals[xLowerIndex][yLowerIndex] << endl;
		
		if (this->instName == "inputDriver") {
			x = 0;
			xLowerIndex = 0;
			xUpperIndex = 1;
			xLower = cellInst->timingArcs[i].riseDelay.loadIndices[xLowerIndex];
			xUpper = cellInst->timingArcs[i].riseDelay.loadIndices[xUpperIndex];
			
			//cout << " indices: " << xLower << " " << xUpper << " " << yLower << " " << yUpper << endl;
			
			weightX = (x-xLower)/(xUpper-xLower);
			weightY = (y-yLower)/(yUpper-yLower);
			
			//cout << " pesos: " << weightX << " " << weightY << endl;
			
			outRiseDelay -= (1.0-weightX)*(1.0-weightY)*(cellInst->timingArcs[i].riseDelay.tableVals[xLowerIndex][yLowerIndex]);
			outRiseDelay -= (weightX)*(1.0-weightY)*(cellInst->timingArcs[i].riseDelay.tableVals[xUpperIndex][yLowerIndex]);
			outRiseDelay -= (1.0-weightX)*(weightY)*(cellInst->timingArcs[i].riseDelay.tableVals[xLowerIndex][yUpperIndex]);
			outRiseDelay -= (weightX)*(weightY)*(cellInst->timingArcs[i].riseDelay.tableVals[xUpperIndex][yUpperIndex]);
			
			//cout << " outRiseDelay: " << outRiseDelay << " " << cellInst->timingArcs[i].riseDelay.tableVals[xLowerIndex][yLowerIndex] << endl;
				
		}
		//calc fall delay time
		
		x = this->actualLoad;
		if (cellInst->isSequential == true) y = cellInst->timingArcs[i].riseDelay.transitionIndices[0];
		else y = this->inputSlews[j].first;
		
		xLowerIndex = xUpperIndex = yLowerIndex = yUpperIndex = 0;
		xLimit = cellInst->timingArcs[i].fallDelay.loadIndices.size()-2;
		yLimit = cellInst->timingArcs[i].fallDelay.transitionIndices.size()-2;
		//no loads viol. are accepted
		
		while ( (xLowerIndex < xLimit) && (cellInst->timingArcs[i].fallDelay.loadIndices[xLowerIndex+1] <= x) ) ++xLowerIndex;
		while ( (xUpperIndex <= xLimit) && (cellInst->timingArcs[i].fallDelay.loadIndices[xUpperIndex] <= x) ) ++xUpperIndex;
		while ( (yLowerIndex < yLimit) && (cellInst->timingArcs[i].fallDelay.transitionIndices[yLowerIndex+1] <= y) ) ++yLowerIndex;
		while ( (yUpperIndex <= yLimit) && (cellInst->timingArcs[i].fallDelay.transitionIndices[yUpperIndex] <= y) ) ++yUpperIndex;
		
		xLower = cellInst->timingArcs[i].fallDelay.loadIndices[xLowerIndex];
		xUpper = cellInst->timingArcs[i].fallDelay.loadIndices[xUpperIndex];
		yLower = cellInst->timingArcs[i].fallDelay.transitionIndices[yLowerIndex];
		yUpper = cellInst->timingArcs[i].fallDelay.transitionIndices[yUpperIndex];
		
		weightX = (x-xLower)/(xUpper-xLower);
		weightY = (y-yLower)/(yUpper-yLower);
		
		outFallDelay = (1.0-weightX)*(1.0-weightY)*(cellInst->timingArcs[i].fallDelay.tableVals[xLowerIndex][yLowerIndex]);
		outFallDelay += (weightX)*(1.0-weightY)*(cellInst->timingArcs[i].fallDelay.tableVals[xUpperIndex][yLowerIndex]);
		outFallDelay += (1.0-weightX)*(weightY)*(cellInst->timingArcs[i].fallDelay.tableVals[xLowerIndex][yUpperIndex]);
		outFallDelay += (weightX)*(weightY)*(cellInst->timingArcs[i].fallDelay.tableVals[xUpperIndex][yUpperIndex]);
		
		if (this->instName == "inputDriver") {
			// [NOTE] Tiago, eu movi novamente pra ca porque parece que
			// potencialmente as variaveis começando com y podem ter valores
			// errados.
			
			x = 0;
			xLowerIndex = 0;
			xUpperIndex = 1;
			xLower = cellInst->timingArcs[i].fallDelay.loadIndices[xLowerIndex];
			xUpper = cellInst->timingArcs[i].fallDelay.loadIndices[xUpperIndex];
			
			//out << " indices: " << xLower << " " << xUpper << " " << yLower << " " << yUpper << endl;
			
			weightX = (x-xLower)/(xUpper-xLower);
			weightY = (y-yLower)/(yUpper-yLower);
			
			outFallDelay -= (1.0-weightX)*(1.0-weightY)*(cellInst->timingArcs[i].fallDelay.tableVals[xLowerIndex][yLowerIndex]);
			outFallDelay -= (weightX)*(1.0-weightY)*(cellInst->timingArcs[i].fallDelay.tableVals[xUpperIndex][yLowerIndex]);
			outFallDelay -= (1.0-weightX)*(weightY)*(cellInst->timingArcs[i].fallDelay.tableVals[xLowerIndex][yUpperIndex]);
			outFallDelay -= (weightX)*(weightY)*(cellInst->timingArcs[i].fallDelay.tableVals[xUpperIndex][yUpperIndex]);
		}
		
		/*
		if ( (this->delays[j].first != outRiseDelay) || (this->delays[j].second != outFallDelay) )
		{
			cout << "delay changed: " << this->instName << endl;
			cout << " rise: " << this->delays[j].first  << " -> " << outRiseDelay << endl;
			cout << " fall: " << this->delays[j].second  << " -> " << outFallDelay << endl;
			cout << " x (Load): " << x << " y (Slew): " << y << endl;
			cout << " indices: " << xLower << " " << xUpper << " " << yLower << " " << yUpper << endl;
			
			cout << " pesos: " << weightX << " " << weightY << endl;
		}
		 */
		
		this->delays[j].first = outRiseDelay;
		this->delays[j].second = outFallDelay;
		
		//assert(this->delays[j].first > 0.0);
		//assert(this->delays[j].second > 0.0);
		
		//cout << "delays: " << tmpCell->delays[j].first << "\t" << tmpCell->delays[j].second << endl;
		
		//calc rise slew time
		
		x = this->actualLoad;
		if (cellInst->isSequential == true) y = cellInst->timingArcs[i].riseDelay.transitionIndices[0];
		else y = this->inputSlews[j].second;
		
		xLowerIndex = xUpperIndex = yLowerIndex = yUpperIndex = 0;
	    xLimit = cellInst->timingArcs[i].riseTransition.loadIndices.size()-2;
	    yLimit = cellInst->timingArcs[i].riseTransition.transitionIndices.size()-2;
		
		//no loads viol. are accepted
		while ( (xLowerIndex < xLimit) && ( cellInst->timingArcs[i].riseTransition.loadIndices[xLowerIndex+1] <= x) ) ++xLowerIndex;
		while ( (xUpperIndex <= xLimit) && ( cellInst->timingArcs[i].riseTransition.loadIndices[xUpperIndex] <= x) ) ++xUpperIndex;
		while ( (yLowerIndex < yLimit) && ( cellInst->timingArcs[i].riseTransition.transitionIndices[yLowerIndex+1] <= y) ) ++yLowerIndex;
		while ( (yUpperIndex <= yLimit) && ( cellInst->timingArcs[i].riseTransition.transitionIndices[yUpperIndex] <= y) ) ++yUpperIndex;
		
		xLower = cellInst->timingArcs[i].riseTransition.loadIndices[xLowerIndex];
		xUpper = cellInst->timingArcs[i].riseTransition.loadIndices[xUpperIndex];
		yLower = cellInst->timingArcs[i].riseTransition.transitionIndices[yLowerIndex];
		yUpper = cellInst->timingArcs[i].riseTransition.transitionIndices[yUpperIndex];
		
		weightX = (x-xLower)/(xUpper-xLower);
		weightY = (y-yLower)/(yUpper-yLower);
		
		//cout << " indices: " << xLower << " " << xUpper << " " << yLower << " " << yUpper << endl;
		
		//cout << " pesos: " << weightX << " " << weightY << endl;
		
		outRiseSlew = (1.0-weightX)*(1.0-weightY)*(cellInst->timingArcs[i].riseTransition.tableVals[xLowerIndex][yLowerIndex]);
		outRiseSlew += (weightX)*(1.0-weightY)*(cellInst->timingArcs[i].riseTransition.tableVals[xUpperIndex][yLowerIndex]);
		outRiseSlew += (1.0-weightX)*(weightY)*(cellInst->timingArcs[i].riseTransition.tableVals[xLowerIndex][yUpperIndex]);
		outRiseSlew += (weightX)*(weightY)*(cellInst->timingArcs[i].riseTransition.tableVals[xUpperIndex][yUpperIndex]);
		
		//calc fall slew time
		
		x = this->actualLoad;
		if (cellInst->isSequential == true) y = cellInst->timingArcs[i].riseDelay.transitionIndices[0];
		else y = this->inputSlews[j].first;
		
		xLowerIndex = xUpperIndex = yLowerIndex = yUpperIndex = 0;
		xLimit = cellInst->timingArcs[i].fallTransition.loadIndices.size()-2;
		yLimit = cellInst->timingArcs[i].fallTransition.transitionIndices.size()-2;
		
		//no loads viol. are accepted
		while ( (xLowerIndex < xLimit) && ( cellInst->timingArcs[i].fallTransition.loadIndices[xLowerIndex+1] <= x) ) ++xLowerIndex;
		while ( (xUpperIndex <= xLimit) && ( cellInst->timingArcs[i].fallTransition.loadIndices[xUpperIndex] <= x) ) ++xUpperIndex;
		while ( (yLowerIndex < yLimit) && ( cellInst->timingArcs[i].fallTransition.transitionIndices[yLowerIndex+1] <= y) ) ++yLowerIndex;
		while ( (yUpperIndex <= yLimit) && ( cellInst->timingArcs[i].fallTransition.transitionIndices[yUpperIndex] <= y) ) ++yUpperIndex;
		
		xLower = cellInst->timingArcs[i].fallTransition.loadIndices[xLowerIndex];
		xUpper = cellInst->timingArcs[i].fallTransition.loadIndices[xUpperIndex];
		yLower = cellInst->timingArcs[i].fallTransition.transitionIndices[yLowerIndex];
		yUpper = cellInst->timingArcs[i].fallTransition.transitionIndices[yUpperIndex];
		
		weightX = (x-xLower)/(xUpper-xLower);
		weightY = (y-yLower)/(yUpper-yLower);
		
		outFallSlew = (1.0-weightX)*(1.0-weightY)*(cellInst->timingArcs[i].fallTransition.tableVals[xLowerIndex][yLowerIndex]);
		outFallSlew += (weightX)*(1.0-weightY)*(cellInst->timingArcs[i].fallTransition.tableVals[xUpperIndex][yLowerIndex]);
		outFallSlew += (1.0-weightX)*(weightY)*(cellInst->timingArcs[i].fallTransition.tableVals[xLowerIndex][yUpperIndex]);
		outFallSlew += (weightX)*(weightY)*(cellInst->timingArcs[i].fallTransition.tableVals[xUpperIndex][yUpperIndex]);
		
		this->outputRiseSlew = max(this->outputRiseSlew,outRiseSlew);
		this->outputFallSlew = max(this->outputFallSlew,outFallSlew);
		
		this->outputRiseDelay = max(this->outputRiseDelay,outRiseDelay);
		this->outputFallDelay = max(this->outputFallDelay,outFallDelay);
		
		//assert(this->outputRiseSlew > 0.0);
		//assert(this->outputFallSlew > 0.0);
		//if (tmpCell->instName == "inputDriver") cout << " slews: " << outFallSlew << " " << outRiseSlew << endl;
		//cout << " slews: " << outRiseSlew << " " << outFallSlew << endl;
		//cout << " done" <<endl;
		//if (tmpCell->instName != "inputDriver") cin >> nada;
		//cout << " delays: " << tmpCell->delays[j].first << " " << tmpCell->delays[j].second << "\t" << j << endl;
		
	}	
	/*
	if ( (this->outputRiseSlew != prevRiseSlew) || (this->outputFallSlew != prevFallSlew) )
	{
		cout << "slew changed: " << this->instName << endl;
		cout << " rise: " << prevRiseSlew  << " -> " << this->outputRiseSlew << endl;
		cout << " fall: " << prevFallSlew  << " -> " << this->outputFallSlew << endl;
		cout << " x (Load): " << x << " y (Slew): " << y << endl;
		cout << " indices: " << xLower << " " << xUpper << " " << yLower << " " << yUpper << endl;
		
		cout << " pesos: " << weightX << " " << weightY << endl;
	}
	*/
}

// -----------------------------------------------------------------------------


