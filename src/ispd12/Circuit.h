/*
 *  Circuit.h
 *  sizer
 *
 *  Created by Tiago Reimann on 11/01/12.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 */

#include <vector>
#include <set>
#include <queue>
#include <stack>
#include <deque>
#include <bitset>
#include <string>
#include <float.h>
#include <climits>
#include <cstdlib>
#include "Vcell.h"
#include "global.h"
#include <cmath>
#include "Stopwatch.h"

#ifdef PARALLEL
#include <Poco/ThreadPool.h>
#include <Poco/Runnable.h>
#include <Poco/Environment.h>
#endif

using namespace std;

#define MAKE_SELF_OPERATOR( OP ) \
friend void operator OP ( EdgeArray<T> &v0, const EdgeArray<T> v1 ) { v0[RISE] OP v1[RISE], v0[FALL] OP v1[FALL]; } \
friend void operator OP ( EdgeArray<T> &v0, const T            v1 ) { v0[RISE] OP v1; v0[FALL] OP v1; }

#define MAKE_OPERATOR( OP ) \
friend EdgeArray<T> operator OP ( const EdgeArray<T> v0, const EdgeArray<T> v1 ) { return EdgeArray<T>(v0[RISE] OP v1[RISE], v0[FALL] OP v1[FALL]); } \
friend EdgeArray<T> operator OP ( const T            v0, const EdgeArray<T> v1 ) { return EdgeArray<T>(v0       OP v1[RISE], v0       OP v1[FALL]); } \
friend EdgeArray<T> operator OP ( const EdgeArray<T> v0, const T            v1 ) { return EdgeArray<T>(v0[RISE] OP v1      , v0[FALL] OP v1      ); }

#ifndef NDEBUG
template<typename T>
class array : public vector<T> {
	public:
    T&  operator[](int __n) { return vector<T>::at(__n); }
    const T& operator[](int __n) const { return vector<T>::at(__n); }		
};
#else
template<typename T>
class array : public vector<T> {};
#endif
	
typedef double Number;	

class Circuit {
		
public:		

	//
	// Baseline flow (ISVLSI 2013) and configuration.
	//
	
	bool optIgnoreLeakagePower;
	bool optEnableLambdaDelaySensitivities;
	bool optEnableSlackFiltering;
	bool optEnableGamma;
	bool optEnableLoadViolationFiltering;
	
	bool optUseTennakoon;
	
	bool optRunInitialSizing;
	bool optRunInitialSizingForLoadAndSlewViolationRemoval;
	bool optRunPowerRecovery;
	bool optRunTimingRecovery;

	void runBaselineFlow();
	
	Stopwatch stopwatchUpdateTiming;
	
private:	
	
	template<typename T>
	class EdgeArray {
		
		friend ostream &operator<<(ostream &out, const EdgeArray<T> array ) {
			return out << "(" << array[RISE] << ", " << array[FALL] << ")";
		} // end operator
		
		MAKE_OPERATOR(+);
		MAKE_OPERATOR(-);
		MAKE_OPERATOR(*);
		MAKE_OPERATOR(/);
		
		MAKE_SELF_OPERATOR(+=);
		MAKE_SELF_OPERATOR(-=);
		MAKE_SELF_OPERATOR(*=);
		MAKE_SELF_OPERATOR(/=);
		
		friend EdgeArray<T> operator-( const EdgeArray<T> &v0 ) { return EdgeArray<T>(-v0[RISE], -v0[FALL]); }
		
		friend EdgeArray<T> max( const EdgeArray<T> v0, const EdgeArray<T> v1 ) { return EdgeArray<T>(max(v0[RISE],v1[RISE]),max(v0[FALL],v1[FALL])); }
		friend EdgeArray<T> min( const EdgeArray<T> v0, const EdgeArray<T> v1 ) { return EdgeArray<T>(min(v0[RISE],v1[RISE]),min(v0[FALL],v1[FALL])); }
		
		friend EdgeArray<T> abs( const EdgeArray<T> v ) { return EdgeArray<T>(fabs(v[RISE]),fabs(v[FALL])); }
		friend EdgeArray<T> pow( const EdgeArray<T> v, const double exp ) { return EdgeArray<T>(pow(v[RISE], exp),pow(v[FALL], exp)); }
		friend EdgeArray<T> sqrt( const EdgeArray<T> v) { return EdgeArray<T>(sqrt(v[RISE]),sqrt(v[FALL])); }
        
	private:
		T clsValue[2];
	public:
		T &operator[]( const EdgeType edgeType ) {return clsValue[edgeType];}
		T  operator[]( const EdgeType edgeType ) const {return clsValue[edgeType];}
		
		EdgeArray &operator=(const EdgeArray &array) {
			clsValue[RISE] = array[RISE];
			clsValue[FALL] = array[FALL];
			return *this;
		} // end operator
		
		EdgeArray(const T rise, const T fall ) {
			clsValue[RISE] = rise;
			clsValue[FALL] = fall;
		} // end constructor
		
		EdgeArray(){};
		
		void set( const T rise, const T fall ) {
			clsValue[RISE] = rise;
			clsValue[FALL] = fall;
		} // end method
		
		T getMax() const { return max(clsValue[RISE], clsValue[FALL]); }
		T getMin() const { return min(clsValue[RISE], clsValue[FALL]); }
        T getRise() const { return clsValue[RISE]; }
        T getFall() const { return clsValue[FALL]; }
		
		EdgeArray<T> getReversed() const { return EdgeArray(getFall(), getRise()); }
		
		void reverse() { return swap(clsValue[RISE], clsValue[FALL]); }
		
		double aggregate() const { return clsValue[RISE] + clsValue[FALL]; }
	}; // end class
    
    void setDeltaLoadSlew(const int netIndex);

    //--------------------------------------------------------------------------
	// Reporting
	//--------------------------------------------------------------------------	
	
	ofstream reportFile_LogicalEffortDiscretization;
	ofstream reportFile_LogicalEffortPropagation;
	ofstream reportFile_LogicalEffortCharacterization;
	
	static const bool enableReport_LogicalEffortDiscretization;
	static const bool enableReport_LogicalEffortPropagation;
	static const bool enableReport_LogicalEffortCharacterization;

    //--------------------------------------------------------------------------
	// Analysis
	//--------------------------------------------------------------------------	
public:		

	void analysis_FanoutDistribution();
	void analysis_Cells();

private:	
    //--------------------------------------------------------------------------
	// Library Analysis and Characterization
	//--------------------------------------------------------------------------	
	
	// [NOTE] HARD CODED: This library characterization is rightly tied to the
	// library from the ISPD contests. Note however, that it can be ported to
	// handle other libraries in a generic way.
	//
	// Inside the library (orgCells.oCells), cells are organized by footprint
	// (function). Every footprint may have several cell types or versions which
	// are not organized in any specific ordering.
	//
	// The library characterization (timingLibraryCharacterization) allows to
	// get the vth index of a cell as well as its size index i.e. relative size
	// w.r.t. the other cells with the same footprint.
	//
	// Similarly one can get the cell type index i.e. the cell index inside its
	// footprint group.
	
	struct LogicalEffort {
		EdgeArray<double> cin; // input capacitance
		EdgeArray<double> g; // logical effort
		EdgeArray<double> p; // parasitic delay
		
		// For arcs, this flag is set to one if the arc is valid and to zero
		// otherwise. For arc groups, this flag counts the number of valid arcs.
		// Note that if the number of valid arcs is zero, this flag is also zero
		// so that is safe to ask if(valid).
		// An invalid arc may happen for some non-timing arcs (e.g. CK->D).
		EdgeArray<int> valid;
		
		// Residuum from least squares.
		EdgeArray<double> residuum;
		
		LogicalEffort() : cin(0, 0), g(0, 0), p(0, 0), valid(0, 0), residuum(0, 0) {
		} // end constructor
	}; // end struct
	
	struct LibArcCharacterization : LogicalEffort{
	}; // end struct
	
	struct LibCellCharacterization {
		// Vth = 0 --> higher Vth, slowest
		// Vth = n --> lower Vth, fastest
		int vth;
		int size;
		
		vector<LibArcCharacterization> arcs;

		LibCellCharacterization() {
			vth = -1;
			size = -1;
		} // end constructor
	}; // end struct
	
	struct LibCellGroupCharacterization {
		// Get the cell type index inside its footprint group by its vth and
		// size.
		// Usage: [vth][size] -> library cell type index
		vector< vector<int> > mappingVthSizeToCellTypeIndex;
		
		// Get the cell characterization by its index inside its footprint
		// group.
		// Usage: [library cell type index] -> characterization
		vector< LibCellCharacterization > cellCharacterization;
		
		// Get logical effort by vth and arc.
		// Usage: [vth][arc] -> logical effort
		vector< vector<LogicalEffort> > arcs;
	}; // end struct
	
	vector<LibCellGroupCharacterization> timingLibraryCharacterization;	
	
	LibParserCellInfo &getLibCell(const int footprint, const int vth, const int size) {
		return orgCells.oCells[footprint].cells[
				timingLibraryCharacterization[footprint].mappingVthSizeToCellTypeIndex[vth][size]];
	} // end method

	LibCellCharacterization &getLibCellCharacterization(const int footprint, const int type) {
		return timingLibraryCharacterization[footprint].cellCharacterization[type];
	} // end method

	LibCellCharacterization &getLibCellCharacterization(Vcell * cell) {
		return getLibCellCharacterization(cell->footprintIndex, cell->actualInstTypeIndex);
	} // end method

	LibCellGroupCharacterization &getLibCellGroupCharacterization(Vcell * cell) {
		return timingLibraryCharacterization[cell->footprintIndex];
	} // end method

	class TimingArc;
	LogicalEffort &getLibArcLogicalEffortForCellGroup(const TimingArc &arc) {
		LibCellGroupCharacterization &groupChar = getLibCellGroupCharacterization(arc.cell);
		LibCellCharacterization &cellChar = getLibCellCharacterization(arc.cell);
		return groupChar.arcs[cellChar.vth][arc.lut];
	} // end method
	
	LogicalEffort &getLibArcLogicalEffortForCellSpecific(const TimingArc &arc) {
		LibCellCharacterization &cellChar = getLibCellCharacterization(arc.cell);
		return cellChar.arcs[arc.lut];
	} // end method	
	
    //--------------------------------------------------------------------------
	// Logical Effort
	//--------------------------------------------------------------------------
	
	struct LogicalEffortPath {
		double pathLogicalEffort;  // G
		double pathBranching;      // B
		double pathParasiticDelay; // P
		double pathCin;            // input capacitance at the path driver
		int pathDepth;		
		
		LogicalEffortPath() {
			pathLogicalEffort = numeric_limits<double>::quiet_NaN();
			pathBranching = numeric_limits<double>::quiet_NaN();
			pathParasiticDelay = numeric_limits<double>::quiet_NaN();
			pathCin = numeric_limits<double>::quiet_NaN();
			pathDepth = -1;
		}
	};
	
    //--------------------------------------------------------------------------
	// Timing
	//--------------------------------------------------------------------------
	
	// In order to understand better this structure, imagine the circuit netlist
	// represented as a directed graph where nets are the graph nodes (that's
	// right!) and timing arcs are the directed edges. Note that cells are not
	// directly represented in this graph.
	
	// [NOTE] Do not include any unnecessary property for timing calculation
	// purposes in TimingArc and TimingNet structures. They must be kept as
	// small as possible to ensure the update timing methods will run fast.
	// Those methods were in some extent designed to be cache-aware so that
	// smaller data structures means a great chance that the next required data
	// will be already at cache. Also less memory access is required too.
	
	struct TimingArc {
		Vcell * cell; // pointer to the node to which this arc belongs to
		int lut;      // index of the timing arc inside the LibParserCellInfo
		int pin;      // index of the pin inside the LibParserCellInfo
		int driver;   // index of the net which drives this arc
		int sink;     // index of the net which is driven by this arc
		
		TimingArc() {
			cell = NULL;
			lut = -1;
			pin = -1;
			driver = -1;
			sink = -1;
		} // end constructor
	}; // end struct
	
	struct TimingNet {
		Vcell * driver;
		int fanout;      // number of timing arcs driven by this net (used in the slew violation calculation)
		int depth;       // logical depth of this net
		
		int ipaths;      // number of paths from inputs to this net
		int opaths;      // number of paths from outputs to this net

		float wipaths;   // similar to ipaths, but endpoints may get a weight != 1
		float wopaths;   // similar to opaths, but endpoints may get a weight != 1

		float centrality;
		
		TimingNet() {
			driver = NULL;
			depth = -1;
			fanout = 0;
			ipaths = 0;
			opaths = 0;
			wipaths = 0;
			wopaths = 0;
			centrality = 0;
		} // end constructor
	}; // end struct
	
	struct TimingArcState {
		EdgeArray<double> islew; //  input slew at this arc
		EdgeArray<double> oslew; // output slew at this arc
		
		EdgeArray<double> delay; // delay of this arc
        EdgeArray<double> lambda; // Lagrange Relaxation
		
		EdgeArray<double> arrivalTime; // arrival time at the input of this arc
		
		//Reimann
		EdgeArray<double> slack; // slack (m_u->v) (worst arrival time + arc delay + worst delay to endpoint) [TODO] consider moving to outside
		
		// Logical Effort
		EdgeArray<LogicalEffortPath> logicalEffortPath;
		double logicalEffortCin;
		EdgeArray<double> logicalEffortArrivalTime; // arrival time at the input of this arc
		EdgeArray<double> logicalEffortDelay; // current delay of this arc
		
		TimingArcState() {
            lambda.set(-1.0,-1.0);
			logicalEffortCin = numeric_limits<double>::quiet_NaN();
			logicalEffortArrivalTime.set(numeric_limits<double>::quiet_NaN(),
					numeric_limits<double>::quiet_NaN());
			logicalEffortDelay.set(numeric_limits<double>::quiet_NaN(),
					numeric_limits<double>::quiet_NaN());			
		} // end constructor		
	};	
	
	struct TimingNetState {
		double load;                    // total net load (due to wires, input pins and primary outputs)
		EdgeArray<double> slew;         // slew in this net
		EdgeArray<double> arrivalTime;  // arrival time at this net
		EdgeArray<int> backtrack;       // index of net driving the driver input with greatest arrival time
		EdgeArray<int> backtrackSlew;   // index of net driving the driver input with greatest slew
		
		//Reimann
		EdgeArray<double> lambdaDelay;   // sum(lambda_arcs)
		
		// Delay Scaling
		double delayScalingFactor;
		
		// Logical Effort
		EdgeArray<LogicalEffortPath> logicalEffortPath;
		EdgeArray<int> logicalEffortBacktrack; // direct index of the timing arc that propagated the logica effort
		double logicalEffortLoad;
		EdgeArray<double> logicalEffortArrivalTime;  // logical effort arrival time at this net
		
		TimingNetState() {
			load = 0;
			slew.set(0,0);
			logicalEffortLoad = numeric_limits<double>::quiet_NaN();
			logicalEffortArrivalTime.set(numeric_limits<double>::quiet_NaN(),
					numeric_limits<double>::quiet_NaN());
			logicalEffortBacktrack.set(-1, -1);
		} // end constructor		
	};
	
	struct State {
		vector<TimingArcState> arcs;
		vector<TimingNetState> nets;
	};
	
	State timingStateCurrent;
	State timignStateStored;
	
	// Timing nets are sorted in ascendent order by their logical depth
	vector<TimingNet> timingNets;
	
	// Timing arcs are sorted by the logical depth of their sink nets (the net
	// which is driven by the arc).
	vector<TimingArc> timingArcs;
    
	// Stores the net sinks of nets. For a net n stores all net which are
	// connected via a timing arc where n is the driver of the timing arc. Note
	// that for a net n driving a flip-flop, the flip-flops sink net is not
	// a sink net of net n.
	vector<int> timingSinkNets;
	vector<int> timingSinkArcs;
	
	// [TODO] Consider storing inside timing net structure.
	vector<int> timingArcPointers; // timing arcs which may drive a net i
	vector<int> timingSinkArcPointers; // timing arcs which are driven by net i
	vector<int> timingSinkNetPointers; // nets which may be driven by net i

	// Stores the net drivers of nets. Do not include duplicates.
	vector<int> timingDriverNets;
	vector<int> timingDriverNetPointers;
	
	// Net names for net i of timingNet array.
	vector<string> timingNetName;
	
	// Pin names for arc i of timingArc array.
	vector<string> timingArcDriverPinName;
	vector<string> timingArcSinkPinName;
	
	vector<int> timingTailNets;
	vector<int> timingTailNetMultiplicities;
    
	int timingOffsetToCombinationArcs;
	int timingOffsetToExtraSequentialArcs;
	int timingOffsetToExtraPrimaryOutputArcs;
	
	vector<double> timingMaxLoad;
	
	// Stores the required time at each net. We left it out of TimingNet
	// structure as it is not necessary to SA, which is the most timing
	// consuming part of our tool.
	vector< EdgeArray<double> > timingRequiredTime;
	
	// Number of dummy nets. Use this as an offset to skip dummy nets as all
	// dummy nets have a index less than timingNumDummyNets. The clock net is
	// also considered as a dummy net. Clock net has index equal to
	// (timingNumDummyNets - 1).
	
	int timingNumDummyNets;
	int timingOffsetToLevelOneNets; // skip dummy and level=0 (primary inputs and ffs driven nets)
	
	vector<int> timingOffsetToNetLevel;
	
	
	// Stores the worst arrival time for both edges and the (tail) net with
	// this worst arrival time.
	EdgeArray<int>    timingWorstArrivalTimeNet;
	EdgeArray<double> timingWorstArrivalTime;
	
	// Timign violations.
	double timingViolationSlew;
	double timingTotalNegativeSlack;
	double timingTotalPositiveSlack;
	double timingTotalAbsoluteSlack;
	
	// Used for thread-safe slew calculation.
	vector<double> timingViolationSlewVector;
	
	int timingNumPathsWithNegativeSlack;
	
	// Span
	vector<Vcell *> timingSortedCellsBySpan;
	vector<int> timingSpan;
	vector<int> timingReverseSpan;
	
	// Lambda-Delay Sensitivities
	vector< EdgeArray<double> > timingArcLambdaDelaySensitivity;
	vector< EdgeArray<double> > timingNetLambdaDelaySensitivity;
	
	// Stores sets of pseudo-independent cells. When a cell is sized, the other
	// cells in the same pseudo-independent set will not have their timing
	// context (arrival time, input slews) changed significantly. Also the load
	// seen by the cells inside a same set will never be affected by changes
	// in other cells in the same set.
	
	int timingNumPseudoIndependentSets;
	vector<int> timingPseudIndependentSets;
	vector<int> timingPseudIndependentSetPointers;
	
	// Stores locals nets of a net. Local nets include the net itself. Also
	// local nets are stored in topological order.
	vector<int> timingLocalNets;
	vector<int> timingLocalNetPointers;
	
	vector<int> timingLocalNetsIncludingSideNets;
	vector<int> timingLocalNetPointersIncludingSideNets;	
	
	vector<int> timingLocalArcs;
	vector<int> timingLocalArcPointers;

	vector<int> timingSideArcs;
	vector<int> timingSideArcPointers;

	vector<int> timingSideNets;
	vector<int> timingSideNetPointers;	
	
	vector<int> timingForwardNets;
	vector<int> timingForwardNetPointers;
	
	// Count the number of critical path passing through each net.
	vector<int> timingCriticalPathCounter;
	vector< pair<int,int> > timingNetSortedByCriticalPathCounter;
	
	// Previous net slews.
	vector< EdgeArray<double> > timingPreviousNetSlew;
	
	void buildTimingStructure();
	void updateTiming_Net( const int n, const int threadId = 0 );
	void updateTiming_Nets( const int n0, const int n1, const int threadId = 0 ); // [n0,n1)
	void updateTiming_WorstArrivalTime();
	void updateTiming_SlewViolation();
	void updateTiming_Deprecated();
	void updateTiming_Debug();
    
	void updateTiming_NetEndpoints(const int n);
	
	// Update only relevant timing information on driver nets and the net n
	// itself. Note that this method does not keep the whole timing information
	// consistent.
	void updateTimingLocally(const int n);
	void updateTimingLocallyIncludingSideNets(const int n);
	
	double computeSizingEffectOnDriverCellDelay(const int netIndex);
    double computeSizingEffectOnDriverCellFanoutDelay(const int netIndex);
    double updateCellAndComputeSizingEffectOnDriverCellDelay(Vcell * cell, const int newSize);
    void updateTimingDriverCell(const int n);
    
    // Refer to Li Li paper.
	double computeSizingEffectOnLambdaDelay(const int n);
	double computeSizingEffectOnLambdaDelaySensitivities(const int n);
    double computeSizingEffectOnDelayWithoutLambda(const int n);
	
	// [TODO] Explain them :)
	double computeSizingEffectOnLocalNegativeSlack(const int n);
	double computeSizingEffectOnAbsoluteSlack(const int n);
	double computeSizingEffectOnExpectedArrivalTime(const int n);
	double computeSizingEffectOnSlack(const int n);
    double computeSizingEffectOnSlackLiLi(const int n);
	double computeSizingEffectOnArcSaturatedSlack(const int n);
	double computeSizingEffectOnDepthSlack(const int n);
	double computeSizingEffectOnArrivalTimeVariance(const int n);
	
	double computeNorm();
	double computeStepSize(const double UB, const double L);
	double computeLagrange();
	
	//--------------------------------------------------------------------------
    //SA iterators
	vector<Vcell*>::iterator cellIterator;
	vector<Vcell*>::iterator criticalIterator;
	//--------------------------------------------------------------------------
    
    // LR iteration index
    double kIndex;
    double referenceLeakage;
    // Pointers to the respective cost/size vector
    vector<int> lagrangianMapPointers;
	// The cost/size vector. Stores the information needed by DP to traverse
    // the circuit minimizing the LR function. For each net (cell) stores the
    // accumulated cost and the corresponding size of previous cells. The second
    // index represents the cell size index
    vector< vector<double> > lagrangianMap;

	double sensitivityOffsetInputSlew;
	double sensitivityOffsetOutputLoad;
    
	// Used to compare to doubles.
	static const double EPSILON;
	
	SDCInfo sdcInfos;
	vector< LibParserCellInfo > lib_cells;
	OrgCells orgCells;
	
	vector< string > inputs;
	vector< pair<string,Vcell*> > outputs;
	vector< pair<string,string> > bestSolution;
	
	set< Wire > wires;
	//set< AddrCell > icells_addr;
	map< string, Net > nets;
	
	Vcell graphRoot;
	
	bool timingOk;
	
	multimap <double, pair<int,EdgeType> > arrivalTimingNet;
	
	// Stores the critical path computed by updateCriticalPath() method. The
	// cells are stored from the head (not including sequential/input driver) to
	// path tail (not including the sequential element);
	vector<Vcell *> criticalPath;
	vector<Vcell *> criticalPathEnlarged;
	
	EdgeType criticalPathEdgeType;
	double criticalPathArrivalTime;
	
	// Stores all combinational cells at top most critical paths computed by
	// updateTopCriticalPaths() method. Note that some cell may appear more
	// than once as they can participate in more than one top critical path.
	// However this may be a good side-effect when picking up randomly a cell
	// from this vector since a cell which appears multiple times is likely to
	// be more critical than others and also it has more chance to be chosen
	// randomly.
	vector<Vcell *> topCriticalCells;
	
	// All combinational cells which at least one input is being driven by a
	// primary input or a sequential element.
	vector<Vcell *> pathHeads;
    
	// All combinational cells which drive a primary output or a sequential
	// element.
	vector<Vcell *> pathTails;
    
	// Stores sizing solution.
	vector<int> storedSolution;
	vector<int> storedPastSolution;
	
	// Assign logical effort for each cell type
	vector< pair<string, double> > footprintLEPairs; //footprint and logical effort
    
	// Assign minimum slew for each cell type
	vector< pair<string, double> > footprintMinimumSlewPairs; //footprint and logical effort
	
	// Cload order
	vector<Vcell *> cloadOrdered;
	
	// Max circuit depth.
	int maxLogicalDepth;
	int maxReverseLogicalDepth;
    
	// Average number of sinks cells drive.
	double avgNumberOfSinks;
	
	//Alpha value used to calculate the slew target
	double alpha;
	
	void readVerilog();
	void readLib();
	void sortCells();
	void readSDC();
	void readSPEF();
    
    void setInitialLambda();
    void setInitialLambdaKKT();
    
	long unsigned int accepts, rejects;
	long unsigned int acceptsChanges, rejectsChanges;
	
#ifdef PARALLEL
	
	class Worker : public Poco::Runnable {
	private:
		Circuit *circuit;
		int n0;
		int n1;
		int threadId;
		
	public:
		void setup( Circuit  * circuit, const int n0, const int n1, const int threadId) {
			this->circuit = circuit;
			this->n0 = n0;
			this->n1 = n1;
			this->threadId = threadId;
		} // end method
		
		virtual void run() {
			circuit->updateTiming_Nets(n0, n1, threadId);
		} // end method
	}; // end class
	
	Poco::ThreadPool myThreadPool;
	
	void setupMultithreading(const int numThreads = 0);
	
	int threadNumThreads;
	
	vector< vector<int> > threadNetPointers;
	vector< vector<Worker> > threadWorkers;
	
	
#endif
	
	//--------------------------------------------------------------------------
	// Logical Effort
	//--------------------------------------------------------------------------	

	// Used to compare to doubles.
	static const double defaultH;	
	static const double defaultCin;
	
	EdgeArray<double> referenceSlew;
	vector<double> gains;
	
	// Stores the worst arrival time for both edges and the (tail) net with
	// this worst arrival time.
	EdgeArray<int>    logicalEffortWorstArrivalTimeNet;
	EdgeArray<double> logicalEffortWorstArrivalTime;
	
	// Timing violations.
	int logicalEffortNumPathsWithNegativeSlack;
	double logicalEffortTotalNegativeSlack;
	double logicalEffortTotalPositiveSlack;
	double logicalEffortTotalAbsoluteSlack;

	// D = k(GBH)^(1/k) + P
	// where
	// product of logical efforts   : G = g1*g2*...*gk
	// product of branching efforts : B = b1*b2*...*bk
	// sum of parasitic delays      : P = p1+p2+...+pk
	// path gain                    : H = Cout/Cin
	double computeLogicalEffortDelay(
			const int k, 
			const double G, 
			const double B, 
			const double P, 
			const double H) const;
		
	void logicalEffort_UpdateNetLoad(const int n);
	void logicalEffort_Reset_Cin(const int k);
	void logicalEffort_Reset();
	void logicalEffort_Propagate();
	void logicalEffort_Sizing();
	void logicalEffort_UpdateTiming_Wns();
	void logicalEffort_UpdateTiming_Net_Endpoint(const int n);
	void logicalEffort_UpdateTiming_Net(const int n);
	void logicalEffort_UpdateTiming();
	
	void logicalEffort_Report_Delays();
	void logicalEffort_Report_ArrivaTimes();
	
	void logicalEffort_Discretize_NetDriver(const int n);
	void logicalEffort_Discretize();
	
public:
	void logicalEffort_Test();
private:
	
	double logicalEffort_GetWns() const { return getT() -
			logicalEffortWorstArrivalTime.getMax(); }
	
	double logicalEffort_GetGain(const TimingNetState &netstate, const TimingArcState &arcstate) const {
		return netstate.logicalEffortLoad / 
				(arcstate.logicalEffortCin > 0? arcstate.logicalEffortCin : defaultCin);
	} // end method
	
	EdgeArray<double> logicalEffort_MeanMSE(const EdgeArray<double> SSE) {
		// see: http://www.stat.purdue.edu/~xuanyaoh/stat350/xyApr6Lec26.pdf
		// SSE: sum of squared errors
		const int N = gains.size();
		return sqrt(SSE/(N-2));
	}
	
	//--------------------------------------------------------------------------
	// Delay Scaling
	//--------------------------------------------------------------------------	
	
	int  updateDelayScalingFactors_PathTrace(const int nstart);
	void updateDelayScalingFactors();
	
public:
	string rootDir, benchName;
	double maxLeakage, maxTimingViol, maxWorstSlack;
	
	// Stores cell pointers not including input drivers in no specific order.
	vector< Vcell* > icells;
	
	// Stores cell pointers sorted in ascending order by their logical depth
	// including input. As input drivers and flip-flops have zero logical
	// depth they come at the beginning of the vector. By definition input
	// drivers are placed first than flip-flops.
	vector< Vcell* > depthSortedCells;
	int offsetSequential;
	int offsetCombinational;
	
	double maxTransition, worstSlack, worstSlew, worstDelay, worstNSDelay, avgDelay, minSlack, minSlew;
	bool hasViol;
	double timingViol, slewViol, minViol, loadViol, totalLeakage, bestCost;
	
	Vcell * worstDelayCell;
	
	void readInputFiles();
	
	// Creates cells vector ordered by Cload
	void createCloadVector();
	
	//Create input slew estimates
	void createInputSlewEstimates(Vcell* cell);
	
	// Update timing violation of the whole circuit.
	void updateTimingViol();
	
	// Stores on memory the current sizing solution.
	void storeSolution();
	
	// Stores on memory the past sizing solution.
	void storePastSolution();
	
	//Restiores from memory the stored sizing solution.
	void restoreFirstSolution();
	
	// Restiores from memory the stored sizing solution.
	void restoreSolution();
	
	// Sorts cells by higher load capacitance
	void cloadOrder();
	
	// Save current solution to .sizes file
	void saveSizes();
	void saveSizesDetailed(); // include size index and vth index
	void saveSizeTable();
	
	// Save best solution to .sizes file
	void saveBestSizes();
	
	// Store solution sizes in bestSolution vector
	void storeBestSolution();
	
	void saveSolution(const string &name);
	void readSolution(const string &name);
	
	void saveHistogramPathSlacks(const int iteration);
	void saveHistogramLambdas(const int iteration);
	
	/*
	 // Save .timing file (different from PT .timing file)
	 void saveTiming();
	 
	 // Save .mytiming file (different from PT .timing file)
	 void saveMyTiming();
	 
	 // Save .myPTtiming file (equal to PT .timing file)
	 void savePTTiming();
	 
	 // Read .timing file (no changes on statistics)
	 void readTiming();
	 
	 // Read .timing file (changes statistics)
	 void readTiming1();
	 
	 // call PrimeTime
	 void callPT();
	 */
	// Copy cells names and sizes to vector
	vector<pair<string,string> > copySizes();
	
	// Tiago's method
	int changeCellsReimann(const double temp);
	int changeCellsFastTimingSA(const double temp);
    void changeCellsFastLoadSA(const double temp);
	void changeCellsFastLoadSAUsedUFSC(const double temp);
	void changeCellsFastLeakageSA(const double temp);
	void criticalImprovement(const int temp);
    
    
    void loadSizes();
    void printCells();
	void lowTemperatureAnneal();
	void lowTemperatureAnneal2();
	void lowTemperatureAnnealTestedWithUFSC();
	void anneal();
    //LR
    void LR();
    void LRDPInitialization();
    void LR_PT();
    bool ignoreArc(TimingArc &arc, TimingArcState &arcstate);
	bool ignoreArcOzdal(const int arcIndex, TimingArcState &arcstate);
	void criticalTreeExtraction();
	void runDP();
	double evaluateLRS();
	double bestLowerBound;
    vector< EdgeArray<double> > deltaDelay_deltaLoad,deltaSlew_deltaLoad,deltaDelay_deltaSlew;
    double slewImpact(TimingArc &arc);
	double delayImpact(const int netIndex);
	bool isRoot(const int netIndex);
    void backTrack(const int netIndex);
    void backTrack2(const int netIndex);
    bool hasNegativeSlack(const int netIndex);
    
	
	// Build circuit graph and set initial solution (false = min leakage, true = random solution)
	void initialSolution(bool random);
	void initialSolutionNoChange();
	void printTree(Vcell* root);
	
	// Try to solve all load violations
	void solveLoadViol();
	
	// Update loads for whole circuit
	void calcLoads();
	
	// Update totalLeakage
	void calcLeakage();
	
	// Update loadViol
	void calcLoadViol();
	
	// Update slewViol
	void calcSlewViol();
	
	// Update loads for a changed cell (update previous cells loads)
	void updateLoads(Vcell * changedCell);
	
	// Set new cell size for a overloaded cell
	void setNewCellSize(Vcell* tmpCell);
	
	void changeCellsFlach1();
	void changeCellsFlach2();
	void changeCellsFlach3();
	
	void changeCellsGraci();
	
	void updateCellLoad(Vcell * cell);
	void updateCellTiming(Vcell * cell);
	void updateCellType(Vcell * cell, int typeIndex );
	void updateTopCriticalTailNets();
	
	// [TODO] Explain...
	bool updateCellTypeByLiLi(Vcell * cell);
    bool reducingLeakageByLiLi();
	
	// [TODO] Explain...
	bool updateCellTypeLagrangeRelaxation(Vcell * cell, const double alpha = 1.0);
	bool updateCellTypeLagrangeRelaxationSensitivities(Vcell * cell, const double gamma, const double alpha = 1.0);
    bool updateCellTypeLagrangeRelaxationSensitivitiesTestingChanges(Vcell * cell, const double gamma, const double alpha = 1.0);
	bool updateCellTypeConsideringLocalSlack(Vcell * cell, const double alpha = 1.0);
	bool updateCellTypeConsideringAbsoluteSlack(Vcell * cell);
	bool updateCellTypeToReducePowerConsideringLocalSlack(Vcell * cell);
	
	// Try to resize the cell such that the cell's sink net arrival time best
	// match the expected arrival time at that net.
	bool updateCellTypeToMeetExpectedArrivalTime(Vcell * cell);
	
	// Try to resize the cell such that the cell's sink net arrival time best
	// match the required time at that net.
	bool updateCellTypeToMeetRequiredTime(Vcell * cell);
	
	// Helper methods to size a cell.
	bool upsize(Vcell * cell, const bool allowLoadViolation = false);
	bool downsize(Vcell * cell, const bool allowLoadViolation = false);
	bool increaseVth(Vcell * cell, const bool allowLoadViolation = false);
	bool decreaseVth(Vcell * cell, const bool allowLoadViolation = false);
	
	// Helper methods to size neighbourhood cells.
	void downsizeSinks(const int n);
	void downsizeSiblings(const int n);
	void upsizeDrivers(const int n);
	
	// Print current timing information.
	void printTiming(const string &label = "");
	void printTimingDigest(const bool header = false, const string &title = "");
	
	void printDelaysVersusExpectedDelays();
	
	// Print current timing information.
	void printPathTails();
	
	// Print critical path timing.
	void printCriticalPathTimingReport(ostream &out);
	
	// Print some library stats.
	void printLibraryReport(ostream &out);
	
	// Generate a sight visualization file.
	void printSigth(const string &filename);

	// Dump vt stats for machine learning training.
	void dumpVtStats(const std::string &filename);

	// Compute cell depths.
	void computeCellDepths();
	void computeCellReverseDepths();
	
	// Compute pseudo-independent set.
	void computePseudoIndependentSets();
	
	// Compute spans.
	void computeSpans();
	void computeReverseSpans();
	
	// Compute local nets.
	void computeLocalNets();
	void computeLocalArcs();
	
	// Compute combinational clouds.
	void computeClouds();
	
	// Compute expected delays.
	void computeExpectedDelays();
	
	// Compute cell depths.
	void computePathBoundaries();
	void computePathBoundaries_CheckCell(Vcell * cell);
    
	// Compute the average number of sinks a cell drives.
	void computeAvgNumberOfSinks();
	
	// Create a helper structure to donw and up size cell as well to increase
	// and decreate Vth.
	void computeCellSizingOptions();
	
	// Compute the number of paths passing through nets.
	void computeNumberOfPathsPassingTrhuNets();
	void computeWeightedNumberOfPathsPassingTrhuNets();
	
	// Compute logical effort of timing arcs.
	void computeLogicalEffort_ReferenceSlew();
	void computeLogicalEffort_LinearLeastSquares(
			const vector<double> &x, 
			const vector<double> &y, 
			double &a, 
			double &b);
	void computeLogicalEffort_LinearLeastSquaresError(
			const vector<double> &x, 
			const vector<double> &y, 
			const double a, 
			const double b,
			double &residuum);
	void computeLogicalEffort_TimingArc(
			const LibParserLUT &lut, 
			const double Cin, 
			const double inputSlew, 
			const vector<double> &gains,
			double &g, 
			double &p,
			double &residuum);
	void computeLogicalEffort_Report(
			ostream &out);
	void computeLogicalEffort();
		
	// Find the current critical path.
	void updateCriticalPath();
    
	// Find the current critical path.
	void updateTopCriticalPaths(const int numPaths);
	
	// Find the current non-critical paths.
	void updateLeastCriticalPaths(const int numPaths);
	
	// Calculate the number of critical path passing through nets.
	void updateCriticalPathCounter();
	
	//
	void updateTopLoadPaths();
	
	double calculateInnerSlack();
	
	// Return a color representing a temperature. Parameter weight must be in
	// the [0,1] range.
	void colorTemperature( const double weight, int &r, int &g, int &b ) const;
	
	// Randomly choose a combinational cell.
	Vcell * chooseCellRandomlyFromVector(vector<Vcell *> &v);
	
	// Compute fanout/fanin ratio of a cell.
	double computeFanoutFaninRatio(const LibParserCellInfo &cellinfo, const double load);
	double computeFanoutFaninRatio(Vcell * cell);
	double computeFanoutFaninRatioLE(const LibParserCellInfo &cellinfo, double load);
	double computeFanoutFaninRatioLE(Vcell * cell);
	double computeFanoutFaninRatioLEPosContest(const LibParserCellInfo &cellinfo, double load);
	double computeFanoutFaninRatioLEPosContest(Vcell * cell);
	
	// Compute the slew violation caused by a slew transition and the number of
	// pins (fanout) driven such slew transition.
	double computeSlewViolation( const double slew, const int fanout ) const;
	
	// Compute the average/sum input capacitance of a cell. [TODO] These method
	// return a constant value for a given cell, so store this value at some
	// place instead recalculating it every time.
	double computeInputCapacitance(const LibParserCellInfo &cellinfo);
	double computeAvgInputCapacitance(const LibParserCellInfo &cellinfo);
	double computeSumInputCapacitance(const LibParserCellInfo &cellinfo);
	
	//To assign the minimum slew for each type of the cell
	double assignMinimumSlew(int footprintIndex );
	
	// Check load violation on driving (previous) cels.
	bool hasLoadViolationOnPreviousCells( Vcell * cell );
	
	// [TODO] Explain it.
	void findCriticalBest();
	
	bool isClockNet(const int n) const { return n == timingNumDummyNets - 1; }
	
	// -------------------------------------------------------------------------
	// Sizing Methods
	// -------------------------------------------------------------------------
	
    // Simply choose a random size for every cell.
	void sizingRandom();
	
	// Randomly choose cells and set its size to the one which best fits the
	// fanout-of-4 rule.
	void sizingRandomFo4();
    
	// Sizes cells based on fanout-of-n, walking from cell with largest logical
	// depth to lowest depth.
	void sizingDepthFanout();
	void sizingDepthFanoutXLogicalEffort();
	void sizingSlewTarget();
	void sizingDepthFanoutXLogicalEffortCopy();
	void sizingDepthFanoutXLogicalEffortPosContest();
	void sizingGreedy();
    
	void sizingByLiLi();
	void sizingLagrangeRelaxation(const double resetLambdas = true);
	void sizingLagrangeRelaxationSensitivities(const double resetLambdas = true);
    void sizingLagrangeRelaxationTestingChanges(const double resetLambdas = true);
	void sizingLambdaGreedyStepSize();
    
	// Tries to solve load and slew violations selecting less leakage cells.
	void sizingForNoLoadAndSlewViolationByLiLi();
	
	// [TODO] Explain :)
	void sizingForNoLoadAndSlewViolationByLivramento();
	
	// Tries to solve load violations without taking into account any aspect of
	// the current solution.
	void sizingForNoLoadViolation();
	
	//
	void sizingDecreaseVthOfHighSpanCells();

	// Decreases vth of cells that drives several paths.
	void sizingDescreaseVthOfBottleneckCells();
	
	// [TODO] Explain it.
	void sizingProportionalDelay();
	
	// [TODO] Explain it.
	void sizingVaiVem1();
	
	// Simulated annealing
	void simulatedAnnealing(const long double temp, bool critical);
    void simulatedAnnealingTestedWithUFSC(const long double temp, bool critical);
    
	// [TODO] Explain it.
	void sizingDepthFanoutX();
	
	//sizing using Logical Effort
	void sizingLogicalEffort();
    
	// The idea was to change the VTh of cell on the critical path. Worked, but
	// not so well.
	void sizingCriticalPathVTh();
	
	// [TODO] Explain it.
	void sizingCriticalPathSensitivity();
	
	// [TODO] Explain it.
	void sizingGreedyPower();
	
	// [TODO]
	void sizingToMeetExpectedArrivalTimes();
	void sizingToMeetRequiredTimes();
	
	void timingRecovery();
	void timingRecoveryPathCounter();
	void powerRecovery();
	void powerRecoveryByDecreasingVth();
	
	
	// -------------------------------------------------------------------------
	// Walker Methods
	// -------------------------------------------------------------------------
	// Allow one to iterate over circuit nodes in some order. Each time a node
	// is visited the callback function parameter is called to such node.
	
	typedef void (Circuit::*WalkerCallback) ( Vcell * );
	
	// Walk from outputs to inputs in a BFS manner.
	void walkBackward( WalkerCallback callback );
    
	// Walk from inputs to outputs in a BFS manner.
	void walkForward( WalkerCallback callback );
    
	// Walk from seed to a temporal barrier (inputs, sequential elements)
	// following cells with worst slack.
	void walkThroughWorstSlackPath( Vcell * seed, EdgeType edgeType, WalkerCallback callback );
	
	// Walk over critical path. Don't forget to call updateCriticalPath() before
	// walking.
	void walkCriticalPathForward( WalkerCallback callback );
	void walkCriticalPathBackward( WalkerCallback callback );
	
	// Walk from primary output and flip-flops to temporal barriers (or vice-versa).
	void walkBackwardFromTemporalBarriers( WalkerCallback callback );
	void walkForwardFromTemporalBarriers( WalkerCallback callback );
	
	// Walk following logical depths.
	void walkBackwardDepth( WalkerCallback callback );
	void walkForwardDepth( WalkerCallback callback );
	
	// Walk from largest load capacitance to smallest load capacitance and smallest to largest
	void walkCloadLargeToSmall( WalkerCallback callback );
	void walkCloadSmallToLarge( WalkerCallback callback );
	
	// Walk over top critical cells
	void walkTopCriticalCells( WalkerCallback callback ) ;
	
	// Walk following backtrack pointers. Backtrack pointers point to the input
	// net of a cell which dominate (worst arrival time + arc delay) other
	// inputs.
	template< class Stepper>
	void walkFollowingBacktrackPointers( const int netIndex, const EdgeType edgeType, Stepper &stepper );
	void walkFollowingBacktrackPointers( const int netIndex, const EdgeType edgeType, WalkerCallback callback );
	
	// -------------------------------------------------------------------------
	// Stepper Methods
	// -------------------------------------------------------------------------
	// Methods intended to be used as callback in walker methods.
	// [TODO] Replace stepper methods with stepper functors as functors are
	//        more flexible.
	
	void stepperPrint(Vcell *cell);
	void stepperFo4(Vcell *cell);
	void stepperFoX(Vcell *cell);
	void stepperSlewTarget(Vcell *cell);
	void stepperFoXLogicalEffort(Vcell *cell);
	void stepperFoXLogicalEffortCopy(Vcell *cell);
	void stepperFoXLogicalEffortPosContest(Vcell *cell);
	void stepperCriticalPathUpdater(Vcell *cell);
	void stepperNoLoadViolation(Vcell *cell);
	void stepperDepthProportionalDelay(Vcell *cell);
	void stepperVaiVem1(Vcell* cell);
	void stepperVaiVem2(Vcell* cell);
	void stepperSetSlewTarget(Vcell* cell);
	void stepperSetLocalCriticality(Vcell* cell);
    void stepperSetLocalCriticalityGC(Vcell* cell);
	void stepperSetCellSlewTarget(Vcell* cell);
    void stepperRefiningSlewTargets(Vcell* cell);
	void stepperRefineSlewTargets(Vcell* cell);
	void stepperSetGlobalCriticality(Vcell* cell);
	void stepperGreedySizing(Vcell* cell);
    void stepperForNoLoadAndSlewViolationByLiLi(Vcell* cell);
	void stepperForNoLoadAndSlewViolationByLivramento(Vcell *cell);
    void stepperReducingLeakageByLiLi(Vcell* cell);
	
	struct StepperPrint { void operator()(Vcell *cell); };
	
	struct StepperPathLogicalEffort {
		Circuit * super;
		Vcell * previousCell;
		Vcell * currentCell;
		double pathLogicalEffort;
		double pathBranchingEffort;
		double pathFanoutEffort;
		int numPathCells;
		
		StepperPathLogicalEffort(Circuit * const circuit) : super(circuit) {
			pathLogicalEffort = 1;
			pathBranchingEffort = 1;
			pathFanoutEffort = 1;
			currentCell = NULL;
			previousCell = NULL;
			numPathCells = 0;
		} // end constructor
		
		void operator()(int arcIndex);
	};
	
    struct StepperAssignSizeLogicalEffort {
		Circuit * super;
        double gainPerStage;
		double fanoutPerStage;
		double fanoutEffortPerStage;
		int numPathCells;
		
		StepperAssignSizeLogicalEffort (Circuit * const circuit) : super(circuit) {
			gainPerStage = 0;
			fanoutPerStage = 0;
			numPathCells = 0;
			fanoutEffortPerStage = 0;
		} // end constructor
		
		void operator()(int arcIndex);
	};
	
    struct StepperAssignSizeLogicalEffortFanout {
		Circuit * super;
		double fanoutPerStage;
		double fanoutEffortPerStage;
		
		StepperAssignSizeLogicalEffortFanout (Circuit * const circuit) : super(circuit) {
			fanoutPerStage = 0;
			fanoutEffortPerStage = 0;
		} // end constructor
		
		void operator()(int arcIndex);
	};
	
    
	// -------------------------------------------------------------------------
	// Look-Up Table
	// -------------------------------------------------------------------------
	
	// [TODO] Put these methods in a separated class, maybe a static one.
    
	static double lookup( const LibParserLUT &lut, const double x, const double y );
    
	static double lookupDelay( const LibParserTimingInfo &timingInfo, const EdgeType edgeType, const double inputSlew, const double loadCapacitance );
	static double lookupOutputSlew( const LibParserTimingInfo &timingInfo, const EdgeType edgeType, const double inputSlew, const double loadCapacitance );
    
	// -------------------------------------------------------------------------
	// Timing
	// -------------------------------------------------------------------------
    
	// Updates the whole circuit timing. Call this method at startup or when
	// more than one cell has been changed at same time.
	void updateTiming();
	void updateTimingSingleThreaded();
	void updateTimingMultiThreaded();
	
	void updateTimingLR();
	void updateTimingLR_KKT();
	void updateTimingLR_PT();
	void updateTimingLR_PT_KKT();
	
	void updateLambdasByLiLi();
	void updateLambdasByTennakoon();
	void updateLambdasByTennakoonFig4(const double stepSize);
	void updateLambdasByFlach();
	void updateLambdasByFlachNew();
	void updateLambdasCriticality();
	
    void updateLambdas_Subgradient();
	void updateLambdas_Normalization();
	void updateLambdas_KKT();
	void updateLambdas_Nets();
	
    void updateLambdasByOzdal();
    
	void resetTimingArcLambdas( const double value );
	
	void computeLambdaSTA();
	
    void callPT();
	void readTimingLR();
	
	// Updates the circuit timing when ONLY ONE cell has been changed. Note that
	// this method KEEP the whole circuit TIMING CORRECT, but it is pretty much
	// efficient than updating the whole circuit timing as it propagates only
	// the changes caused by the changed cell. Also known as cone timing.
	void updateTiming( Vcell * cell);
	
	// Updates required times (slacks) of nets. Must be called after the circuit
	// timing (arrival times) has been computed. Note that calling this method
	// is not necessary if only slacks at path tails are required.
	void updateRequiredTime();
	// Also updates lambdas used in Lagrangian Relaxation
	void updateRequiredTimeLR();
	void updateRequiredTimeLR_KKT();
	// Updates lambdas used in Lagrangian Relaxation after reading slacks from PT
	void updateRequiredTimeLR_PT();
	void updateRequiredTimeLR_PT_KKT();
	
	// Lambda-Delay Sensitivities.
	void updateLambdaDelaySensitivities();
	
	// Computes the cell timing information based on current cell context (input
	// slews and output load).
	static void computeArcTiming(
		const LibParserTimingInfo &timingInfo,
		const EdgeArray<double> inputSlew,
		const double load,
		EdgeArray<double> &outDelay,
		EdgeArray<double> &outSlew
	);
    
	// Compute delay and slew sensitivities of timing arc k. Note that 
	// sensitivities are calculated based on the current arc's cell context:
	// output load, input slew.
	void computeArcTimingDelayAndSlewSensitivityToOutputLoad( const int k, const int size,  EdgeArray<double> &sensitivityDelay, EdgeArray<double> &sensitivitySlew );
	void computeArcTimingDelayAndSlewSensitivityToInputSlew( const int k, const int size, EdgeArray<double> &sensitivityDelay, EdgeArray<double> &sensitivitySlew );
	
	// -------------------------------------------------------------------------
	// Auxiliaries
	// -------------------------------------------------------------------------
    
	double myRound(const double v, const int digits) const {
		const double m = pow(10.0,digits);
		return floor(v*m + 0.5)/m;
	} // end method
	
	// Checks if two numbers (double) are approximately equal.
	// Source: http://java2s.com/Tutorial/Cpp/0040__Data-Types/Testswhethertwofloatingpointnumbersareapproximatelyequal.htm
	static bool nearlyEqual( const double x, const double y, const double precision = EPSILON ) { 
		if (x == 0) return fabs(y) <= precision;
		if (y == 0) return fabs(x) <= precision;
		return fabs(x - y) / max(fabs(x), fabs(y)) <= precision;		
	} // end method
	
	// Checks if a number is approximately zero.
	static bool nearZero( const double v, const double precision = EPSILON  ) { return fabs(v) < precision; }
	
	// Source: The Art of Computer Programming by Knuth
	// Use this functions instead of nearlyEqual() and nearlyZero() as they are
	// more precise and efficient. The other functions are kept in order to not
	// change the results w.r.t. the ISVLSI 2013 paper.
	
	template<typename T>
	static bool approximatelyEqual(const T a, const T b, const T precision = EPSILON) {
		return std::abs(a - b) <= 
				((std::abs(a) < std::abs(b) ? std::abs(b) : std::abs(a)) * precision);
	} // end method	

	template<typename T>
	static bool approximatelyZero(const T a, const T precision = EPSILON) {
		return std::abs(a) <= 
				((std::abs(a) < 0 ? 0 : std::abs(a) * precision));
	} // end method	
	
	template<typename T>
	static bool definitelyGreaterThan(const T a, const T b, const T precision = 1e-6) {
		return (a - b) > 
				((std::abs(a) < std::abs(b) ? std::abs(b) : std::abs(a)) * precision);
	} // end method

	template<typename T>
	static bool definitelyLessThan(const T a, const T b, const T precision = 1e-6) {
		return (b - a) > 
				((std::abs(a) < std::abs(b) ? std::abs(b) : std::abs(a)) * precision);
	} // end method	
	
	// Assign logical effort for each cell type.
	// [WARNING][TODO] This method is hard coded! Only works with the ISPD 2012
	// library.
	void assignLogicalEffort();
	void assignLogicalEffortContestLibrary();
	void assignMinimumSlew();
	
	// -------------------------------------------------------------------------
	// Get Methods
	// -------------------------------------------------------------------------
	
	const TimingArcState &getTimingArcState(const int k) const { return timingStateCurrent.arcs[k]; }
	const TimingNetState &getTimingNetState(const int n) const { return timingStateCurrent.nets[n]; } 	
	
	const TimingArcState &getTimingArcState(const int k, const State &state) const { return state.arcs[k]; }
	const TimingNetState &getTimingNetState(const int n, const State &state) const { return state.nets[n]; } 	
	
	TimingArcState &getTimingArcState(const int k) { return timingStateCurrent.arcs[k]; }
	TimingNetState &getTimingNetState(const int n) { return timingStateCurrent.nets[n]; } 	
	
	TimingArcState &getTimingArcState(const int k, State &state) { return state.arcs[k]; }
	TimingNetState &getTimingNetState(const int n, State &state) { return state.nets[n]; } 		
	
	double getAvgNumberOfSinks() const { return avgNumberOfSinks; }
	double getClkPeriod() const { return this->sdcInfos.clk_period; }
	int getPathTailsSize() const { return this->pathTails.size(); }
	int getSize() const { return this->icells.size(); }
	
	double getWorstSlack() const {
		return min(
                   this->sdcInfos.clk_period - timingWorstArrivalTime[RISE],
                   this->sdcInfos.clk_period - timingWorstArrivalTime[FALL]);
	} // end method
	
	double getSlackSlack(double slack) const {
		const double T = this->sdcInfos.clk_period;
		
		return (-(min(0.0, slack))/T + 1);
	} // end method

	double getTotalNegativeSlackSlack() const {
		const double T = this->sdcInfos.clk_period;
		
		return (-(min(0.0, -timingTotalNegativeSlack))/T + 1);
	} // end method
	
	int getNumCellSizes( Vcell * cell ) const {	
		return orgCells.oCells[cell->footprintIndex].cells.size();
	} // end method
		
	EdgeArray<double> getNetSlack(const int n) const { return timingRequiredTime[n] - getTimingNetState(n).arrivalTime; }
	EdgeArray<double> getNetNegativeSlack(const int n) const { return min(EdgeArray<double>(0,0), timingRequiredTime[n] - getTimingNetState(n).arrivalTime); }
	
	double getArcInputCapacitance(const int k, const int size) {
		const TimingArc &arc = timingArcs[k];
		const LibParserCellInfo &cellinfo = orgCells.oCells[arc.cell->footprintIndex].cells[size];	
		return cellinfo.pins[arc.pin].capacitance;
	} // end method
	
	EdgeArray<double> getNetArrivalTimeStandardDeviation(const int n) const {
		EdgeArray<double> x1(0,0);
		EdgeArray<double> x2(0,0);
		
		const int k0 = timingArcPointers[n];
		const int k1 = timingArcPointers[n+1];
		for ( int k = k0; k < k1; k++ ) {
			const TimingArcState &arc = getTimingArcState(k);
			
			const EdgeArray<double> arrivalTime = (arc.arrivalTime + arc.delay.getReversed());
			x1 += arrivalTime;
			x2 += pow(arrivalTime, 2.0);
		} // end for
		
		const int numTimingArcs = k1 - k0;
		x1 /= numTimingArcs;
		x2 /= numTimingArcs;
		
		return sqrt(x2 - pow(x1, 2.0));
	} // end method
	
	double getTimingViolation() const { return timingTotalNegativeSlack; }
	double getSlewViolation() const { return timingViolationSlew; }
	double getLoadViolation() const { return loadViol; }
	double getDeltaLoad(Vcell * cell, Vcell * nextCell, int newIndex);
	
	const vector< pair<string,string> > &getBestSolution() const { return bestSolution; }
	
	double getT() const { return sdcInfos.clk_period; };
	void setT(const double T) { sdcInfos.clk_period = T; }
	
	int getLibPinIndexByName(const LibParserCellInfo &lcell, const string &name) {
		for ( int i = 0; i < lcell.pins.size(); i++ )
			if ( lcell.pins[i].name == name )
				return i;
		return -1;
	} // end method
	
	double getLibPinCapacitance(const int footprint, const int type, const int pin) {
		return orgCells.oCells[footprint].cells[type].pins[pin].capacitance;
	} // end method	
	
	double getLibPinCapacitanceFromTimingArc(const TimingArc &arc) {
		return arc.cell->actualInstType->pins[arc.pin].capacitance;
	} // end method
	
	string getLibPinName(const TimingArc &arc) {
		return (arc.pin >= 0)?
				arc.cell->actualInstType->pins[arc.pin].name : "<null>";
	} // end method

	double getMaxLoad(Vcell * cell) {
		// [HARD CODED] Assuming the output pin is the first one.
		return cell->actualInstType->pins[0].maxCapacitance;
	} // end method

	double getMaxLoad(const int footprint, const int type) {
		// [HARD CODED] Assuming the output pin is the first one.
		return orgCells.oCells[footprint].cells[type].pins[0].maxCapacitance;
	} // end method	
	
	double getLoadViolation(Vcell * cell) {
		const double maxLoad = getMaxLoad(cell);
		return (cell->actualLoad > maxLoad)?
			(cell->actualLoad - maxLoad) : 0.0;
	} // end method

	double getLoadViolationForCellType(Vcell * cell, const int type) {
		const double maxLoad = getMaxLoad(cell->footprintIndex, type);
		return (cell->actualLoad > maxLoad)?
			(cell->actualLoad - maxLoad) : 0.0;
	} // end method
	
	// -------------------------------------------------------------------------
	// Constructor
	// -------------------------------------------------------------------------
	bool stopWalking;
	
	Circuit () :
	maxLeakage (DBL_MAX),
	maxTimingViol (DBL_MAX),
	maxWorstSlack (-DBL_MAX),
	maxTransition (0.0),
	worstSlack (DBL_MAX),
	worstSlew (0.0),
	worstDelay (0.0),
	worstNSDelay (0.0),
	avgDelay (0.0),
	minSlack (-DBL_MAX),
	minSlew (DBL_MAX),
	minViol (DBL_MAX),
	bestCost (DBL_MAX),
	timingOk (false),
	timingViolationSlew(0),
	timingTotalNegativeSlack(0),
	loadViol(0),
	slewViol(0),
	timingViol(0),
	accepts(0),
	rejects(0),
	acceptsChanges(0),
	rejectsChanges(0),
	kIndex(0)
	{
	
		// Default option values for baseline (ISVLSI 2013) flow.
		optIgnoreLeakagePower = false;	
		
		optEnableLambdaDelaySensitivities = true;
		optEnableSlackFiltering = true;
		optEnableGamma = true;
		optEnableLoadViolationFiltering = true;

		optUseTennakoon = false;
		
		optRunInitialSizing = true;
		optRunInitialSizingForLoadAndSlewViolationRemoval = true;
		optRunPowerRecovery = true;
		optRunTimingRecovery = true;	
		
	}
};

// -----------------------------------------------------------------------------

template<typename Stepper>
void Circuit::walkFollowingBacktrackPointers( const int netIndex, const EdgeType edgeType, Stepper &stepper ) {
	int n = netIndex;
	EdgeType e = edgeType;
	
	while (true) {
		const TimingNet &net = timingNets[n];
		if ( net.depth == 0 )
			break;

		TimingNetState &state = getTimingNetState(n);
		stepper(state.backtrack[e]);
		
		n = timingArcs[state.backtrack[e]].driver;
		e.reverse();
	} // end while
} // end method