/*
 *  Circuit.h
 *  sizer
 *
 *  Created by Tiago Reimann on 11/01/12.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 */

#include <cmath>
#include <climits>
#include <cstdlib>
#include <cfloat>
#include <cassert>

#include <vector>
#include <set>
#include <queue>
#include <stack>
#include <deque>
#include <bitset>
#include <string>
#include <iomanip>
using std::setw;
#include <sstream>
using std::ostringstream;
#include <limits>
using std::numeric_limits;

#include "global.h"

#include "Vcell.h"
#include "Stopwatch.h"
#include "EdgeArray.h"
#include "RCTree.h"

#ifdef REMOTE_PRIMETIME
#include "PracticalSocket.h"
#endif

#ifdef PARALLEL
#include <Poco/ThreadPool.h>
#include <Poco/Runnable.h>
#include <Poco/Environment.h>
#endif

using namespace std;

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
	
class Circuit {
	
	class RCTreeDriver {
	private:
		const LibParserTimingInfo &clsTimingInfo;
		const EdgeArray<double> &clsInputSlew;
	public:
		RCTreeDriver(const LibParserTimingInfo &timingInfo, const EdgeArray<double> &inputSlew ) :
			clsTimingInfo(timingInfo), clsInputSlew(inputSlew) {
		} // end constructor

		EdgeArray<double> computeSlew(const EdgeArray<double> &ceff) const {
			return EdgeArray<double>(
				lookup(clsTimingInfo.riseTransition, ceff[RISE]*1e15, clsInputSlew[FALL]) * 1e-12,
				lookup(clsTimingInfo.fallTransition, ceff[FALL]*1e15, clsInputSlew[RISE]) * 1e-12
			);
		} // end method
		
		EdgeArray<double> computeDelay(const EdgeArray<double> &ceff) const {
			return EdgeArray<double>(
				lookup(clsTimingInfo.riseDelay, ceff[RISE]*1e15, clsInputSlew[FALL]) * 1e-12,
				lookup(clsTimingInfo.fallDelay, ceff[FALL]*1e15, clsInputSlew[RISE]) * 1e-12
			);
		} // end method		
		
		EdgeArray<double> getInputSlew() const { return clsInputSlew * 1e-12; }
		
	}; // end class	
	
	// -------------------------------------------------------------------------
	
	class DigestDescriptor {
	private:
		const Circuit &clsCircuit;
		const string clsTitle;
				
		struct Column {
			string label;
			int width;
			const double * valueptr;
			
			Column() : label(""), valueptr(0), width(0) {}
			Column(const string &label, const double &value, const int width) :
				label(label),
				valueptr(&value),
				width(width) 
			{}
		};
		
		Stopwatch clsStopwatch;
		
		int clsIteration;
		int clsTotalExtraWidth;
		
		vector<Column> clsExtraColumnValues;
		map<string, int> clsExtraColumnMapping;
		
		void printHeader() {
			const int M = 80 + clsTotalExtraWidth;
			const int N = 12;
			
			cout << setw(M) << setfill('=') << "\n" << setfill(' ');

			cout << setw((M + clsTitle.length()) / 2) << clsTitle << "\n";
			cout << setw(M) << setfill('=') << "\n" << setfill(' ');

			cout
				<< setw(4) << "#"
				<< setw(N) << "Worst"
				<< setw(N) << "Leakage"
				<< setw(N) << "Timing"
				<< setw(N) << "Slew"
				<< setw(N) << "Load"
				<< setw(N) << "Path"
			;

			for ( int i = 0; i < clsExtraColumnValues.size(); i++ ) {
				const Column &c = clsExtraColumnValues[i];
				cout << setw(c.width) << c.label;
			} // end for
			cout << "\n";
			
			cout
				<< setw(4) << ""
				<< setw(N) << "Slack"
				<< setw(N) << "(W)"
				<< setw(N) << "Violation"
				<< setw(N) << "Violation"
				<< setw(N) << "Violation"
				<< setw(N) << "Violation"
				<< "\n";

			cout << setw(M) << setfill('=') << "\n" << setfill(' ');
		} // end method		
		
		double myRound(const double v, const int digits) const {
			const double m = pow(10.0,digits);
			return floor(v*m + 0.5)/m;
		} // end method

		
	public:	
		
		DigestDescriptor(const Circuit &circuit, const string &title) : 
			clsCircuit(circuit),
			clsTitle(title), 
			clsIteration(0), 
			clsTotalExtraWidth(0) 
		{
			clsStopwatch.start();
		}
		
		void print() {
			const int N = 12;
			
			if ( clsIteration == 0 ) {
				printHeader();
				cout << setw(4) << "-";
			} else {
				cout << setw(4) << clsIteration;
			} // end else
			
			const double roundedLoadViol = myRound(clsCircuit.getLoadViolationUsingEffectiveCap(), 2);
			const double roundedSlewViol = myRound(clsCircuit.getSlewViolation(), 2);

			const double roundedTimingViol = myRound(clsCircuit.getTimingViolation(), 2);
			const double roundedWorstSlack = myRound(clsCircuit.getWorstSlack(), 2);			

			cout
				<< setw(N) << roundedWorstSlack
				<< setw(N) << (clsCircuit.getLeakagePower() / 1e6)
				<< setw(N) << roundedTimingViol
				<< setw(N) << roundedSlewViol
				<< setw(N) << roundedLoadViol
				<< setw(N) << clsCircuit.getNumPathsWithNegativeSlack();

			for ( int i = 0; i < clsExtraColumnValues.size(); i++ ) {
				const Column &c = clsExtraColumnValues[i];
				cout << setw(c.width) << (*c.valueptr);
			} // end for

			//cout << "\t" << myRound(clsCircuit.getLoadViolationUsingEffectiveCap(), 2) << " / " << myRound(clsCircuit.getLoadViolationUsingDownstreamCap(), 2) << " / " << myRound(clsCircuit.getLoadViolationUsingEffectiveCapCellWise(), 2);

			cout << "\n";
			
			clsIteration++;
		} // end method
		
		void addExtraColumn(const string &label, const double &value, const int width = 12 ) {
			clsExtraColumnValues.push_back(Column(label, value, width));
			clsTotalExtraWidth += width;
		} // end method
		
		double getElapsedTime() {
			return clsStopwatch.getElapsedTime();
		} // end method
		
		~DigestDescriptor() {
			const double slack = clsCircuit.getWorstSlack();

			const int x0 = clsCircuit.getNumTailNets()*2;
			const double x1 = myRound(100 * (clsCircuit.getNumPathsWithNegativeSlack() / double(x0)), 2);
			const double T = clsCircuit.getT();

			const double roundedLoadViol = myRound(clsCircuit.getLoadViolationUsingEffectiveCap(), 2);
			const double roundedSlewViol = myRound(clsCircuit.getSlewViolation(), 2);

			const double roundedTimingViol = myRound(clsCircuit.getTimingViolation(), 2);
			const double roundedWorstSlack = myRound(clsCircuit.getWorstSlack(), 2);			

			cout << "\n";
			cout << "Summary (T=" << T << "): " << clsTitle << "\n";
			cout << "\tRuntime (h:mm:ss): " << clsStopwatch.getFormattedTime() << "\n";
			cout << "\tWorst Slack......: " << roundedWorstSlack << " [" << slack << "]" /*<< " (" << minSlack << ")"*/ << "\n";
			cout << "\tLeakage..........: " << (clsCircuit.getLeakagePower() / 1e6) << "\n";
			cout << "\tTiming Violation.: " << roundedTimingViol << " [" << clsCircuit.getTimingViolation() << "] #Paths: " << clsCircuit.getNumPathsWithNegativeSlack()  << "/" << (x0) << " (" << x1 << "%)" << "\n";
			cout << "\tSlew Violation...: " << roundedSlewViol << " [" << clsCircuit.getSlewViolation() << "]" << "\n";
			cout << "\tLoad Violation...: " << roundedLoadViol << "\n";
			cout << endl;	
			
			const_cast<Circuit&>(clsCircuit).primeTimeReport();
		} // end descrutctor
		
	}; // end class
	
	// -------------------------------------------------------------------------
	
	// [NOTE] HARD CODED: Number of Vth = 3.
	struct CellSizingOption {
		// Vth = 0 --> higher Vth, slowest
		// Vth = n --> lower Vth, fastest
		
		// [Vth][Size]
		vector< vector<int> > option;
		
		// Map instance type index to option:
		// mapping[inst type index] -> (Vth index, Size Index)
		vector< pair<int,int> > mapping;
	};
	
	vector<CellSizingOption> timingCellSizingOptions;
    
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
		int node;     // index of the tree node which drives this arc
		
		TimingArc() {
			cell = NULL;
			lut = -1;
			pin = -1;
			driver = -1;
			sink = -1;
			node = -1;
		} // end constructor
	}; // end struct
	
	struct TimingNet {
		Vcell * driver;
		int fanout;                     // number of timing arcs driven by this net (used in the slew violation calculation)
		int depth;                      // logical depth of this net
		
		TimingNet() {
			driver = NULL;
			depth = -1;
			fanout = 0;
		} // end constructor
	}; // end struct
	
	struct TimingArcState {
		EdgeArray<double> islew; //  input slew at this arc
		EdgeArray<double> oslew; // output slew at this arc
		
		EdgeArray<double> delay; // delay of this arc
        EdgeArray<double> lambda; // Lagrange Relaxation
		
		EdgeArray<double> rcdelay; // delay at input of this arc due to rc tree
		EdgeArray<double> ceff; // effective capacitance driven by this arc
		
		EdgeArray<double> arrivalTime; // arrival time at the input of this arc
		EdgeArray<double> requiredTime; // required time at the input of this arc

		double loadViolation; // load violation based on ceff

		//Reimann
		EdgeArray<double> slack; // slack (m_u->v) (worst arrival time + arc delay + worst delay to endpoint) [TODO] consider moving to outside
		
		TimingArcState() {
            lambda.set(-1.0,-1.0);
		} // end constructor		
	};	
	
	struct TimingNetState {
		double load;                    // total net load (due to wires, input pins and primary outputs)
		double loadViolation;           // load violation based on ceff
		EdgeArray<double> slew;         // slew in this net
		EdgeArray<double> arrivalTime;  // arrival time at this net
		EdgeArray<double> requiredTime; // required time at the driver of this net
		EdgeArray<int> backtrack;       // index of net driving the driver input with greatest arrival time
		EdgeArray<int> backtrackSlew;   // index of net driving the driver input with greatest slew

		EdgeArray<double> worstRCDelay; // worst rc tree delay at sinks

		//Reimann
		EdgeArray<double> lambdaDelay;   // sum(lambda_arcs)
		
		TimingNetState() {
			load = 0;
			slew.set(0,0);
		} // end constructor		
	};
	
	struct TreeNodePointer {
		int arc;  // timing arc node driving by the tree node
		int node; // tree node index inside the tree
		
		TreeNodePointer(const int arcIndex, const int nodeIndex ) {
			arc = arcIndex;
			node = nodeIndex;
		} // end constructor
	}; // end class
	
	struct State {
		vector<TimingArcState> arcs;
		vector<TimingNetState> nets;
	};
	
	State timingStateCurrent;
	State timingStateStored;
	
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
    
	int timingOffsetToSequentialArcs;
	int timingOffsetToCombinationArcs;
	int timingOffsetToExtraSequentialArcs;
	int timingOffsetToExtraPrimaryOutputArcs;
	
	vector<double> timingMaxLoad;
	
	// Number of dummy nets. Use this as an offset to skip dummy nets as all
	// dummy nets have a index less than timingNumDummyNets. The clock net is
	// also considered as a dummy net. Clock net has index equal to
	// (timingNumDummyNets - 1).
	
	int timingNumDummyNets;
	int timingOffsetToLevelOneNets; // skip dummy and level=0 (primary inputs and ffs driven nets)
	
	vector<int> timingOffsetToNetLevel;
	
	// Stores the worst arrival time for both edges and the (tail) net with
	// this worst arrival time.
	EdgeArray<int>    timingWorstArrivalTimeArc;
	EdgeArray<double> timingWorstArrivalTime;
	
	// Timign violations.
	double timingViolationSlew;
	double timingViolationLoad;
	double timingViolationLoadCellWise;
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
	
	// RC Trees
	vector< RCTree > timingTrees;
	vector< RCTreeDescriptor > timingTreeDescriptors; 
	
	vector<int> timingTreeNodePointers;
	vector<TreeNodePointer> timingTreeNodes;
		
	void buildTimingStructure();
	void buildTreeStructure();
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

	// Update timing regardless RC Tree. Suppose RC Tree timing is linearly 
	// dependent of load.
	void updateTiming_Net_LinearApproximation( const int n );
	void updateTimingLocally_LinearApproximation(const int n);
	
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
	double computeSizingEffectOnLocalPositiveSlack(const int n);
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

	bool initialized;

    // LR iteration index
    double kIndex;
    double upperBound;
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
	
    multimap<double,int> pathMappedCells;
    multimap<int, double> slackMappedCells;
    set<int> slewViolCells;
    
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
	double timingViol, slewViol, minViol, loadViol, totalLeakage, bestCost, totalArea;
	
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
	
	//Restores from memory the stored sizing solution.
	void restoreFirstSolution();
	
	// Restores from memory the stored sizing solution.
	void restoreSolution();

	// Stores current net and arc states.
	void storeState();
	
	// Sorts cells by higher load capacitance
	void cloadOrder();
	
	// Save ceff file
	void saveCeff();
	
	// Save current solution to .sizes file
	void saveSizes();
	
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
	
	// Update totalArea
	void calcArea();
	
	// Update loadViol
	void calcLoadViol();
	
	// Update slewViol
	void calcSlewViol();
	
	// Report usage by Size and Vth
    void reportCellUsage();
	
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
	bool updateCellTypeLagrangeRelaxationLinearApproximation(Vcell * cell, const double alpha = 1.0);
	bool updateCellTypeLagrangeRelaxationSensitivitiesLinearApproximation(Vcell *cell, const double gamma, const double alpha = 1.0);
	bool updateCellTypeLagrangeRelaxationSensitivities(Vcell * cell, const double gamma, const double alpha = 1.0);
    bool updateCellTypeLagrangeRelaxationSensitivitiesDefault(Vcell * cell, const double gamma, const double alpha = 1.0);
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
	void printTiming(const string &label = "") const;
	
	void printDelaysVersusExpectedDelays();
	
	// Print current timing information.
	void printPathTails();
	
	// Print critical path timing.
	void printCriticalPathTimingReport(ostream &out);
	
	// Print some library stats.
	void printLibraryReport(ostream &out);
	
	// Generate a sight visualization file.
	void printSigth(const string &filename);
	
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

	// Compute load violation using effective capacitance.
	double computeLoadViolationUsingEffectiveCap( const TimingArc &arc, const TimingArcState &state ) const;
	double computeLoadViolationUsingEffectiveCapCellWise( const TimingArc &arc, const TimingArcState &state ) const;
	double computeLoadViolationUsingDownstreamCap( const Vcell * cell ) const;

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
	void sizingLagrangeRelaxation(const bool resetLambdas = true);
	void sizingLagrangeRelaxationLinearApproximation(const bool resetLambdas = true);
	void sizingLagrangeRelaxationSensitivitiesLinearApproximation(const bool resetLambdas = true);
    void sizingLagrangeRelaxationLinearApproximationTestingChanges(const bool resetLambdas = true);
	void sizingLagrangeRelaxationSensitivities(const bool resetLambdas = true);
    void sizingLagrangeRelaxationSensitivitiesDefault(const bool resetLambdas = true, const int iterations = 150);
    void sizingLagrangeRelaxationSensitivitiesNew(const bool resetLambdas = true);
	void sizingLagrangeRelaxationSensitivitiesSomeChanges(const bool resetLambdas = true);
    void sizingLagrangeRelaxationTestingChanges(const bool resetLambdas = true);
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
	void timingRecovery(const int limit);
	void timingRecoveryPathCounter();
	void timingRecoveryPrimeTime(const int iterationLimit);
	void timingRecoveryPrimeTime2();
	void timingRecoveryPathCounter(const int limit, const bool highEffort = false);
    void timingRecoveryPrimeTimeTestingChanges();
	void timingRecoveryPathCounterLimited(const int limit);
	void timingRecoveryPathCounterTestingChanges();
    void powerRecovery();
	void powerRecovery(const double limit);
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
	void updateLambdasByFlachAlpha(const double alpha);
	void updateLambdasByFlachReimann();
	void updateLambdasByFlachReimannAlpha(const double alpha);
	void updateLambdasByFlachWeighted();
	void updateLambdasCriticality();
	
    void updateLambdas_Subgradient();
	void updateLambdas_Normalization();
	void updateLambdas_KKT();
	void updateLambdas_Nets();
	
    void updateLambdasByOzdal();
    
	void resetTimingArcLambdas( const double value );
	
	void computeLambdaSTA();
	
    void callPT();
	void callPTNegSlackOnly();
	void callPTCeffOnly();
	void callPTNoReport();
	void callPT_TR();
	void callPT_PR();
	void callPTNonBlocking();
	void callPTNonBlockingNoReport();
	void readTimingLR();
	void readTimingFromPT();
    void readTimingFromPTForTimingRecovery();
    void legalizeLoadViolPrimeTime();
    void compareTimingEngines();
    
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
		const EdgeArray<double> ceff,
		EdgeArray<double> &outDelay,
		EdgeArray<double> &outSlew);

	static void computeArcTiming(
		const LibParserTimingInfo &timingInfo,
		const EdgeArray<double> inputSlew,
		const double load,
		EdgeArray<double> &outDelay,
		EdgeArray<double> &outSlew
	) {computeArcTiming(timingInfo, inputSlew, EdgeArray<double>(load,load), outDelay, outSlew);}
    
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
	static bool nearlyZero( const double v, const double precision = EPSILON  ) { return fabs(v) < precision; }
	
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

	const string &getBenchmarkName() const { return benchName; }
	int getNumTailNets() const { return timingTailNets.size(); }
	int getNumPathsWithNegativeSlack() const { return timingNumPathsWithNegativeSlack; }
	
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
		
	EdgeArray<double> getNetSlack(const int n) const { const TimingNetState &s = getTimingNetState(n); return s.requiredTime - s.arrivalTime; }
	EdgeArray<double> getNetNegativeSlack(const int n) const { const TimingNetState &s = getTimingNetState(n); return min(EdgeArray<double>(0,0), s.requiredTime - s.arrivalTime); }
	EdgeArray<double> getNetPositiveSlack(const int n) const { const TimingNetState &s = getTimingNetState(n); return max(EdgeArray<double>(0,0), s.requiredTime - s.arrivalTime); }
	
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
	
	string getTimingArcInputNodeName(const int k) {
		const TimingArc &arc = timingArcs[k];
		if ( k < 0 ) {
			assert(false); 
			return "bug";			
		} else if ( k < timingOffsetToSequentialArcs ) {
			// Timing Arc: Dummy driving primary inputs (path head).
			assert(false);
			return "bug";
		} else if ( k < timingOffsetToCombinationArcs ) {
			// Timing Arc: Dummy driving output of sequential elements (path head).
			assert(false);
			return "bug";
		} else if ( k < timingOffsetToExtraSequentialArcs ) {
			// Timing Arc: Combinational.
			return arc.cell->instName + ":" + arc.cell->actualInstType->pins[arc.pin].name;
		} else if ( k < timingOffsetToExtraPrimaryOutputArcs ) {
			// Timing Arc: Dummy driving sequential elements (path tail).
			return arc.cell->instName + ":" + arc.cell->actualInstType->pins[arc.pin].name;
		} else if ( k < timingArcs.size() ) {
			// Timing Arc: Dummy driven by primary outputs (path tail).
			return timingNetName[arc.driver];
		} else {
			assert(false);
			return "out-of-bounds";
		} // end else
	} // end method
	
	string getTimingArcOutputNodeName(const int k) {
		const TimingArc &arc = timingArcs[k];
		if ( k < 0 ) {
			assert(false); 
			return "bug";			
		} else if ( k < timingOffsetToSequentialArcs ) {
			// Timing Arc: Dummy driving primary inputs (path head).
			return timingNetName[arc.sink];
		} else if ( k < timingOffsetToCombinationArcs ) {
			// Timing Arc: Dummy driving output of sequential elements (path head).
			return timingNetName[arc.sink];
		} else if ( k < timingOffsetToExtraSequentialArcs ) {
			// Timing Arc: Combinational.
			return arc.cell->instName + ":o";
		} else if ( k < timingOffsetToExtraPrimaryOutputArcs ) {
			// Timing Arc: Dummy driving sequential elements (path tail).
			assert(false);
			return "bug";
		} else if ( k < timingArcs.size() ) {
			// Timing Arc: Dummy driven by primary outputs (path tail).
			assert(false);
			return "out-of-bounds";
		} else {
			assert(false);
			return "out-of-bounds";
		} // end else
	} // end method

	double getTimingViolation() const { return timingTotalNegativeSlack; }
	double getSlewViolation() const { return timingViolationSlew; }
	double getLoadViolationUsingEffectiveCap() const { return timingViolationLoad; }
	double getLoadViolationUsingEffectiveCapCellWise() const { return timingViolationLoadCellWise; }
	double getLoadViolationUsingDownstreamCap() const { return loadViol; }
	double getDeltaLoad(Vcell * cell, Vcell * nextCell, int newIndex);
	
	double getLeakagePower() const { return totalLeakage; }
	double getArea() const { return totalArea; }
	
	const vector< pair<string,string> > &getBestSolution() const { return bestSolution; }
	
	double getT() const { return sdcInfos.clk_period; };
	
    double runTimeLimit;
    
	// -------------------------------------------------------------------------
	
	void setT(const double T) { sdcInfos.clk_period = T; }
	void setDeltaLoadSlew(const int netIndex);
	
	// -------------------------------------------------------------------------
	// PrimeTime Remote Simulation
	// -------------------------------------------------------------------------	
	
#ifdef REMOTE_PRIMETIME	
	TCPSocket sock;
#endif
	
	bool primeTimeConnectionEstablished;
	
	void primeTimeReport();
	
	void primeTimeConnect(const string &serverAddress, const unsigned short serverPort);
	void primeTimeExec(const string &cmd);
	void primeTimeUpdateTiming();
	void primeTimeWait(const string &match = "", ostream * out = &cout);
	void primeTimeDisconnect();
	
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
		kIndex(0),
		initialized(false),
        totalArea(-1),
        runTimeLimit(-1)
    {}
	
	~Circuit() {
		primeTimeDisconnect();
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

