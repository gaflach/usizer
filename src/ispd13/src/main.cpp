/*
 *  main_cl.cpp
 *  gr
 *
 *  Created by Tiago Reimann on 27/06/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 
 
 to do:
 
 calcular delay das celulas? viavel? pq PT?
 algoritmo de troca d sizes
 
 
 */

#include <vector>
#include <set>
#include <queue>

#include "main.h"
#include <climits>
#include <cstring>
#include "Circuit.h"
#include "parser_helper.h"
#include "global.h"
#include "Stopwatch.h"

App::OptionMap App::clsOptions;
App app;

int main(int argc, char *argv[]) {
	Stopwatch runtime;
	runtime.start();
	
    if (argc < 3) {
        cerr << "[USAGE] " << argv[0] << " <benchmark dir> <benchmark name>\n";
        exit(1);
    } // end if
    
    app.parseCommandLineArguments(argc, argv, 4);
    //app.printCommandLineArguments(cerr);
    
    string benchmarkName, dirRoot;
    
    benchmarkName = argv[2];
    dirRoot = argv[1];
    
	const string teamName = "UFRGS-Brazil Gate Sizer";
	
	cout << setw(80) << setfill('*') << "\n" << setfill(' ');
	cout << setw((80 + teamName.length()) / 2) << teamName << "\n";
	cout << setw(80) << setfill('*') << "\n" << setfill(' ');
	
#ifdef PARALLEL
	cout << "Version: Multi-Threaded\n";
#else
	cout << "Version: Single-Threaded\n";
#endif
	
	cout << "Benchmark: " << benchmarkName << "\n";
	
	cout << setw(80) << setfill('*') << "\n" << setfill(' ');
	cout << "\n";
	
    int flowIndex = 3;
    const bool fast = ( (argc >= 4) && !strcmp(argv[3],"-fast") );
    const bool load = ( (argc >= 4) && !strcmp(argv[3],"-load") );
    if (load) {
        flowIndex++;
        cout << "\nLOAD MODE ENABLED!\n" << endl;
        flowDefaultLoad(benchmarkName, dirRoot);
    }
    else if (fast) {
        flowIndex++;
        cout << "\nFAST MODE ENABLED!\n" << endl;
        flowDefaultFast(benchmarkName, dirRoot);
    }
    else {
        const string flow = (argc >= flowIndex+1) ? argv[flowIndex] : "default";
        if (flow == "tiago")
            flowTiago(benchmarkName, dirRoot);
        else if (flow == "graci")
            flowGraci(benchmarkName, dirRoot);
        else if (flow == "flach")
            flowFlach(benchmarkName, dirRoot);
        else if (flow == "johann")
            flowJohann(benchmarkName, dirRoot);
        else if (flow == "default")
            flowDefault(benchmarkName, dirRoot);
        else {
            cout << "[WARNING] Invalid flow name '" << flow << "'. Using default flow.\n";
            flowDefault(benchmarkName, dirRoot);
        } // end else
    }
	runtime.stop();
	cout << "Total runtime (m): " << (runtime.getElapsedTime()/60) << "\n";
    
} // end method


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void flowTiago(string benchmarkName, string dirRoot) {
    
    Circuit myCircuit;
    
    cout << benchmarkName << endl;
    
    myCircuit.benchName = benchmarkName;
    myCircuit.rootDir = dirRoot;
	
    
    myCircuit.readInputFiles();
    
	const double T = myCircuit.getT();
    const int numGates = myCircuit.icells.size();
    
    //myCircuit.printLibraryReport(cerr);
	
#ifdef PARALLEL
	/*
	 myCircuit.initialSolution(false); //set initial solution (min leakage cells ou random cells)
	 
	 const int N = 1;
	 
	 Stopwatch timer0;
	 Stopwatch timer1;
	 
	 timer1.start();
	 for ( int i = 0; i < N; i++ ) {
	 myCircuit.updateTimingMultiThreaded();
	 } // end for
	 timer1.stop();
	 myCircuit.printTiming("Multi Threaded");
	 cerr << "Multi Threaded Runtime (s): " << timer1.getElapsedTime() << "\n";
	 
	 timer0.start();
	 for ( int i = 0; i < N; i++ ) {
	 myCircuit.updateTimingSingleThreaded();
	 } // end for
	 timer0.stop();
	 myCircuit.printTiming("Single Threaded");
	 cerr << "Single Threaded Runtime (s): " << timer0.getElapsedTime() << "\n";
	 
	 return;
	 */
#endif
	
    cout << "*** Initial Solution ***" << endl;
    myCircuit.initialSolution(false); //set initial solution (min leakage cells ou random cells)
    //myCircuit.calcTiming();
    myCircuit.updateTiming();
    myCircuit.printTiming("Initial Solution");
	
	myCircuit.sizingForNoLoadAndSlewViolationByLiLi();
	//myCircuit.sizingForNoLoadAndSlewViolationByLivramento();
	//myCircuit.sizingForNoLoadViolation();
	//myCircuit.sizingByLiLi();
	
	myCircuit.updateTiming();
	myCircuit.printTiming("Sizing for No Load/Slew Violations");
	
    myCircuit.timingRecoveryPrimeTime(15);
    return;
    
	myCircuit.sizingLagrangeRelaxationSensitivitiesNew();
    
	myCircuit.saveSolution("after-lagrange");
	
	myCircuit.timingRecoveryPathCounter();
    
	//myCircuit.timingRecovery();
	//myCircuit.sizingCriticalPathSensitivity();
	myCircuit.saveSolution("after-lagrange-tr");
	
	//myCircuit.timingRecovery();
	//myCircuit.sizingGreedyPower();
	
	//myCircuit.lowTemperatureAnneal2();
	//myCircuit.saveSolution("after-sa");
	
	myCircuit.powerRecovery();
	myCircuit.saveSolution("after-sa-pr");
	
	myCircuit.saveSizes();
	
    //myCircuit.calcTiming();
    myCircuit.updateTiming();
    myCircuit.printTiming("Final Solution");
	
	return;
    
    cout << "Benchmark: " << benchmarkName << endl;
    
    myCircuit.benchName = benchmarkName;
    myCircuit.rootDir = dirRoot;
    
    myCircuit.readInputFiles();
	
    srand(2012); //time(NULL));
    myCircuit.initialSolution(false); //set initial solution (min leakage cells ou random cells)
    
    const long double minLeakage = myCircuit.totalLeakage;
    cout << "Initial Leakage: " << minLeakage << endl;
    cout << myCircuit.getSize() << " " << (myCircuit.getSize() / 35000.0) << " " << ceil(myCircuit.getSize() / 35000.0) << endl;
    
    const double timeLimit = 5 * 60 * 60 + 1 * 60 * 60 * ceil(myCircuit.getSize() / 35000.0);
    cout << " Time limit: " << timeLimit << endl;
    
    myCircuit.LR();
    /*
	 cout << "\nSizingNoLoadViolation..." << endl;
	 myCircuit.sizingForNoLoadViolation();
	 cout << "SizingNoLoadViolation...done\n" << endl;
	 
	 cout << "\nCalculating timing..." << endl;
	 //myCircuit.calcTiming();
	 myCircuit.updateTimingLR_KKT();
	 myCircuit.printTiming();
	 cout << "Calculating timing...done\n" << endl;
	 
	 cout << "\nSizingDepthFanoutXLogicalEffort..." << endl;
	 myCircuit.sizingDepthFanoutXLogicalEffort();
	 //    myCircuit.updateTimingLR();
	 //myCircuit.printTiming();
	 cout << "SizingDepthFanoutXLogicalEffort...done\n" << endl;
	 
	 myCircuit.saveSizes();
	 
	 myCircuit.anneal();
	 */
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void flowGraci(string benchmarkName, string dirRoot) {
    
    cout << "Graciiiiiii" << endl;
    
    Circuit myCircuit;
    
    cout << benchmarkName << endl;
    
    myCircuit.benchName = benchmarkName;
    myCircuit.rootDir = dirRoot;
	
    
    myCircuit.readInputFiles();
    
	const double T = myCircuit.getT();
    const int numGates = myCircuit.icells.size();
	
    cout << "*** Initial Solution ***" << endl;
    myCircuit.initialSolution(false); //set initial solution (min leakage cells ou random cells)
    //myCircuit.calcTiming();
    myCircuit.updateTiming();
    myCircuit.printTiming("Initial Solution");
	
	myCircuit.sizingForNoLoadAndSlewViolationByLiLi();
	//myCircuit.sizingForNoLoadAndSlewViolationByLivramento();
	//myCircuit.sizingForNoLoadViolation();
	//myCircuit.sizingByLiLi();
    
    myCircuit.updateTiming();
	myCircuit.printTiming("Sizing for No Load/Slew Violations");
    //myCircuit.callPTNonBlocking();
    myCircuit.callPTNonBlockingNoReport();
    
	myCircuit.sizingLagrangeRelaxationLinearApproximationTestingChanges();
    //myCircuit.sizingLagrangeRelaxationSensitivitiesLinearApproximation();
	//myCircuit.sizingLagrangeRelaxationSensitivitiesDefault();
    //myCircuit.sizingLagrangeRelaxationSensitivities();
    myCircuit.callPTNoReport();
    
	//myCircuit.saveSolution("after-lagrange");
	
	//myCircuit.primeTimeUpdateTiming();
	//myCircuit.primeTimeExec("run_timing");
	//myCircuit.primeTimeExec("report_timing");
	//myCircuit.primeTimeWait("slack (");
	
    /*
     myCircuit.timingRecovery();
     
     myCircuit.primeTimeUpdateTiming();
     myCircuit.primeTimeExec("run_timing");
     myCircuit.primeTimeExec("report_timing");
     myCircuit.primeTimeWait("slack (");
     */
    
	myCircuit.timingRecoveryPathCounterTestingChanges();
	//myCircuit.timingRecoveryPathCounterLimited(myCircuit.getSize()/100);
	myCircuit.callPTNoReport();
    
    //	myCircuit.primeTimeUpdateTiming();
    //	myCircuit.primeTimeExec("run_timing");
    //	myCircuit.primeTimeExec("report_timing");
    //	myCircuit.primeTimeWait("slack (");
	
	//myCircuit.timingRecovery();
	//myCircuit.sizingGreedyPower();
	
	//myCircuit.lowTemperatureAnneal2();
	//myCircuit.saveSolution("after-sa");
	
	myCircuit.powerRecovery();
	//myCircuit.saveSolution("after-sa-pr");
	myCircuit.callPTNoReport();
    
    //	myCircuit.primeTimeUpdateTiming();
    //	myCircuit.primeTimeExec("run_timing");
    //	myCircuit.primeTimeExec("report_timing");
    //	myCircuit.primeTimeWait("slack (");
    
    myCircuit.timingRecoveryPrimeTimeTestingChanges();
	
	myCircuit.saveSizes();
    
    //myCircuit.calcTiming();
    myCircuit.updateTiming();
    myCircuit.printTiming("Final Solution");
    //myCircuit.callPTNegSlackOnly();
    myCircuit.callPTNoReport();
    
	return;
    
	
	myCircuit.sizingByLiLi();
    myCircuit.updateTiming();
    myCircuit.printTiming("Li Li");
	
	myCircuit.sizingLagrangeRelaxation();
    myCircuit.updateTiming();
    myCircuit.printTiming("Lambda Greedy");
	
	return;
	
	
	myCircuit.setT(T);
	
	//myCircuit.sizingByLiLi();
    myCircuit.sizingForNoLoadAndSlewViolationByLiLi();
    //myCircuit.sizingDepthFanoutXLogicalEffort();
	
	myCircuit.sizingLagrangeRelaxation();
	
    
    myCircuit.updateTiming();
    myCircuit.printTiming("LR");
	
	myCircuit.setT(T*0.95);
	
	myCircuit.sizingCriticalPathSensitivity();
	
    
    myCircuit.updateTiming();
    myCircuit.printTiming("TR");
	
	myCircuit.powerRecovery();
	
    myCircuit.updateTiming();
    myCircuit.printTiming("PR");
	
	
	return;
	
	//myCircuit.sizingByLiLi();
    
    //    myCircuit.sizingLagrangeRelaxationTestingChanges();
    //    myCircuit.updateTiming();
    //    myCircuit.printTiming("Lambda Greedy testing changes");
	
	//    myCircuit.printTimingDigest(true, "Power Recovery");
	//	for ( int i = 0; i < 100; i++ )
	//		if ( !myCircuit.powerRecovery() )
	//			break;
    
	
	
	//myCircuit.updateTiming();
	//myCircuit.printTiming("Sizing critical Path Sensitivity");
	
	//myCircuit.lowTemperatureAnnealTestedWithUFSC();
	
	//	myCircuit.printTimingDigest(true, "Reducingleakage LiLi");
	//	for ( int i = 0; i < 100; i++ )
	//		if ( !myCircuit.reducingLeakageByLiLi() )
	//			break;
	
	//myCircuit.lowTemperatureAnnealTestedWithUFSC();
	
	
	//myCircuit.printTimingDigest(true, "Power Recovery");
	//myCircuit.powerRecovery();
	//     for ( int i = 0; i < 100; i++ )
	//		 if ( !myCircuit.powerRecovery() )
	//			 break;
	
	//    myCircuit.timingRecovery();
	
	
    //myCircuit.lowTemperatureAnneal();
	
	//myCircuit.powerRecovery();
	
	//myCircuit.sizingCriticalPathSensitivity();
	
	//myCircuit.lowTemperatureAnneal();
	
	myCircuit.saveSizes();
	
    //myCircuit.calcTiming();
    myCircuit.updateTiming();
    myCircuit.printTiming("Final Solution");
	
	return;
    
}

////////////////////////////////////////////////////////////////////////////////

void flowFlach(string benchmarkName, string dirRoot) {
    Circuit myCircuit;
    
    cout << benchmarkName << endl;
    
    myCircuit.benchName = benchmarkName;
    myCircuit.rootDir = dirRoot;
	
    
    myCircuit.readInputFiles();
    
	const double T = myCircuit.getT();
    const int numGates = myCircuit.icells.size();
    
    //myCircuit.printLibraryReport(cerr);
	
#ifdef PARALLEL
	/*
	 myCircuit.initialSolution(false); //set initial solution (min leakage cells ou random cells)
	 
	 const int N = 1;
	 
	 Stopwatch timer0;
	 Stopwatch timer1;
	 
	 timer1.start();
	 for ( int i = 0; i < N; i++ ) {
	 myCircuit.updateTimingMultiThreaded();
	 } // end for
	 timer1.stop();
	 myCircuit.printTiming("Multi Threaded");
	 cerr << "Multi Threaded Runtime (s): " << timer1.getElapsedTime() << "\n";
	 
	 timer0.start();
	 for ( int i = 0; i < N; i++ ) {
	 myCircuit.updateTimingSingleThreaded();
	 } // end for
	 timer0.stop();
	 myCircuit.printTiming("Single Threaded");
	 cerr << "Single Threaded Runtime (s): " << timer0.getElapsedTime() << "\n";
	 
	 return;
	 */
#endif
    
    cout << "*** Initial Solution ***" << endl;
    myCircuit.initialSolution(false); //set initial solution (min leakage cells ou random cells)
    //myCircuit.calcTiming();
    myCircuit.updateTiming();
    myCircuit.printTiming("Initial Solution");
	
    myCircuit.sizingForNoLoadAndSlewViolationByLiLi();
	//myCircuit.sizingForNoLoadAndSlewViolationByLivramento();
	//myCircuit.sizingForNoLoadViolation();
	//myCircuit.sizingByLiLi();
    
	myCircuit.updateTiming();
	myCircuit.printTiming("Sizing for No Load/Slew Violations");
    
	// #########################################################################
	
    //	myCircuit.primeTimeUpdateTiming();
    //	myCircuit.primeTimeExec("report_timing");
    //	myCircuit.primeTimeWait("slack (");
    //
    //	myCircuit.sizingRandom();
    //	myCircuit.primeTimeUpdateTiming();
    //	myCircuit.primeTimeExec("run_timing");
    //	myCircuit.primeTimeExec("report_timing");
    //	myCircuit.primeTimeWait("slack (");
    //
    //	exit(8);
	
	// #########################################################################
    
    //	myCircuit.sizingByLiLi();
    //	myCircuit.primeTimeUpdateTiming();
    //	myCircuit.primeTimeExec("run_timing");
    //	myCircuit.primeTimeExec("report_timing");
    //	myCircuit.primeTimeWait("slack (");
    
	//myCircuit.sizingLagrangeRelaxationSensitivities();
	myCircuit.sizingLagrangeRelaxationLinearApproximation();
	//myCircuit.sizingLagrangeRelaxationSensitivitiesLinearApproximation();
	//myCircuit.sizingLagrangeRelaxationSensitivities();
    
	//myCircuit.saveSolution("after-lagrange");
    
	myCircuit.timingRecovery();
	
	//myCircuit.timingRecovery();
	//myCircuit.sizingCriticalPathSensitivity();
	//myCircuit.saveSolution("after-lagrange-tr");
	
	//myCircuit.timingRecovery();
	//myCircuit.sizingGreedyPower();
	
	//myCircuit.lowTemperatureAnneal2();
	//myCircuit.saveSolution("after-sa");
	
	myCircuit.powerRecovery();
	//myCircuit.saveSolution("after-sa-pr");
    
	myCircuit.saveSizes();
	
    //myCircuit.calcTiming();
    myCircuit.updateTiming();
    myCircuit.printTiming("Final Solution");
    
	return;
	
	myCircuit.sizingByLiLi();
    myCircuit.updateTiming();
    myCircuit.printTiming("Li Li");
	
	myCircuit.sizingLagrangeRelaxation();
    myCircuit.updateTiming();
    myCircuit.printTiming("Lambda Greedy");
	
	
	//myCircuit.printDelaysVersusExpectedDelays();
	
	//myCircuit.anneal();
	
	//myCircuit.printDelaysVersusExpectedDelays();
	
    return;
	
	/*
	 myCircuit.sizingByLiLi();
	 myCircuit.updateTiming();
	 myCircuit.printTiming("LiLi");
	 
	 myCircuit.saveSizes();
	 
	 return;
	 
	 for ( int i = 0; i < 1; i++) {
	 myCircuit.updateRequiredTime();
	 myCircuit.updateLambdasByTennakoon();
	 
	 myCircuit.sizingWorstLambda();
	 myCircuit.updateTiming();
	 myCircuit.printTiming("Solution");
	 
	 //myCircuit.enableLambdaAwareSTA = true;
	 //myCircuit.updateTiming();
	 //myCircuit.printTiming("Solution");
	 //myCircuit.enableLambdaAwareSTA = false;
	 //myCircuit.updateTiming();
	 
	 } // end for
	 
	 return;
	 */
	
    /*
     myCircuit.sizingRandom();
     myCircuit.calcTiming();
     myCircuit.updateTiming();
     myCircuit.printTiming("Random Solution");
     //*/
    
    /*
     myCircuit.sizingForNoLoadViolation();
     myCircuit.calcTiming();
     myCircuit.updateTiming();
     myCircuit.printTiming("No Load Solution");
     //*/
    
    /*
     myCircuit.assignLogicalEffort();
     myCircuit.sizingDepthFanoutXLogicalEffort();
     myCircuit.calcTiming();
     myCircuit.updateTiming();
     myCircuit.printTiming("Logical Effort");
     //*/
    
    /*
     myCircuit.sizingDepthFanoutXLogicalEffort();
     myCircuit.calcTiming();
     myCircuit.updateTiming();
     myCircuit.printTiming("sizingDepthFanoutXLogicalEffort");
     //*/
    
    /*
     myCircuit.sizingDepthFanout();
     myCircuit.calcTiming();
     myCircuit.updateTiming();
     myCircuit.printTiming("Logical Effort");
     //*/
    
    /*
     myCircuit.sizingWorstPath();
     myCircuit.printTiming("Worst Path");
     //*/
    
    //myCircuit.printLibraryReport(cerr);
    
    //myCircuit.sizingForLeakageReduction();
    
    //myCircuit.sizingForward();
    //myCircuit.printTiming("Sizing Forward");
    
    // *************************************************************************
    // Timing Engines Comparison
    // *************************************************************************
    /*
     const int N = 1000;
     vector< pair<Vcell*,int> > values( N );
     for ( int i = 0; i < N; i++ ) {
     values[i] = make_pair(
     myCircuit.chooseCellRandomlyFromVector(myCircuit.icells),
     (int)floor((rand()/double(RAND_MAX))*(30-1) + 0.5)
     );
     } // end for
     
     // ---
     
     for ( int i = 0; i < myCircuit.icells.size(); i++ )
     myCircuit.updateCellType( myCircuit.icells[i], 0 );
     myCircuit.calcTiming();
     myCircuit.updateTiming();
     
     Stopwatch timer1;
     cerr << "updateTiming() started\n";
     timer1.start();
     for ( int i = 0; i < N; i++ ) {
     myCircuit.updateCellType( values[i].first, values[i].second );
     myCircuit.updateTiming(values[i].first);
     }
     timer1.stop();
     myCircuit.printTiming();
     cerr << "updateTiming() ended: " << timer1.getElapsedTime() << "\n";
     
     // ---
     
     for ( int i = 0; i < myCircuit.icells.size(); i++ )
     myCircuit.updateCellType( myCircuit.icells[i], 0 );
     myCircuit.calcTiming();
     myCircuit.updateTiming();
     
     Stopwatch timer0;
     cerr << "calcTiming() started\n";
     timer0.start();
     for ( int i = 0; i < N; i++ ) {
     myCircuit.updateCellType( values[i].first, values[i].second );
     myCircuit.coneTiming2(values[i].first);
     } // end for
     timer0.stop();
     myCircuit.printTiming();
     cerr << "calcTiming() ended: " << timer0.getElapsedTime() << "\n";
     
     // ---
     
     cerr << "Ratio: " << (timer0.getElapsedTime()/ double(timer1.getElapsedTime())) << "\n";
     
     myCircuit.calcTiming();
     myCircuit.updateTiming();
     myCircuit.printTiming();
     
     exit(8);
     //*/
    
    /*
     for ( int i = 0; i < 10000; i++ ) {
     cout << "*** Cone Timing Test " << i << " ***" << endl;
     
     Vcell * cell = myCircuit.chooseCellRandomlyFromVector(myCircuit.icells);
     const int randomTypeIndex =  (rand()/double(RAND_MAX))*(30-1) + 0.5;
     
     
     //Vcell * cell = myCircuit.icells[121];
     //const int randomTypeIndex =  10;
     
     
     myCircuit.updateCellType(cell, randomTypeIndex);
     
     cerr << "\t*** " << cell->vectorIndex << "\t" << randomTypeIndex << "\n";
     
     myCircuit.calcTiming();
     myCircuit.updateTiming(cell);
     
     //myCircuit.printTiming();
     }// end for
     
     exit(8);
     //*/
    
    /*
     myCircuit.printSigth("out-initial.sight");
     
     cout << "*** No Load Violation Solution ***" << endl;
     myCircuit.sizingForNoLoadViolation();
     myCircuit.calcTiming();
     
     myCircuit.updateTiming();
     myCircuit.printTiming();
     
     myCircuit.printSigth("out-no-load-violation.sight");
     
     exit(8);
     */
    
    //myCircuit.calcTiming();
    myCircuit.updateTiming();
    myCircuit.printTiming("Final Solution");
    myCircuit.printSigth("out.sight");
    myCircuit.saveSizes();
    
}

////////////////////////////////////////////////////////////////////////////////

void flowJohann(string benchmarkName, string dirRoot) {
    // Never happens :)
}

////////////////////////////////////////////////////////////////////////////////

void flowDefault(string benchmarkName, string dirRoot) {
    
    Stopwatch runtime;
	runtime.start();
	
    Circuit myCircuit;
    
    
    myCircuit.benchName = benchmarkName;
    myCircuit.rootDir = dirRoot;
	
    
    myCircuit.readInputFiles();
    
	const double T = myCircuit.getT();
    const int numGates = myCircuit.getSize();
    
    const double runTimeLimit = 3.0*60.0*60.0 + 60.0*60.0*ceil(numGates/40000.0);
    myCircuit.runTimeLimit = runTimeLimit;
    cout << benchmarkName << endl;
    cout << runTimeLimit << endl;
    
    cout << "*** Initial Solution ***" << endl;
    myCircuit.initialSolution(false); //set initial solution (min leakage cells ou random cells)
    myCircuit.updateTiming();
    myCircuit.printTiming("Initial Solution");
	
    myCircuit.sizingForNoLoadAndSlewViolationByLiLi();
	
	myCircuit.updateTiming();
	myCircuit.printTiming("Sizing for No Load/Slew Violations");
    myCircuit.callPTNonBlockingNoReport();
    //myCircuit.compareTimingEngines();
    
    //myCircuit.sizingLagrangeRelaxationLinearApproximation();
	myCircuit.sizingLagrangeRelaxationSensitivitiesDefault();
    //myCircuit.callPTNoReport();
	/*
     myCircuit.compareTimingEngines();
     //myCircuit.timingRecoveryPathCounter(900+numGates/1000.0);
     myCircuit.timingRecoveryPathCounter((runTimeLimit - runtime.getElapsedTime())/20.0);
     //myCircuit.timingRecovery((runTimeLimit - runtime.getElapsedTime())/3.0);
     myCircuit.timingRecoveryPathCounter((runTimeLimit - runtime.getElapsedTime())/3.0,true);
     myCircuit.callPTNoReport();
     myCircuit.timingRecoveryPrimeTime(5);
     myCircuit.timingRecoveryPrimeTime2();
     //myCircuit.callPTNonBlockingNoReport();
     myCircuit.updateTiming();
     //myCircuit.powerRecovery(900+numGates/1000.0);
     */
	//myCircuit.callPT_TR();
    //myCircuit.loadSizes();
    myCircuit.powerRecovery((runTimeLimit - runtime.getElapsedTime())/2.0);
    /*	myCircuit.callPTNoReport();
     */
    
	myCircuit.callPT_PR();
    myCircuit.loadSizes();
    myCircuit.legalizeLoadViolPrimeTime();
    myCircuit.timingRecoveryPrimeTime(2);
    myCircuit.callPT_TR();
    myCircuit.loadSizes();
    myCircuit.legalizeLoadViolPrimeTime();
    myCircuit.timingRecoveryPrimeTime(5);
    myCircuit.timingRecoveryPrimeTime2();
    
	myCircuit.saveSizes();
	//myCircuit.saveCeff();
    
    myCircuit.updateTiming();
    myCircuit.printTiming("Final Solution");
    //myCircuit.callPTNoReport();
    //myCircuit.compareTimingEngines();
    
	return;
    
} // end method
////////////////////////////////////////////////////////////////////////////////

void flowDefaultFast(string benchmarkName, string dirRoot) {
    
    Stopwatch runtime;
	runtime.start();
	
    Circuit myCircuit;
    
    myCircuit.benchName = benchmarkName;
    myCircuit.rootDir = dirRoot;
	
    myCircuit.readInputFiles();
    
	const double T = myCircuit.getT();
    const int numGates = myCircuit.getSize();
    
    const double runTimeLimit = 60.0*ceil(((3.0 + ceil(numGates/40000.0))/5.0)*60.0);
    myCircuit.runTimeLimit = runTimeLimit;
    cout << benchmarkName << endl;
    cout << runTimeLimit << endl;
    
    cout << "*** Initial Solution ***" << endl;
    myCircuit.initialSolution(false); //set initial solution (min leakage cells ou random cells)
    myCircuit.updateTiming();
    myCircuit.printTiming("Initial Solution");
	
    myCircuit.sizingForNoLoadAndSlewViolationByLiLi();
	
	myCircuit.updateTiming();
	myCircuit.printTiming("Sizing for No Load/Slew Violations");
    myCircuit.callPTNonBlockingNoReport();
    //myCircuit.compareTimingEngines();
    
    //myCircuit.sizingLagrangeRelaxationLinearApproximation();
	myCircuit.sizingLagrangeRelaxationSensitivitiesDefault();
    //myCircuit.callPTNoReport();
	/*
     myCircuit.compareTimingEngines();
     //myCircuit.timingRecoveryPathCounter(900+numGates/1000.0);
     myCircuit.timingRecoveryPathCounter((runTimeLimit - runtime.getElapsedTime())/20.0);
     //myCircuit.timingRecovery((runTimeLimit - runtime.getElapsedTime())/3.0);
     myCircuit.timingRecoveryPathCounter((runTimeLimit - runtime.getElapsedTime())/3.0,true);
     myCircuit.callPTNoReport();
     myCircuit.timingRecoveryPrimeTime(5);
     myCircuit.timingRecoveryPrimeTime2();
     //myCircuit.callPTNonBlockingNoReport();
     myCircuit.updateTiming();
     //myCircuit.powerRecovery(900+numGates/1000.0);
     */
	//myCircuit.callPT_TR();
    //myCircuit.loadSizes();
    myCircuit.powerRecovery((runTimeLimit - runtime.getElapsedTime())/2.0);
    /*	myCircuit.callPTNoReport();
     */
    
	myCircuit.callPT_PR();
    myCircuit.loadSizes();
    myCircuit.legalizeLoadViolPrimeTime();
    myCircuit.timingRecoveryPrimeTime(2);
    myCircuit.callPT_TR();
    myCircuit.loadSizes();
    myCircuit.legalizeLoadViolPrimeTime();
    myCircuit.timingRecoveryPrimeTime(5);
    myCircuit.timingRecoveryPrimeTime2();
    
	myCircuit.saveSizes();
	//myCircuit.saveCeff();
    
    myCircuit.updateTiming();
    myCircuit.printTiming("Final Solution");
    //myCircuit.callPTNoReport();
    //myCircuit.compareTimingEngines();
    
	return;
    
} // end method
////////////////////////////////////////////////////////////////////////////////

void flowDefaultLoad(string benchmarkName, string dirRoot) {
    Circuit myCircuit;
    
    
    myCircuit.benchName = benchmarkName;
    myCircuit.rootDir = dirRoot;
	
    
    myCircuit.readInputFiles();
    
	const double T = myCircuit.getT();
    const int numGates = myCircuit.getSize();
    
    const double runTimeLimit = 3.0*60.0*60.0 + 60.0*60.0*ceil(numGates/40000.0);
    cout << benchmarkName << endl;
    cout << runTimeLimit << endl;
    
    cout << "*** Initial Solution ***" << endl;
    myCircuit.initialSolution(false); //set initial solution (min leakage cells ou random cells)
    myCircuit.updateTiming();
    myCircuit.printTiming("Initial Solution");
	/*
     myCircuit.sizingForNoLoadAndSlewViolationByLiLi();
     
     myCircuit.updateTiming();
     myCircuit.printTiming("Sizing for No Load/Slew Violations");
     myCircuit.callPTNonBlockingNoReport();
     
     //myCircuit.sizingLagrangeRelaxationLinearApproximation();
     myCircuit.sizingLagrangeRelaxationSensitivitiesDefault(true, 50);
     //myCircuit.callPTNoReport();
     
     myCircuit.timingRecoveryPathCounter(900+numGates/1000.0);
     myCircuit.callPTNonBlockingNoReport();
     
     myCircuit.powerRecovery(900+numGates/1000.0);
     //myCircuit.callPTNoReport();
     */
    myCircuit.loadSizes();
    
    myCircuit.legalizeLoadViolPrimeTime();
    myCircuit.timingRecoveryPrimeTime(15);
    
	myCircuit.saveSizes();
	//myCircuit.saveCeff();
    
    myCircuit.updateTiming();
    myCircuit.printTiming("Final Solution");
    //myCircuit.callPTNoReport();
    return;
    
} // end method
