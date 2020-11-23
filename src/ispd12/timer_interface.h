//////////////////////////////////////////////////////////////////
//
//
//  Timing Analysis Interface helper class to interface with
//  the timer.
//
//  This code is provided for description purposes only. The contest
//  organizers cannot guarantee that the provided code is free of
//  bugs or defects. !!!! USE THIS CODE AT YOUR OWN RISK !!!!!
//
//
//  The contestants are free to use these functions as-is or make
//  modifications. If the contestants choose to use the provided
//  code, they are responsible for making sure that it works as
//  expected.
//
//  The code provided here has no real or implied warranties.
//
//
////////////////////////////////////////////////////////////////////

#ifndef _TIMERINTERFACE_H_
#define _TIMERINTERFACE_H_
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <dirent.h>
#include <cassert>
#include <cstdlib>

class TimerInterface {
	// This class contains functions for the timing analysis interface.
	// To use any function belonging to this class, call TimerInterface::<function_name>(<argument_list>);
	
	// Declarations
	// LOOK AT THIS PUBLIC SECTION - FUNCTION IMPLEMENTATIONS BELOW IN THIS FILE
public:
	// Status (outside of this class, you must use TimerInterface::Status to define variables of this type)
	enum Status { TIMER_NOT_STARTED = 0,   // Timing analysis has not been started
		TIMER_BUSY,              // Timer is busy (is reading design or performing timing analysis)
		TIMER_FINISHED_SUCCESS,  // Timing analysis finished successfully
		TIMER_FINISHED_ERROR,    // Error occured during timing analysis
		TIMER_INTERFACEERROR     // Error indicating that the program could not get timer status (could not read status file)
	};
	
	// Get timer status
	// Inputs: contest root directory (string)
	//         benchmark name (string)
	// Return: status (see enum Status above)
	static Status getTimerStatus(const std::string &contest_root, const std::string &benchmark);
	
	// Write sizes and run timing analysis in blocking mode
	// 1. Write sizes
	// 2. Starts timing analysis
	// 3. Waits for timing analysis to be completed
	// Inputs: vector of pairs where first value is instance name (string) and second value is cell name (string)
	//         contest root directory (string)
	//         benchmark name (string)
	//         polling time (number of seconds that the function should wait before polling timer status to check whether timer is done)
	// Return: timer status
	static Status runTimingAnalysisBlocking(const std::vector<std::pair<std::string, std::string> > &sizes, const std::string &contest_root, const std::string &benchmark, const unsigned pollingTime);
	
	// Start timing analysis in non-blocking mode
	// 1. Write sizes
	// 2. Starts timing analysis and returns (does not wait for timing analysis to be completed)
	// Inputs: vector of pairs where first value is instance name (string) and second value is cell name (string)
	//         contest root directory (string)
	//         benchmark name (string)
	// Return: timer status
	static Status startTimingAnalysisNonBlocking(const std::vector<std::pair<std::string, std::string> > &sizes, const std::string &contest_root, const std::string &benchmark);
	
	// Wait for given number of seconds (useful function if you want to wait before checking timer status after calling startTimingAnalysisNonBlocking)
	// Input: seconds to wait
	static void wait(int seconds);
	
	
	
	// PRIVATE SECTION --------------------------------------------------------------------------------
	// DO NOT LOOK AT THIS PRIVATE SECTION, YOU SHOULD ONLY LOOK AT FUNCTIONS DEFINED IN PUBLIC SECTION  
private:
	// Get timer status (helper function for isTimerDone)
	// Input:  vector of file names (returned by getFiles)
	// Return: string indicating timer status
	static std::string getTimerStatusString(const std::vector<std::string> &files);
	
	// Checks if a file exists (returns true if it does, false otherwise)
	// Input:  filename including path (string)
	// Return: true if the file exists and is readable, false otherwise
	static bool doesFileExist(const std::string &file);
	
	// Get a list of files from given directory (used by getTimerStatus to check if timer is done)
	// Input: directory name (string)
	// Output: vector of file names (strings), argument passed by reference
	// Return: true if directory could be read, false otherwise
	static bool getFiles(std::vector<std::string> &files, const std::string &dir);
	
	// Remove a file from the given directory (helper function used by startTimingAnalysis)
	// Inputs: name of the file without directory name (string)
	//         directory name (string)
	// Return: true if file was removed successfully, false otherwise
	static bool removeFile(const std::string &dir, const std::string &file);
	
	// Write sizes to a file for timing analysis call
	// Inputs: vector of pairs where first value is instance name (string) and second value is cell name (string)  
	//         contest root directory (string)
	//         benchmark name (string)
	// Return: true if sizes were written successfully to .int.sizes file, false otherwise
	static bool writeSizesForTimer(const std::vector<std::pair<std::string, std::string> > &sizes, const std::string &contest_root, const std::string &benchmark);
	
	// Start timing analysis (does not wait for it to finish)
	// Input:  contest root directory (string)
	//         benchmark name (string)
	// Return: true if successfully wrote command to start timing analysis, false otherwse
	static bool startTimingAnalysis(const std::string &contest_root, const std::string &benchmark);
	
	// Run timing analysis in blocking mode
	// 1. Starts timing analysis
	// 2. Waits for timing analysis to be completed
	// Input:  contest root directory (string)
	//         benchmark name (string)
	//         polling time (number of seconds that the function should wait before polling timer status to check whether timer is done)
	// Return: timer status
	static Status runTimingAnalysisBlocking(const std::string &contest_root, const std::string &benchmark, const unsigned pollingTime);
	
	// Start timing analysis in non-blocking mode
	// 1. Starts timing analysis and returns (does not wait for timing analysis to be completed)
	// Input:  contest root directory (string)
	//         benchmark name (string)
	// Return: timer status
	static Status startTimingAnalysisNonBlocking(const std::string &contest_root, const std::string &benchmark);
	
	// END PRIVATE SECTION ----------------------------------------------------------------------------
}; // END class TimerInterface



#endif // _TIMERINTERFACE_H_

