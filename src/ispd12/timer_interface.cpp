#include "timer_interface.h"

// Function Definitions -----------------------------------------------------------------------------
// Get timer status
TimerInterface::Status TimerInterface::getTimerStatus(const std::string &contest_root, const std::string &benchmark) {
	const std::string dir = contest_root + "/" + benchmark;
	
	// Get files from directory to see if there are any timer status
	std::vector<std::string> files;
	if (!getFiles(files,dir)) {
		return TIMER_INTERFACEERROR;
	}
	
	// Get command and perform action
	std::string cmd = getTimerStatusString(files);
	if ("__SIZERCMD_TIMERERROR_" == cmd) {
		return TIMER_FINISHED_ERROR;
	} else if ("__SIZERCMD_TIMERDONE_" == cmd) {
		return TIMER_FINISHED_SUCCESS;
	} else if ("__TCMD_RUNTIMER_" == cmd) {
		return TIMER_BUSY;
	} else {
		return TIMER_NOT_STARTED;
	}
}

// Write sizes and run timing analysis in blocking mode
TimerInterface::Status TimerInterface::runTimingAnalysisBlocking(const std::vector<std::pair<std::string, std::string> > &sizes, const std::string &contest_root, const std::string &benchmark, const unsigned pollingTime) {
	if (!writeSizesForTimer(sizes, contest_root, benchmark)) {
		std::cout << "-E- runTimingAnalysisBlocking: problem writing sizes" << std::endl;
		return TIMER_INTERFACEERROR;
	}
	return runTimingAnalysisBlocking(contest_root,benchmark,pollingTime);
}

// Start timing analysis in non-blocking mode
TimerInterface::Status TimerInterface::startTimingAnalysisNonBlocking(const std::vector<std::pair<std::string, std::string> > &sizes, const std::string &contest_root, const std::string &benchmark) {
	if (!writeSizesForTimer(sizes, contest_root, benchmark)) {
		std::cout << "-E- startTimingAnalysisNonBlocking: problem writing sizes" << std::endl;
		return TIMER_INTERFACEERROR;
	}
	return startTimingAnalysisNonBlocking(contest_root,benchmark);  
}

// Wait for given number of seconds (useful function if you want to wait before checking timer status after calling startTimingAnalysisNonBlocking)
void TimerInterface::wait(int seconds) {
	std::ostringstream ostr;
	ostr << seconds;
	system(("sleep " + ostr.str()).c_str());
}

// PRIVATE SECTION --------------------------------------------------------------------------------
// Get timer status (helper function for isTimerDone)
std::string TimerInterface::getTimerStatusString(const std::vector<std::string> &files) {
	std::string cmd = "";
	for (unsigned i=0; i<files.size(); ++i) {
		if ("__SIZERCMD_TIMERERROR_" == files[i] ||
			"__SIZERCMD_TIMERDONE_" == files[i] ||
			"__TCMD_RUNTIMER_" == files[i]) {
			if (cmd != "") {
				std::cout << "-Error- getTimerStatusString: multiple status found" << std::endl;
				for (unsigned j=0; j<files.size(); ++j) {
					if ("__SIZERCMD_TIMERERROR_" == files[j] ||
						"__SIZERCMD_TIMERDONE_" == files[j] ||
						"__TCMD_RUNTIMER_" == files[j]) {
						std::cout << "   Status File: " << files[j] << std::endl;
					}
				}
				assert(false);
			}
			cmd = files[i];
		}
	}
	return cmd;
}

// Checks if a file exists (returns true if it does, false otherwise)
bool TimerInterface::doesFileExist(const std::string &file) {
	std::ifstream infile(file.c_str());
	if (!infile) {
		return false;
	}
	infile.close();
	return true;
}

// Get a list of files from given directory (used by getTimerStatus to check if timer is done)
bool TimerInterface::getFiles(std::vector<std::string> &files, const std::string &dir) {
	files.clear();
	DIR *d = opendir(dir.c_str());
	if (NULL == d) {
		std::cout << "-E- getFiles: could not list files in directory '" << dir << "' to get timer status" << std::endl;
		return false;
	}
	files.clear();
	dirent *f = readdir(d);
	while (NULL != f) {    
		files.push_back(f->d_name);
		f = readdir(d);
	}
	closedir(d);
	//for (unsigned i=0; i<files.size(); ++i) {
	//  std::cout << files[i] << std::endl;
	//}
	return true;
}

// Remove a file from the given directory (helper function used by startTimingAnalysis)
bool TimerInterface::removeFile(const std::string &dir, const std::string &file) {
	if (file != "__SIZERCMD_TIMERERROR_" &&
		file != "__SIZERCMD_TIMERDONE_") {
		std::cout << "-E-: removeFile: You can't use this to remove any files other than __SIZERCMD_TIMERERROR_ and __SIZERCMD_TIMERDONE_" << std::endl;
		assert(false);
	}
	std::string filename = dir + "/" + file;
	if (doesFileExist(filename) && remove(filename.c_str())) {
		perror(("-E- removeFile: could not remove '" + filename + "'").c_str());
		return false;
	}
	return true;
}

// Write sizes to a file for timing analysis call
bool TimerInterface::writeSizesForTimer(const std::vector<std::pair<std::string, std::string> > &sizes, const std::string &contest_root, const std::string &benchmark) {
	const std::string filename = contest_root + "/" + benchmark + "/" + benchmark + ".int.sizes";
	std::ofstream ofile(filename.c_str());
	if (!ofile) {
		std::cout << "-E- could not open file '" << filename << "' for output" << std::endl;
		return false;
	}
	for (unsigned i=0; i<sizes.size(); ++i)
		ofile << sizes[i].first << " " << sizes[i].second << std::endl;
	ofile.close();
	return true;
}

// Start timing analysis (does not wait for it to finish)
bool TimerInterface::startTimingAnalysis(const std::string &contest_root, const std::string &benchmark) {
	const std::string filename = contest_root + "/" + benchmark + "/__TCMD_RUNTIMER_";
	if (doesFileExist(filename)) {
		return false;
	}
	
	// Delete previous status files
	removeFile(contest_root+"/"+benchmark,"__SIZERCMD_TIMERERROR_");
	removeFile(contest_root+"/"+benchmark,"__SIZERCMD_TIMERDONE_");
	
	// Instruct the timer to start timing analysis
	std::ofstream ofile(filename.c_str());
	if (!ofile) {    
		std::cout << "-E- startTimingAnalysis: problem instructing timer to run timing, could not write out '" << filename << "'" << std::endl;  
		assert(false);
	}
	ofile.close();
	
	return true;
}

// Run timing analysis in blocking mode
TimerInterface::Status TimerInterface::runTimingAnalysisBlocking(const std::string &contest_root, const std::string &benchmark, const unsigned pollingTime) {
	// Write a file out to instruct timer loop to run timing analysis
	if (!startTimingAnalysis(contest_root,benchmark))
		return TIMER_INTERFACEERROR;
	
	// Wait till timer is done
	Status status = getTimerStatus(contest_root,benchmark);
	while (status == TIMER_BUSY) {
		status = getTimerStatus(contest_root, benchmark);
		wait(pollingTime);
	}
	return status;
}

// Start timing analysis in non-blocking mode
TimerInterface::Status TimerInterface::startTimingAnalysisNonBlocking(const std::string &contest_root, const std::string &benchmark) {
	if (!startTimingAnalysis(contest_root,benchmark))
		return TIMER_INTERFACEERROR;
	return TIMER_BUSY;
}

// END PRIVATE SECTION ----------------------------------------------------------------------------

// Function to pretty-print Status
std::ostream& operator<<(std::ostream &o, const TimerInterface::Status &s) {
	switch (s) {
		case TimerInterface::TIMER_NOT_STARTED:
			o << "TIMER_NOT_STARTED"; break;
		case TimerInterface::TIMER_BUSY:
			o << "TIMER_BUSY"; break;
		case TimerInterface::TIMER_FINISHED_SUCCESS:
			o << "TIMER_FINSIHED_SUCCESS"; break;
		case TimerInterface::TIMER_FINISHED_ERROR:
			o << "TIMER_FINISHED_ERROR"; break;
		case TimerInterface::TIMER_INTERFACEERROR:
			o << "TIMER_INTERFACEERROR"; break;
		default:
			break;
	}
	return o;
}
