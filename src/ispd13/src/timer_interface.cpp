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
    } else if ("__TCMD_RUNTIMER_NEG_SLACK_ONLY_" == cmd) {
        return TIMER_BUSY;
    } else if ("__TCMD_RUNTIMER_NO_REPORT_" == cmd) {
        return TIMER_BUSY;
    } else if ("__TCMD_RUNTIMER_TR_" == cmd) {
        return TIMER_BUSY;
    } else if ("__TCMD_RUNTIMER_PR_" == cmd) {
        return TIMER_BUSY;
    } else {
        return TIMER_NOT_STARTED;
    }
}

// Write sizes and run timing analysis in blocking mode
TimerInterface::Status TimerInterface::runTimingAnalysisBlocking(const std::vector<std::pair<std::string, std::string> > &sizes,
                                                                 const bool timing_request, const std::vector<std::string> &timing_pins,
                                                                 const bool ceff_request, const std::vector<std::string> &ceff_pins,
                                                                 const std::string &contest_root, const std::string &benchmark, const unsigned pollingTime) {
    
    const std::string filename = contest_root + "/" + benchmark + "/__TIMER_STARTED_";
    if (!doesFileExist(filename)) {
        return TIMER_INTERFACEERROR;
    }
    
    Status status = getTimerStatus(contest_root,benchmark);
    while (status == TIMER_BUSY) {
        status = getTimerStatus(contest_root, benchmark);
        std::cout << "-I- timer BUSY. Waiting..." << std::endl;
        wait(pollingTime*2);
        status = getTimerStatus(contest_root, benchmark);
    }
    
    if (!(timing_request || ceff_request)) {
        std::cout << "-E- runTimingAnalysisBlocking: both timing_request and ceff_request are false, there is nothing for the timer to do, don't call this function in such cases" << std::endl;
        return TIMER_INTERFACEERROR;
    }
    if (!writeSizesForTimer(sizes, contest_root, benchmark)) {
        std::cout << "-E- runTimingAnalysisBlocking: problem writing sizes" << std::endl;
        return TIMER_INTERFACEERROR;
    }
    if (!writePinsForTimer(timing_request, timing_pins, ceff_request, ceff_pins, contest_root, benchmark)) {
        std::cout << "-E- runTimingAnalysisBlocking: problem writing pins for timing/ceff request" << std::endl;
        return TIMER_INTERFACEERROR;
    }
    return runTimingAnalysisBlocking(contest_root,benchmark,pollingTime);
}

// Write sizes and run timing analysis in blocking mode reporting only negative slack paths in .timing file
TimerInterface::Status TimerInterface::runTimingAnalysisBlockingNegSlackOnly(const std::vector<std::pair<std::string, std::string> > &sizes,
                                                                 const bool timing_request, const std::vector<std::string> &timing_pins,
                                                                 const bool ceff_request, const std::vector<std::string> &ceff_pins,
                                                                 const std::string &contest_root, const std::string &benchmark, const unsigned pollingTime) {
    
    const std::string filename = contest_root + "/" + benchmark + "/__TIMER_STARTED_";
    if (!doesFileExist(filename)) {
        return TIMER_INTERFACEERROR;
    }
    
    Status status = getTimerStatus(contest_root,benchmark);
    while (status == TIMER_BUSY) {
        status = getTimerStatus(contest_root, benchmark);
        std::cout << "-I- timer BUSY. Waiting..." << std::endl;
        wait(pollingTime*2);
        status = getTimerStatus(contest_root, benchmark);
    }
    
    if (!(timing_request || ceff_request)) {
        std::cout << "-E- runTimingAnalysisBlocking: both timing_request and ceff_request are false, there is nothing for the timer to do, don't call this function in such cases" << std::endl;
        return TIMER_INTERFACEERROR;
    }
    if (!writeSizesForTimer(sizes, contest_root, benchmark)) {
        std::cout << "-E- runTimingAnalysisBlocking: problem writing sizes" << std::endl;
        return TIMER_INTERFACEERROR;
    }
    /*if (!writePinsForTimer(timing_request, timing_pins, ceff_request, ceff_pins, contest_root, benchmark)) {
        std::cout << "-E- runTimingAnalysisBlocking: problem writing pins for timing/ceff request" << std::endl;
        return TIMER_INTERFACEERROR;
    }*/
    return runTimingAnalysisBlockingNegSlackOnly(contest_root,benchmark,pollingTime);
}

// Start timing analysis in non-blocking mode
TimerInterface::Status TimerInterface::startTimingAnalysisNonBlocking(const std::vector<std::pair<std::string, std::string> > &sizes,
                                                                      const bool timing_request, const std::vector<std::string> &timing_pins,
                                                                      const bool ceff_request, const std::vector<std::string> &ceff_pins,
                                                                      const std::string &contest_root, const std::string &benchmark) {
    
    const std::string filename = contest_root + "/" + benchmark + "/__TIMER_STARTED_";
    if (!doesFileExist(filename)) {
        return TIMER_INTERFACEERROR;
    }
    
/*
    while (status == TIMER_BUSY) {
        status = getTimerStatus(contest_root, benchmark);
        std::cout << "-I- timer BUSY. Waiting..." << std::endl;
        wait(2);
        status = getTimerStatus(contest_root, benchmark);
    }
 */   
    if (!(timing_request || ceff_request)) {
        std::cout << "-E- startTimingAnalysisNonBlocking: both timing_request and ceff_request are false, there is nothing for the timer to do, don't call this function in such cases" << std::endl;
        return TIMER_INTERFACEERROR;
    }
    if (!writeSizesForTimer(sizes, contest_root, benchmark)) {
        std::cout << "-E- startTimingAnalysisNonBlocking: problem writing sizes" << std::endl;
        return TIMER_INTERFACEERROR;
    }
    if (!writePinsForTimer(timing_request, timing_pins, ceff_request, ceff_pins, contest_root, benchmark)) {
        std::cout << "-E- startTimingAnalysisNonBlocking: problem writing pins for timing/ceff request" << std::endl;
        return TIMER_INTERFACEERROR;
    }

    Status status = getTimerStatus(contest_root,benchmark);
    if (status == TIMER_BUSY) {
        return TIMER_INTERFACEERROR;
    }

    return startTimingAnalysisNonBlocking(contest_root,benchmark);
}

// Start timing analysis in non-blocking mode
TimerInterface::Status TimerInterface::startTimingAnalysisNonBlockingNoReport(const std::vector<std::pair<std::string, std::string> > &sizes,
                                                                      const bool timing_request, const std::vector<std::string> &timing_pins,
                                                                      const bool ceff_request, const std::vector<std::string> &ceff_pins,
                                                                      const std::string &contest_root, const std::string &benchmark) {
    
    const std::string filename = contest_root + "/" + benchmark + "/__TIMER_STARTED_";
    if (!doesFileExist(filename)) {
        return TIMER_INTERFACEERROR;
    }
    
/*
    while (status == TIMER_BUSY) {
        status = getTimerStatus(contest_root, benchmark);
        std::cout << "-I- timer BUSY. Waiting..." << std::endl;
        wait(2);
        status = getTimerStatus(contest_root, benchmark);
    }
 */   
    if (!(timing_request || ceff_request)) {
        std::cout << "-E- startTimingAnalysisNonBlocking: both timing_request and ceff_request are false, there is nothing for the timer to do, don't call this function in such cases" << std::endl;
        return TIMER_INTERFACEERROR;
    }
    if (!writeSizesForTimer(sizes, contest_root, benchmark)) {
        std::cout << "-E- startTimingAnalysisNonBlocking: problem writing sizes" << std::endl;
        return TIMER_INTERFACEERROR;
    }
    /*if (!writePinsForTimer(timing_request, timing_pins, ceff_request, ceff_pins, contest_root, benchmark)) {
        std::cout << "-E- startTimingAnalysisNonBlocking: problem writing pins for timing/ceff request" << std::endl;
        return TIMER_INTERFACEERROR;
    }*/
    Status status = getTimerStatus(contest_root,benchmark);
    if (status == TIMER_BUSY) {
        return TIMER_INTERFACEERROR;
    }

    return startTimingAnalysisNonBlockingNoReport(contest_root,benchmark);
}

// Write sizes and run timing analysis in blocking mode reporting only negative slack paths in .timing file
TimerInterface::Status TimerInterface::runTimingAnalysisBlockingNoReport(const std::vector<std::pair<std::string, std::string> > &sizes,
                                                                             const bool timing_request, const std::vector<std::string> &timing_pins,
                                                                             const bool ceff_request, const std::vector<std::string> &ceff_pins,
                                                                             const std::string &contest_root, const std::string &benchmark, const unsigned pollingTime) {
    
    const std::string filename = contest_root + "/" + benchmark + "/__TIMER_STARTED_";
    if (!doesFileExist(filename)) {
        return TIMER_INTERFACEERROR;
    }
    
    Status status = getTimerStatus(contest_root,benchmark);
    while (status == TIMER_BUSY) {
        status = getTimerStatus(contest_root, benchmark);
        std::cout << "-I- timer BUSY. Waiting..." << std::endl;
        wait(pollingTime*2);
    }
    
    if (!(timing_request || ceff_request)) {
        std::cout << "-E- runTimingAnalysisBlocking: both timing_request and ceff_request are false, there is nothing for the timer to do, don't call this function in such cases" << std::endl;
        return TIMER_INTERFACEERROR;
    }
    if (!writeSizesForTimer(sizes, contest_root, benchmark)) {
        std::cout << "-E- runTimingAnalysisBlocking: problem writing sizes" << std::endl;
        return TIMER_INTERFACEERROR;
    }
    /*if (!writePinsForTimer(timing_request, timing_pins, ceff_request, ceff_pins, contest_root, benchmark)) {
        std::cout << "-E- runTimingAnalysisBlocking: problem writing pins for timing/ceff request" << std::endl;
        return TIMER_INTERFACEERROR;
    }*/
    return runTimingAnalysisBlockingNoReport(contest_root,benchmark,pollingTime);
}

// Write sizes and run timing analysis in blocking mode reporting only negative slack paths in .timing file
TimerInterface::Status TimerInterface::runTimingAnalysisBlockingTR(const std::vector<std::pair<std::string, std::string> > &sizes,
                                                                             const bool timing_request, const std::vector<std::string> &timing_pins,
                                                                             const bool ceff_request, const std::vector<std::string> &ceff_pins,
                                                                             const std::string &contest_root, const std::string &benchmark, const unsigned pollingTime) {
    
    const std::string filename = contest_root + "/" + benchmark + "/__TIMER_STARTED_";
    if (!doesFileExist(filename)) {
        return TIMER_INTERFACEERROR;
    }
    
    Status status = getTimerStatus(contest_root,benchmark);
    while (status == TIMER_BUSY) {
        status = getTimerStatus(contest_root, benchmark);
        std::cout << "-I- timer BUSY. Waiting..." << std::endl;
        wait(pollingTime*2);
    }
    
    if (!(timing_request || ceff_request)) {
        std::cout << "-E- runTimingAnalysisBlocking: both timing_request and ceff_request are false, there is nothing for the timer to do, don't call this function in such cases" << std::endl;
        return TIMER_INTERFACEERROR;
    }
    if (!writeSizesForTimer(sizes, contest_root, benchmark)) {
        std::cout << "-E- runTimingAnalysisBlocking: problem writing sizes" << std::endl;
        return TIMER_INTERFACEERROR;
    }
    /*if (!writePinsForTimer(timing_request, timing_pins, ceff_request, ceff_pins, contest_root, benchmark)) {
        std::cout << "-E- runTimingAnalysisBlocking: problem writing pins for timing/ceff request" << std::endl;
        return TIMER_INTERFACEERROR;
    }*/
    return runTimingAnalysisBlockingTR(contest_root,benchmark,pollingTime);
}

// Write sizes and run timing analysis in blocking mode reporting only negative slack paths in .timing file
TimerInterface::Status TimerInterface::runTimingAnalysisBlockingPR(const std::vector<std::pair<std::string, std::string> > &sizes,
                                                                             const bool timing_request, const std::vector<std::string> &timing_pins,
                                                                             const bool ceff_request, const std::vector<std::string> &ceff_pins,
                                                                             const std::string &contest_root, const std::string &benchmark, const unsigned pollingTime) {
    
    const std::string filename = contest_root + "/" + benchmark + "/__TIMER_STARTED_";
    if (!doesFileExist(filename)) {
        return TIMER_INTERFACEERROR;
    }
    
    Status status = getTimerStatus(contest_root,benchmark);
    while (status == TIMER_BUSY) {
        status = getTimerStatus(contest_root, benchmark);
        std::cout << "-I- timer BUSY. Waiting..." << std::endl;
        wait(pollingTime*2);
    }
    
    if (!(timing_request || ceff_request)) {
        std::cout << "-E- runTimingAnalysisBlocking: both timing_request and ceff_request are false, there is nothing for the timer to do, don't call this function in such cases" << std::endl;
        return TIMER_INTERFACEERROR;
    }
    if (!writeSizesForTimer(sizes, contest_root, benchmark)) {
        std::cout << "-E- runTimingAnalysisBlocking: problem writing sizes" << std::endl;
        return TIMER_INTERFACEERROR;
    }
    /*if (!writePinsForTimer(timing_request, timing_pins, ceff_request, ceff_pins, contest_root, benchmark)) {
        std::cout << "-E- runTimingAnalysisBlocking: problem writing pins for timing/ceff request" << std::endl;
        return TIMER_INTERFACEERROR;
    }*/
    return runTimingAnalysisBlockingPR(contest_root,benchmark,pollingTime);
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
            "__TCMD_RUNTIMER_NEG_SLACK_ONLY_" == files[i] ||
            "__TCMD_RUNTIMER_NO_REPORT_" == files[i] ||
            "__TCMD_RUNTIMER_TR_" == files[i] ||
            "__TCMD_RUNTIMER_PR_" == files[i] ||
            "__TCMD_RUNTIMER_" == files[i]) {
            if (cmd != "") {
                std::cout << "-Error- getTimerStatusString: multiple status found" << std::endl;
                for (unsigned j=0; j<files.size(); ++j) {
                    if ("__SIZERCMD_TIMERERROR_" == files[j] ||
                        "__SIZERCMD_TIMERDONE_" == files[j] ||
                        "__TCMD_RUNTIMER_NEG_SLACK_ONLY_" == files[j] ||
                        "__TCMD_RUNTIMER_NO_REPORT_" == files[j] ||
                        "__TCMD_RUNTIMER_TR_" == files[j] ||
                        "__TCMD_RUNTIMER_PR_" == files[j] ||
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
        std::cout << "-E- removeFile: You can't use this to remove any files other than __SIZERCMD_TIMERERROR_ and __SIZERCMD_TIMERDONE_" << std::endl;
        assert(false);
    }
    std::string filename = dir + "/" + file;
    if (doesFileExist(filename) && remove(filename.c_str())) {
        std::cout << "-E- removeFile: could not remove '" << filename << "'" << std::endl;
        return false;
    }
    return true;
}

// Write sizes to a file for timing analysis call
bool TimerInterface::writeSizesForTimer(const std::vector<std::pair<std::string, std::string> > &sizes, const std::string &contest_root, const std::string &benchmark) {
    const std::string filename = contest_root + "/" + benchmark + "/" + benchmark + ".int.sizes";
    std::ofstream ofile(filename.c_str());
    if (!ofile) {
        std::cout << "-E- writeSizesForTimer: could not open file '" << filename << "' for output" << std::endl;
        return false;
    }
    for (unsigned i=0; i<sizes.size(); ++i)
        ofile << sizes[i].first << " " << sizes[i].second << std::endl;
    ofile.close();
    return true;
}

// Write pins (timing and ceff) to a file for timing analysis call
bool TimerInterface::writePinsForTimer(const bool timing_request, const std::vector<std::string> &timing_pins,
                                       const bool ceff_request, const std::vector<std::string> &ceff_pins,
                                       const std::string &contest_root, const std::string &benchmark) {
    if (!(timing_request || ceff_request)) {
        std::cout << "-E- writePinsForTimer: both timing_request and ceff_request are false" << std::endl;
        return false;
    }
    
    // Timing
    const std::string timing_pins_filename = contest_root + "/" + benchmark + "/" + benchmark + ".timing_pins";
    if (timing_request) {
        std::ofstream ofile(timing_pins_filename.c_str());
        if (!ofile) {
            std::cout << "-E- writePinsForTimer: could not open file '" << timing_pins_filename << "' for output" << std::endl;
            return false;
        }
        for (unsigned i=0; i<timing_pins.size(); ++i)
            ofile << timing_pins[i] << std::endl;
        ofile.close();
    } else {
        if (doesFileExist(timing_pins_filename) && remove(timing_pins_filename.c_str())) {
            std::cout << "-E- writePinsForTimer: could not remove '" << timing_pins_filename << "'" << std::endl;
            return false;
        }
    }
    
    // Ceff
    const std::string ceff_pins_filename   = contest_root + "/" + benchmark + "/" + benchmark + ".ceff_pins";
    if (ceff_request) {
        std::ofstream ofile(ceff_pins_filename.c_str());
        if (!ofile) {
            std::cout << "-E- writePinsForTimer: could not open file '" << ceff_pins_filename << "' for output" << std::endl;
            return false;
        }
        for (unsigned i=0; i<ceff_pins.size(); ++i)
            ofile << ceff_pins[i] << std::endl;
        ofile.close();
    } else {
        if (doesFileExist(ceff_pins_filename) && remove(ceff_pins_filename.c_str())) {
            std::cout << "-E- writePinsForTimer: could not remove '" << ceff_pins_filename << "'" << std::endl;
            return false;
        }
    }
    
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

// Start timing analysis (does not wait for it to finish)
bool TimerInterface::startTimingAnalysisNegSlackOnly(const std::string &contest_root, const std::string &benchmark) {
    const std::string filename = contest_root + "/" + benchmark + "/__TCMD_RUNTIMER_NEG_SLACK_ONLY_";
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

// Start timing analysis (does not wait for it to finish)
bool TimerInterface::startTimingAnalysisNoReport(const std::string &contest_root, const std::string &benchmark) {
    const std::string filename = contest_root + "/" + benchmark + "/__TCMD_RUNTIMER_NO_REPORT_";
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

// Start timing analysis (does not wait for it to finish)
bool TimerInterface::startTimingAnalysisTR(const std::string &contest_root, const std::string &benchmark) {
    const std::string filename = contest_root + "/" + benchmark + "/__TCMD_RUNTIMER_TR_";
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

// Start timing analysis (does not wait for it to finish)
bool TimerInterface::startTimingAnalysisPR(const std::string &contest_root, const std::string &benchmark) {
    const std::string filename = contest_root + "/" + benchmark + "/__TCMD_RUNTIMER_PR_";
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
TimerInterface::Status TimerInterface::runTimingAnalysisBlockingNoReport(const std::string &contest_root, const std::string &benchmark, const unsigned pollingTime) {
    
    // Write a file out to instruct timer loop to run timing analysis
    if (!startTimingAnalysisNoReport(contest_root,benchmark))
        return TIMER_INTERFACEERROR;
    
    // Wait till timer is done
    Status status = getTimerStatus(contest_root,benchmark);
    while (status == TIMER_BUSY) {
        status = getTimerStatus(contest_root, benchmark);
        wait(pollingTime);
    }
    return status;
}

// Run timing analysis in blocking mode
TimerInterface::Status TimerInterface::runTimingAnalysisBlockingTR(const std::string &contest_root, const std::string &benchmark, const unsigned pollingTime) {
    
    // Write a file out to instruct timer loop to run timing analysis
    if (!startTimingAnalysisTR(contest_root,benchmark))
        return TIMER_INTERFACEERROR;
    
    // Wait till timer is done
    Status status = getTimerStatus(contest_root,benchmark);
    while (status == TIMER_BUSY) {
        status = getTimerStatus(contest_root, benchmark);
        wait(pollingTime);
    }
    return status;
}

// Run timing analysis in blocking mode
TimerInterface::Status TimerInterface::runTimingAnalysisBlockingPR(const std::string &contest_root, const std::string &benchmark, const unsigned pollingTime) {
    
    // Write a file out to instruct timer loop to run timing analysis
    if (!startTimingAnalysisPR(contest_root,benchmark))
        return TIMER_INTERFACEERROR;
    
    // Wait till timer is done
    Status status = getTimerStatus(contest_root,benchmark);
    while (status == TIMER_BUSY) {
        status = getTimerStatus(contest_root, benchmark);
        wait(pollingTime);
    }
    return status;
}

// Run timing analysis in blocking mode
TimerInterface::Status TimerInterface::runTimingAnalysisBlockingNegSlackOnly(const std::string &contest_root, const std::string &benchmark, const unsigned pollingTime) {
    
    // Write a file out to instruct timer loop to run timing analysis
    if (!startTimingAnalysisNegSlackOnly(contest_root,benchmark))
        return TIMER_INTERFACEERROR;
    
    Status status = getTimerStatus(contest_root,benchmark);
    // Wait till timer is done
    while (status == TIMER_BUSY) {
        status = getTimerStatus(contest_root, benchmark);
        wait(pollingTime);
    }
    return status;
}

// Run timing analysis in blocking mode
TimerInterface::Status TimerInterface::runTimingAnalysisBlocking(const std::string &contest_root, const std::string &benchmark, const unsigned pollingTime) {
    // Write a file out to instruct timer loop to run timing analysis
    if (!startTimingAnalysis(contest_root,benchmark))
        return TIMER_INTERFACEERROR;
    
    Status status = getTimerStatus(contest_root,benchmark);
    // Wait till timer is done
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

// Start timing analysis in non-blocking mode
TimerInterface::Status TimerInterface::startTimingAnalysisNonBlockingNoReport(const std::string &contest_root, const std::string &benchmark) {
    
    if (!startTimingAnalysisNoReport(contest_root,benchmark))
        return TIMER_INTERFACEERROR;
    return TIMER_BUSY;
}


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


// END PRIVATE SECTION ----------------------------------------------------------------------------
