#ifndef STEPWATCH_H
#define	STEPWATCH_H

#include <iostream>
using std::cout;
using std::endl;

#include "Stopwatch.h"

class Stepwatch {
private:
	Stopwatch clsStopwatch;
	const string clsMsg;
	
	void printMsg() const {
		//for ( int i = 0; i < clsDepth; i++ )
		//	cout << "  ";
		std::cout << clsMsg << "...";
	} // end method
	
public:
	Stepwatch(const string &msg) : clsMsg(msg) {
		clsStopwatch.start();
		printMsg();
		std::cout << std::endl;
	} // end constructor
	
	void finish() {
		clsStopwatch.stop();
		printMsg();
		std::cout << " runtime: " 
				<< clsStopwatch.getElapsedTime() << " s" << std::endl;
	} // end method
	
	~Stepwatch() {
		if ( clsStopwatch.isRunning() )
			finish();
	} // end constructor
}; // end class

#endif

