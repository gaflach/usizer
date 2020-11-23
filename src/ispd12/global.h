/*
 *  global.h
 *  sizer
 *
 *  Created by Tiago Reimann on 11/01/12.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef _GLOBAL_H_
#define _GLOBAL_H_
#include <vector>
#include <set>
#include <map>
#include <queue>
#include <string>
#include "parser_helper.h"
#include "Vcell.h"

using namespace std;

class comp_cellLeakage {//make heap
public:
	bool operator()( const LibParserCellInfo c1, const LibParserCellInfo c2) const
	{
		if ((c2.leakagePower > c1.leakagePower) || ((c2.leakagePower == c1.leakagePower) && (c2.pins[0].maxCapacitance < c1.pins[0].maxCapacitance)))
            return true;
		else return false;
	}
};

class OrdCells {
public:
	string footprint;
	double minimumFallSlew;		//minimum slew value for each footprint (type of cell)
	double minimumRiseSlew;
	
	double maximumFallSlew;		//maximum slew value for each footprint (type of cell)
	double maximumRiseSlew;
	
	vector< LibParserCellInfo > cells;
	
	int findCellTypeIndex(string instType);
};

class OrgCells {
public:
	vector< OrdCells > oCells;
	
	LibParserCellInfo* findCellInst(string instType);
	int findFootprintIndex(string footprint);
};

class Wire {
public:
	string wire_name;
	double cap;
	bool operator<(const Wire& w1) const
	{
		return (wire_name > w1.wire_name);
	}
	bool operator==(const Wire& w1) const
	{
		return (wire_name == w1.wire_name);
	}
	bool operator>(const Wire& w1) const
	{
		return (wire_name < w1.wire_name);
	}
	bool operator!=(const Wire& w1) const
	{
		return (wire_name != w1.wire_name);
	}
};

class InputDelay {
public:
	string port_name;
	double delay;
};

class OutputDelay {
public:
	string port_name;
	double delay;
};

class InputDriver {
public:
	string port_name;
	string driver;
	double rise,fall;
};

class OutputLoad {
public:
	string port_name;
	double load;
};

class SDCInfo {
public:
	string clk_name;
	string clk_port;
	double clk_period;
	vector < InputDelay > inputDelays;
	vector < InputDriver > input_drivers;
	vector < OutputDelay > output_delays;
	vector < OutputLoad > output_loads;
};

class Net {
public:
	string netName;
	set< Vcell* > * cells;
	int numPins;
	bool operator<(const Net& w1) const
	{
		return (netName < w1.netName);
	}
	bool operator==(const Net& w1) const
	{
		return (netName == w1.netName);
	}
	bool operator>(const Net& w1) const
	{
		return (netName > w1.netName);
	}
	bool operator!=(const Net& w1) const
	{
		return (netName != w1.netName);
	}
};

class AddrCell{
public:
	
	string instName;
	int vectorIndex;
	
	AddrCell () : vectorIndex (0) {}
	
	bool operator<(const AddrCell& w1) const
	{
		return (instName < w1.instName);
	}
	bool operator==(const AddrCell& w1) const
	{
		return (instName == w1.instName);
	}
	bool operator>(const AddrCell& w1) const
	{
		return (instName > w1.instName);
	}
	bool operator!=(const AddrCell& w1) const
	{
		return (instName != w1.instName);
	}
};

// -----------------------------------------------------------------------------

#define debug(tag, x) if (app.hasOptionValuePair("tag",tag)) cout << x

class App {
private:
	typedef map< string, set<string> > OptionMap;
	static OptionMap clsOptions;
	
public:
	
	// Returns the number of occurrence of the option specified by key.
	static int hasOption( const string &key ) {
		return clsOptions.count(key);
	} // end method
	
	// Returns true if there are an option specified by key, and one of its
	// values matches value.
	static bool hasOptionValuePair( const string &key, const string &value ) {
		OptionMap::iterator it = clsOptions.find(key);
		if ( it != clsOptions.end() )
			return it->second.count(value);
		else
			return false;
	} // end method
	
	// Return the ith value defined for an option specified by key. If no option
	// exists or an invalid index is used, an empty string is returned.
	static string getOptionValue( const string &key, const int index = 0 ) {
		OptionMap::iterator it = clsOptions.find(key);
		if ( it != clsOptions.end() ) {
			set<string> &optionParameters = it->second;
			
			int counter = 0;
			for (set<string >::iterator its = optionParameters.begin(); its!=optionParameters.end(); ++its) {
				if ( counter++ == index )
					return *its;
			} // end for
		} // end if
		return "";
	} // end method
	
	// Parse the command line arguments.
	static void parseCommandLineArguments( const int argc, char ** argv, int offset ) {
		string currentOption;
		
		for ( int r=offset; r < argc; r++ ) {
			const string arg( argv[r] );

			if ( arg[0] == '-' ) {
				currentOption = arg.substr(1);
				
				OptionMap::iterator it = clsOptions.find(currentOption);
				if ( it == clsOptions.end() )
					clsOptions.insert( make_pair(currentOption, set<string>() ));
			} else {
				OptionMap::iterator it = clsOptions.find(currentOption);
				if ( it != clsOptions.end() )
					it->second.insert(arg);
			} // end else
		} // end for		
	} // end method
	
	
	// Print options.
	static void printCommandLineArguments( ostream &out ) {
		out << "Command line arguments:\n";
		for (OptionMap::iterator it = clsOptions.begin(); it!=clsOptions.end(); ++it) {
			out << "\t" << it->first << "\t[";
			
			bool first = true;
			set<string> &optionParameters = it->second;
			for (set<string >::iterator its = optionParameters.begin(); its!=optionParameters.end(); ++its) {
				if ( !first )
					out << " ";
				out << (*its);
			} // end for
			out << "]\n";
		} // end for
	} // end method
}; // end class

extern App app;

#endif //_GLOBAL_H_
