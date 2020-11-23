//////////////////////////////////////////////////////////////////
//
//
//  Helper functions and classes to parse the ISPD 2012 contest
//  benchmark files.
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

#ifndef _PARSER_HELPER_H
#define _PARSER_HELPER_H

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <limits>

using std::cout ;
using std::endl ;
using std::istream ;
using std::ostream ;
using std::vector ;
using std::string ;


/////////////////////////////////////////////////////////////////////
//
// This class can be used to parse the specific verilog
// format as defined in the ISPD-13 contest benchmarks. It is not
// intended to be used as a generic verilog parser.
//
// See test_verilog_parser () function in parser_helper.cpp for an
// example of how to use this class.
//
/////////////////////////////////////////////////////////////////////
class VerilogParser {
    
    std::ifstream is ;
    
public:
    
    
    // Constructor that opens the given filename
    VerilogParser (string filename): is(filename.c_str()) {}
    
    // The following functions must be issued in a particular order
    // See test_verilog_parser function for an example
    
    // Read the module definition
    bool read_module (string& moduleName) ;
    
    // Read the next primary input.
    // Return value indicates if the last read was successful or not.
    bool read_primary_input (string& primaryInput) ;
    
    // Read the next primary output.
    // Return value indicates if the last read was successful or not.
    bool read_primary_output (string& primaryInput) ;
    
    
    // Read the next net.
    // Return value indicates if the last read was successful or not.
    bool read_wire (string& wire) ;
    
    // Read the next cell instance.
    // Return value indicates if the last read was successful or not.
    bool read_cell_inst (string& cellType, string& cellInstName,
                         vector<std::pair<string, string> >& pinNetPairs) ;
    
    
} ;


/////////////////////////////////////////////////////////////////////
//
// This class can be used to parse the specific sdc
// format as defined in the ISPD-13 contest benchmarks. It is not
// intended to be used as a generic sdc parser.
//
// See test_sdc_parser () function in parser_helper.cpp for an
// example of how to use this class.
//
/////////////////////////////////////////////////////////////////////
class SdcParser {
    
    std::ifstream is ;
    
public:
    
    
    // Constructor that opens the given filename
    SdcParser (string filename): is(filename.c_str()) {}
    
    // The following functions must be issued in a particular order
    // See test_sdc_parser function for an example
    
    // Read clock definition
    // Return value indicates if the last read was successful or not.
    bool read_clock (string& clockName, string& clockPort, double& period) ;
    
    // Read input delay
    // Return value indicates if the last read was successful or not.
    bool read_input_delay (string& portName, double& delay) ;
    
    // Read driver info for the input port
    // Return value indicates if the last read was successful or not.
    bool read_driver_info (string& inPortName, string& driverSize, string& driverPin,
                           double& inputTransitionFall, double& inputTransitionRise) ;
    
    // Read output delay
    // Return value indicates if the last read was successful or not.
    bool read_output_delay (string& portName, double& delay) ;
    
    // Read output load
    // Return value indicates if the last read was successful or not.
    bool read_output_load (string& outPortName, double& load) ;
    
    
} ;



/////////////////////////////////////////////////////////////////////
//
// The following classes can be used to parse the specific spef
// format as defined in the ISPD-13 contest benchmarks. It is not
// intended to be used as a generic spef parser.
//
// See test_spef_parser () function in parser_helper.cpp for an
// example of how to use this class.
//
/////////////////////////////////////////////////////////////////////
struct SpefNodeName {
    string n1 ;
    string n2 ;
    
    // A node in the spef file can be defined in 3 different ways:
    // 1. For the node corresponding to the connection to a port:
    //       nodeName = "portName", i.e. n1 = "portName", n2 = ""
    //
    // 2. For the node corresponding to the connection to a cell pin:
    //       nodeName = "cellName":"pinName", i.e. n1 = "cellName", n2 = "pinName"
    //
    // 3. For an internal node of an RC tree:
    //       nodeName = "netName":"index", i.e. n1 = "netName", n2 = "index"
	
	// Added by Guilherme Flach
	operator string() const {return n1 + ((n2 != "") ? ":" : "") + n2;}
	
} ;

    ostream& operator<< (ostream& os, const SpefNodeName& n) ;

    struct SpefConnection {
        char nodeType ; // either 'P' (port) or 'I' (internal)
        SpefNodeName nodeName ;
        char direction ; // either 'I' (receiver pin) or 'O' (driver pin)
    };
    
    ostream& operator<< (ostream& os, const SpefConnection& c) ;
    
    struct SpefCapacitance {
        SpefNodeName nodeName ;
        double capacitance ;
        
    } ;
    
    ostream& operator<< (ostream& os, const SpefCapacitance& c) ;

    struct SpefResistance {
        SpefNodeName fromNodeName ;
        SpefNodeName toNodeName ;
        double resistance ;
    } ;
    
    ostream& operator<< (ostream& os, const SpefResistance& r) ;

    struct SpefNet {
        string netName ;
        double netLumpedCap ;
        vector<SpefConnection> connections ;
        vector<SpefCapacitance> capacitances ;
        vector<SpefResistance> resistances ;
        
        void clear() {
            netName = "" ;
            netLumpedCap = 0.0 ;
            connections.clear() ;
            capacitances.clear() ;
            resistances.clear() ;
        }
        
    } ;
    
    class SpefParser {
        
        std::ifstream is ;
        
        bool read_connections (vector<SpefConnection>& connections) ;
        void read_capacitances (vector<SpefCapacitance>& capacitances) ;
        void read_resistances (vector<SpefResistance>& resistances) ;
        
    public:
        
        SpefParser (string filename): is(filename.c_str()) {}
        
        // Read the spef data for the next net.
        // Return value indicates if the last read was successful or not.
        bool read_net_data (SpefNet& spefNet) ;
        
        
    } ;
    
    
    /////////////////////////////////////////////////////////////////////
    //
    // This class can be used to parse the specific .timing
    // format as defined in the ISPD-13 contest benchmarks.
    //
    // See test_timing_parser () function in parser_helper.cpp for an
    // example of how to use this class.
    //
    /////////////////////////////////////////////////////////////////////
    class TimingParser {
        
        std::ifstream is ;
        
    public:
        
        TimingParser (string filename): is(filename.c_str()) {}
        
        // Read timing info for the next pin or port
        // Return value indicates if the last read was successful or not.
        // If the line read corresponds to a pin, then name1 and name2 will be set to the cell
        // instance name and the pin name, respectively.
        // If the line read corresponds to a port, then name1 will be set to the port name, and
        // name2 will be set to "".
        bool read_timing_line (string& name1, string& name2, double& riseSlack, double& fallSlack,
                               double& riseTransition, double& fallTransition,
                               double& riseArrival, double& fallArrival) ;
        
        
    } ;
    
    
    
    /////////////////////////////////////////////////////////////////////
    //
    // This class can be used to parse the specific .ceff
    // format as defined in the ISPD-13 contest benchmarks.
    //
    // See test_ceff_parser () function in parser_helper.cpp for an
    // example of how to use this class.
    //
    /////////////////////////////////////////////////////////////////////
    class CeffParser {
        
        std::ifstream is ;
        
    public:
        
        CeffParser (string filename): is(filename.c_str()) {}
        
        // Read ceff values for the next pin or port
        // Return value indicates if the last read was successful or not.
        // If the line read corresponds to a pin, then name1 and name2 will be set to the cell
        // instance name and the pin name, respectively.
        // If the line read corresponds to a port, then name1 will be set to the port name, and
        // name2 will be set to "".
        bool read_ceff_line (string& name1, string& name2, double& riseCeff, double& fallCeff) ;
        
        
    } ;
    
    
    
    
    
    /////////////////////////////////////////////////////////////////////
    //
    // The following classes can be used to parse the specific lib
    // format as defined in the ISPD-13 contest benchmarks. They are not
    // intended to be used as a generic lib parser.
    //
    // See test_lib_parser () function in parser_helper.cpp for an
    // example of how to use these classes.
    //
    /////////////////////////////////////////////////////////////////////
    
    // Look up table to store delay or slew functions
    struct LibParserLUT {
        
        // Look up table is indexed by the output load and the input transition values
        // Example:
        //   Let L = loadIndices[i]
        //       T = transitionIndices[j]
        //   Then, the table value corresponding to L and T will be:
        //       table[i][j]
        //
        vector<double> loadIndices ;
        vector<double> transitionIndices ;
        vector<vector<double> > tableVals ;
        
    } ;
    
    ostream& operator<< (ostream& os, LibParserLUT& lut) ;
    
    struct LibParserTimingInfo {
        
        string fromPin ;
        string toPin ;
        string timingSense ; // "non_unate" or "negative_unate" or "positive_unate".
        // Note that ISPD-13 library will have only negative-unate combinational cells. The clock arcs
        // for sequentials will be non_unate (which can be ignored because of the simplified sequential
        // timing model for ISPD-13).
        
        
        LibParserLUT fallDelay ;
        LibParserLUT riseDelay ;
        LibParserLUT fallTransition ;
        LibParserLUT riseTransition ;
        
    } ;
    
    ostream& operator<< (ostream& os, LibParserTimingInfo& timing) ;
    
    struct LibParserPinInfo {
        
        string name ; // pin name
        double capacitance ; // input pin cap (not defined for output pins)
        double maxCapacitance ; // the max load this pin can drive
        bool isInput ; // whether the pin is input or output pin
        bool isClock ; // whether the pin is a clock pin or not
        
        LibParserPinInfo () : capacitance (0.0), maxCapacitance (std::numeric_limits<double>::max()),
        isInput(true), isClock(false) {}
        
    } ;
    
    ostream& operator<< (ostream& os, LibParserPinInfo& pin) ;
    
    struct LibParserCellInfo {
        
        string name ; // cell name
        string footprint ; // only the cells with the same footprint are swappable
        double leakagePower ; // cell leakage power
        double area ; // cell area (will not be a metric for ISPD-13)
        bool isSequential ; // if true then sequential cell, else combinational
        bool dontTouch ; // is the sizer allowed to size this cell?
        
        vector<LibParserPinInfo> pins ;
        vector<LibParserTimingInfo> timingArcs ;
        
        LibParserCellInfo () : leakagePower (0.0), area (0.0), isSequential (false), dontTouch(false) {}
        
    } ;
    
    ostream& operator<< (ostream& os, LibParserCellInfo& cell) ;
    
    
    // See test_lib_parser () function in parser_helper.cpp for an
    // example of how to use this class.
    class LibParser {
        
        std::ifstream is ;
        
        void _skip_lut_3D () ;
        void _begin_read_lut (LibParserLUT& lut) ;
        void _begin_read_timing_info (string pinName, LibParserTimingInfo& cell) ;
        void _begin_read_pin_info (string pinName, LibParserCellInfo& cell, LibParserPinInfo& pin) ;
        void _begin_read_cell_info (string cellName, LibParserCellInfo& cell) ;
        
    public:
        
        LibParser (string filename): is(filename.c_str()) {}
        
        // Read the default max_transition defined for the library.
        // Return value indicates if the last read was successful or not.
        // This function must be called in the beginning before any read_cell_info function call.
        bool read_default_max_transition (double& maxTransition) ;
        
        
        // Read the next standard cell definition.
        // Return value indicates if the last read was successful or not.
        bool read_cell_info (LibParserCellInfo& cell) ;
        
        
    } ;
    
#endif
