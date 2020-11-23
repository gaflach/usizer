#ifndef RCTREE_H
#define	RCTREE_H

#include <cassert>

#include <vector>
using std::vector;
#include <map>
using std::map;
#include <queue>
using std::queue;

#include "EdgeArray.h"
#include "NewtonRaphson.h"

class Newton_functor {
public:
    /// This is the hidden variable for class Newton_functor, which makes this class a functor, as opposed to a group of function pointers.
    double Rd;
    double R;
    double C1;
    double C2;
    double dt;
    
    /// the creator for Newton_functor which assigns the hidden variable offset
    Newton_functor( double Rd, double R, double C1, double C2, double dt  ) : Rd(Rd), R(R), C1(C1), C2(C2), dt(dt) { }
    
    /// the function to be solved
    double f( double Ceff) {
        //2.282849659*10^(-22)-6.00*10^(-10)*Ceff+9*Ceff^2*(1-exp(-6.666666667*10^(-11)/Ceff))
        //return 2.282849659e-22 - (6.00e-10 * x) + 9*x*x*(1-exp(-6.666666667e-11/x));
        
        double z = (C2+C1)/(C2*R*C1);
        double a = Rd*R*C1*C2;
        double b = (C1 + C2)*Rd + R*C1;
        double inRoot = b*b - 4.0*a;
        double p1 = (b + sqrt(inRoot))/(2.0*a);
        double p2 = (b - sqrt(inRoot))/(2.0*a);
        
        double A = z/(p1*p2);
        double B = (z-p1)/(p1*(p1-p2));
        double D = (z-p2)/(p2*(p2-p1));
        
        double RdCeff = Rd*Ceff;
        double p = dt/(RdCeff);
        double Rddt = Rd*dt;
        
        return A*dt + (B/p1)*(1-exp(-p1*dt)) + (D/p2)*(1-exp(-p2*dt)) - RdCeff*dt + RdCeff*RdCeff*(1-exp(-p));
    }
    
    /// the derivative of the function to be solved
    double df( double Ceff) {
        double RdCeff = Rd*Ceff;
        double p = dt/(RdCeff);
        double Rddt = Rd*dt;
        return -Rddt + 2.0*Rd*RdCeff*(1-exp(-p)) - Rddt*exp(-p);
    }
};

class RCTreeDescriptor {
public:
	struct Node {
		Node() {
			totalCap = 0;
		} // end constructor
		
		string propName;
		vector<int> propResistors;
		double totalCap;
		string propTag;
	}; // end struct
	
	struct Resistor {
		Resistor() {
			propNode0 = -1;
			propNode1 = -1;
			propValue = 0;
		} // end constructor
		
		int propNode0;
		int propNode1;
		double propValue;
		
		int getOtherNode(const int n) const {
			if ( n == propNode0 )
				return propNode1;
			else if ( n == propNode1 )
				return propNode0;
			else
				return -1;
		} // end method
		
	}; // end struct
    
	struct Capacitor {
		Capacitor() {
			propNode = -1;
			propValue = 0;
		} // end constructor
		
		int propNode;
		double propValue;
	}; // end struct
    
private:
	map<string, int> clsNodeMap;
	
	vector<Node> clsNodes;
	vector<Resistor> clsResistors;
	vector<Capacitor> clsCapacitors;
	
	double clsTotalTreeCapacitance;
	
	int createNode(const string &name) {
		map<string,int>::iterator it = clsNodeMap.find(name);
		if ( it != clsNodeMap.end() ) {
			return it->second;
		} else {
			const int index = clsNodes.size();
			clsNodes.resize(clsNodes.size()+1);
			clsNodes.back().propName = name;
			clsNodeMap[name] = index;
			return index;
		} // end if
	} // end method
	
public:
	
	RCTreeDescriptor() {
		clsTotalTreeCapacitance = 0;
	} // end constructor
	
	void addResistor( const string &sourceNode, const string &targetNode, const double resistance) {
		const int index = clsResistors.size();
		clsResistors.resize(clsResistors.size()+1);
		
		Resistor &r = clsResistors.back();
		r.propNode0 = createNode(sourceNode);
		r.propNode1 = createNode(targetNode);
		r.propValue = resistance;
		
		clsNodes[r.propNode0].propResistors.push_back(index);
		clsNodes[r.propNode1].propResistors.push_back(index);
	} // end method
	
	void addCapacitor( const string &node, const double capacitance ) {
		clsCapacitors.resize(clsCapacitors.size()+1);
		
		Capacitor &c = clsCapacitors.back();
		c.propNode = createNode(node);
		c.propValue = capacitance;
		
		clsNodes[c.propNode].totalCap += c.propValue;
		
		clsTotalTreeCapacitance += c.propValue;
	} // end method
	
	double getTotalTreeCapacitance() const { return clsTotalTreeCapacitance; }
	
	int getNumNodes() const { return clsNodes.size(); }
	int getNumResistors() const { return clsResistors.size(); }
	int getNumCapacitors() const { return clsCapacitors.size(); }
	
	const Node &getNode(const int index) const { return clsNodes[index]; }
	const Resistor &getResistor(const int index) const { return clsResistors[index]; }
	const Capacitor &getCapacitor(const int index) const { return clsCapacitors[index]; }
	
	int findNode( const string &name ) const {
		map<string, int>::const_iterator it = clsNodeMap.find(name);
		return ( it != clsNodeMap.end()? it->second : -1 );
	} // end method
	
	void setNodeTag( const int nodeIndex, const string &tag ) {
		clsNodes[nodeIndex].propTag = tag;
	} // end method
	
	void applyDefaultCap(const double cap) {
		const int numCaps = clsCapacitors.size();
		for ( int i = 0; i < numCaps; i++) {
			Capacitor &capacitor = clsCapacitors[i];
			if ( capacitor.propValue == 0.0 ) {
				capacitor.propValue += cap;
				clsTotalTreeCapacitance += cap;
			} // end if
		} // end for
	} // end method
	
}; // end class

// =============================================================================
// Original version of RCTree.
// =============================================================================

class RCTreeBase {
	
public:
	struct Node {
		double propCap;
		double propOriginalCap;
		double propDownstreamCap;
		double propDrivingResistance;
		int propParent;
		bool propEndpoint;
		
		EdgeArray<double> propDelay;
		EdgeArray<double> propSlew;
		EdgeArray<double> propEffectiveCap;
		
		double propY1;
		double propY2;
		double propY3;
	};
	
protected:
    
	struct Ref {
		int propNodeIndex;
		int propDrivingResistorIndex;
		
		Ref() {
			propNodeIndex = -1;
			propDrivingResistorIndex = -1;
		}
		
		Ref(const int n, const int r ) {
			propNodeIndex = n;
			propDrivingResistorIndex = r;
		} // end constructor
		
	};
	
	vector<string> clsNodeNames;
	vector<string> clsNodeTags;
	vector<Node> clsNodes;
	vector< EdgeArray<double> > clsCeffs;
	
	bool clsDirty;
	
	// -------------------------------------------------------------------------
	
	// Source: http://java2s.com/Tutorial/Cpp/0040__Data-Types/Testswhethertwofloatingpointnumbersareapproximatelyequal.htm
	static bool nearlyEqual( const double x, const double y, const double precision ) {
		if (x == 0) return fabs(y) <= precision;
		if (y == 0) return fabs(x) <= precision;
		return fabs(x - y) / max(fabs(x), fabs(y)) <= precision;
	} // end method
	
	// -------------------------------------------------------------------------
	
	static double pow2(const double v) { return v*v; }
	static double pow3(const double v) { return v*v*v; }
	
	// -------------------------------------------------------------------------
	
	// Build tree from tree descriptor given a root node.
	void buildTopology( const RCTreeDescriptor &dscp, const int root ) {
		const int numNodes = dscp.getNumNodes();
		const int numResistors = dscp.getNumResistors();
		
		// Compute topological order.
		vector<bool> visited(numNodes, false); // loop detect
		vector<bool> tackled(numResistors, false); // avoid processing resistors more than one time
		
		vector<int> topology(numNodes, -1);
		vector<int> reverseTopology(numNodes, -1);
		
		queue<Ref> q;
        
		const RCTreeDescriptor::Node &rootNodeDescriptor = dscp.getNode(root);
		
		Node &rootNode = clsNodes[0];
		rootNode.propDelay.set(0.0, 0.0);
		rootNode.propCap = rootNodeDescriptor.totalCap;
		rootNode.propOriginalCap = rootNodeDescriptor.totalCap;
		rootNode.propDrivingResistance = -1;
		rootNode.propParent = -1;
		clsNodeNames[0] = rootNodeDescriptor.propName;
		clsNodeTags[0] = rootNodeDescriptor.propTag;
        
		for ( int k = 0; k < rootNodeDescriptor.propResistors.size(); k++ ) {
			const int r = rootNodeDescriptor.propResistors[k];
			tackled[r] = true;
			q.push(Ref(dscp.getResistor(r).getOtherNode(root), r));
		} // end method
		rootNode.propEndpoint = q.empty();
        
		topology[0] = root;
		reverseTopology[root] = 0;
        
		int counter = 1;
		
		while(!q.empty()) {
			const int n = q.front().propNodeIndex;
			const int r = q.front().propDrivingResistorIndex;
			q.pop();
			
			assert(!visited[n]); // loop detected
			visited[n] = true;
			
			const RCTreeDescriptor::Node &nodeDescriptor = dscp.getNode(n);
			const RCTreeDescriptor::Resistor &resistorDescriptor = dscp.getResistor(r);
			
			Node &node = clsNodes[counter];
			node.propCap = nodeDescriptor.totalCap;
			node.propOriginalCap = nodeDescriptor.totalCap;
			node.propDrivingResistance = resistorDescriptor.propValue;
			node.propParent = reverseTopology[resistorDescriptor.getOtherNode(n)];
			clsNodeNames[counter] = nodeDescriptor.propName;
			clsNodeTags[counter] = nodeDescriptor.propTag;
			
			topology[counter] = n;
			reverseTopology[n] = counter;
			
			int counterNeighbours = 0;
			
			const int numResistors = nodeDescriptor.propResistors.size();
			for ( int k = 0; k < numResistors; k++ ) {
				const int r = nodeDescriptor.propResistors[k];
				if ( !tackled[r] ) {
					tackled[r] = true;
					q.push(Ref(dscp.getResistor(r).getOtherNode(n), r));
					counterNeighbours++;
				} // end if
			} // end method
			
			node.propEndpoint = ( counterNeighbours == 0 );
			
			counter++;
		} // end while
	} // end method
	
	// -------------------------------------------------------------------------
    
	template<class RCTreeDriver>
	void stepForward( const RCTreeDriver &driver ) {
		Node &rootState = clsNodes[0];
		
		rootState.propSlew = driver.computeSlew(rootState.propEffectiveCap);
        
		const int numNodes = clsNodes.size();
		for ( int n = 1; n < numNodes; n++ ) { // 1 => skips root node
			Node &node = clsNodes[n];
			const Node &parent = clsNodes[node.propParent];
			
			const EdgeArray<double> S0 = parent.propSlew;
			const EdgeArray<double> Ceff1 = node.propEffectiveCap;
			const double R1 = node.propDrivingResistance;
            
			const EdgeArray<double> RCeff = R1*Ceff1;
			node.propDelay = parent.propDelay + RCeff;
			node.propSlew = S0 / ( 1.0 - ( (RCeff)/S0 )*(1.0 - exp(-S0/(RCeff) ) ) );
		} // end for
	} // end for
	
	// -------------------------------------------------------------------------
	
	void stepBackward() {
		const int numNodes = clsNodes.size();
        
		for ( int i = 0; i < numNodes; i++ ) {
			Node &node = clsNodes[i];
            
			clsCeffs[i] = node.propEffectiveCap;
			node.propEffectiveCap.set(node.propCap, node.propCap);
		} // end for
		
		for ( int n = numNodes - 1; n > 0; n-- ) { // n > 0 skips root node
			const Node &sink = clsNodes[n];
			Node &driver = clsNodes[sink.propParent];
			
            const EdgeArray<double> S0 = driver.propSlew;
			const EdgeArray<double> Ceff1 = clsCeffs[n];
			const double R1 = sink.propDrivingResistance;
			
			const EdgeArray<double> RCeff = R1*Ceff1;
			const EdgeArray<double> K1 = 1.0 - ( (2.0*RCeff)/S0 )*(1.0 - exp(-S0/(2.0*RCeff) ) );
            driver.propEffectiveCap += K1 * sink.propDownstreamCap;
		} // end for
	} // end for
	
public:
	
	void build( const RCTreeDescriptor &dscp, const string &rootNodeName ) {
		const int numNodes = dscp.getNumNodes();
		
		// Clean up
		clsDirty = false;
		
		clsNodes.resize(numNodes);
		clsNodeNames.resize(numNodes);
		clsNodeTags.resize(numNodes);
		clsCeffs.resize(numNodes);
        
		// Map topological index (e.g. 0 = root) to respective node index in the
		// descriptor.
		buildTopology(dscp, dscp.findNode(rootNodeName));
        
		// Update downstream cap.
		updateDownstreamCap();
	} // end method
    
	// -------------------------------------------------------------------------
	
	void updateDownstreamCap() {
		const int numNodes = getNumNodes();
		for ( int i = 0; i < numNodes; i++ ) {
			Node &node = clsNodes[i];
			node.propDownstreamCap = node.propCap;
		} // end for
		
		for ( int n = numNodes - 1; n > 0; n-- ) { // n > 0 skips root node
			const Node &node = clsNodes[n];
			clsNodes[node.propParent].propDownstreamCap += node.propDownstreamCap;
		} // end for
		
		clsDirty = false;
	} // end method
	
	// -------------------------------------------------------------------------
	
	void updateDrivingPoint() {
		const int numNodes = getNumNodes();
		for ( int i = 0; i < numNodes; i++ ) {
			Node &node = clsNodes[i];
			node.propY1 = 0;
			node.propY2 = 0;
			node.propY3 = 0;
		} // end for
		
		// Compute pi-model of the RC Tree.
		for ( int n = numNodes - 1; n > 0; n-- ) { // n > 0 skips root node
			const Node &node = clsNodes[n];
            
			const double C = node.propCap;
			const double R = node.propDrivingResistance;
			
			const double yD1 = node.propY1;
			const double yD2 = node.propY2;
			const double yD3 = node.propY3;
			
			const double yU1 = yD1 + C;
			const double yU2 = yD2 - R * (pow2(yD1) + C*yD1 + (1.0/3.0)*pow2(C));
			const double yU3 = yD3 - R * (2*yD1*yD2 + C*yD2) +
            pow2(R)*(pow3(yD1) + (4.0/3.0)*C*pow2(yD1) + (2.0/3.0)*pow2(C)*yD1 + (2.0/15.0)*pow3(C) );
			
			Node &parent = clsNodes[node.propParent];
			parent.propY1 += yU1;
			parent.propY2 += yU2;
			parent.propY3 += yU3;
			
			//cout << "Resistor: " << clsNodeNames[node.propParent] << " -> " << clsNodeNames[n] << "\n";
		} // end for
		
		Node &root = clsNodes[0];
		root.propY1 += root.propCap;
	} // end method
	
	// -------------------------------------------------------------------------
	
	// [PAPER] Fast and Accurate Wire Delay Estimation for Physical Synthesis
	// of Large ASICs.
	
	template<class RCTreeDriver>
	void simulate( const RCTreeDriver &driver, const double epsilon = 1e-6, const int maxIterations = 100 ) {
		if ( clsDirty )
			updateDownstreamCap();
		
		const int numNodes = clsNodes.size();
        
		// Initially set effective capacitance equals to downstream capacitance.
		for ( int i = 0; i < numNodes; i++ ) {
			Node &node = clsNodes[i];
			node.propEffectiveCap = EdgeArray<double>(clsNodes[i].propDownstreamCap, clsNodes[i].propDownstreamCap);
		} // end for
		
		bool converged = false;
		EdgeArray<double> previousRootEffectiveCapacitance = clsNodes[0].propEffectiveCap;
		for ( int i = 0; i < maxIterations; i++ ) {
			stepForward(driver);
			stepBackward();
			
			const Node &root = clsNodes[0];
			
			if ( nearlyEqual(previousRootEffectiveCapacitance[RISE], root.propEffectiveCap[RISE], epsilon ) &&
                nearlyEqual(previousRootEffectiveCapacitance[FALL], root.propEffectiveCap[FALL], epsilon )	)
			{
				converged = true;
				break;
			} // end if
            
			previousRootEffectiveCapacitance = root.propEffectiveCap;
		} // end for
		
		if ( !converged )
			cout << "[WARNING] Simulation for RC tree driven by node '" << clsNodeNames[0] << "' did not converge within " << maxIterations << " iterations.\n";
	} // end method
	
	// -------------------------------------------------------------------------
	
	// [PAPER] Modeling the Effective Capacitance for the RC Interconnect of
	// CMOS Gates
	
	template<class RCTreeDriver>
	EdgeArray<double> computeEffectiveCapacitanceBasedOnJessica( const RCTreeDriver &driver, const double epsilon = 1e-6, const int maxIterations = 100 ) {
		bool converged = false;
		
		const EdgeArray<double> tr = driver.getInputSlew();
		
		double C2, R, C1;
		reduceToPiModel(C1, R, C2);
		
		EdgeArray<double> Ceff(getLumpedCap(),getLumpedCap());
		EdgeArray<double> previousCeff = Ceff;
		
		for ( int i = 0; i < maxIterations; i ++ ) {
            
			const EdgeArray<double> td = driver.computeDelay(Ceff);
			const EdgeArray<double> tf = driver.computeSlew(Ceff);
            
			const EdgeArray<double> tD = td + tr/2.0;
			const EdgeArray<double> tx = tD - 0.5*tf;
            
			const EdgeArray<double> shielding  = ( 1 -
                                                  (R*C1)/(tD - tx/2.0) +
                                                  (pow(R*C1, 2.0)/(tx*(tD - tx/2.0))) *
                                                  (exp(-(tD - tx)/(R*C1) ) ) *
                                                  (1 - exp( (-tx)/(R*C1) ) )
                                                  );
            
			Ceff = C2 + C1*(1.0*shielding + 1.0)/2.0;
			
			if ( nearlyEqual(previousCeff[RISE], Ceff[RISE], epsilon ) &&
                nearlyEqual(previousCeff[FALL], Ceff[FALL], epsilon )	)
			{
				converged = true;
				break;
			} // end if
            
			previousCeff = Ceff;
			
		} // end for
        /*
         if ( !converged )
         cout << "[WARNING] Simulation for RC tree driven by node '" << clsNodeNames[0] << "' did not converge within " << maxIterations << " iterations.\n";
         */
		return Ceff;
	} // end method
	
	
	// -------------------------------------------------------------------------
	
	// [PAPER] Performance Computation for Prec aracterized
    // CMOS Gates with RC Loads
	
	template<class RCTreeDriver>
	EdgeArray<double> computeEffectiveCapacitanceBasedOnMenezes( const RCTreeDriver &driver, const double epsilon = 1e-6, const int maxIterations = 10 ) {
		bool converged = false;
		
		const EdgeArray<double> tr = driver.getInputSlew();
		
		double C2, R, C1;
		reduceToPiModel(C1, R, C2);
		
		EdgeArray<double> Ceff(getLumpedCap(),getLumpedCap());
		EdgeArray<double> Rd, previousCeff = Ceff;
		double t0, dt;
        
		for ( int i = 0; i < maxIterations; i ++ ) {
            
			const EdgeArray<double> td = driver.computeDelay(Ceff);
			const EdgeArray<double> deltatd = (driver.computeDelay(Ceff) - driver.computeDelay(Ceff*0.99));
			const EdgeArray<double> tf = driver.computeSlew(Ceff);
            
			const EdgeArray<double> tD = td + tr/2;
			const EdgeArray<double> tx = tD - 0.5*tf;
            
            Rd =  deltatd/(Ceff*0.01);
			Newton_functor functorRise( Rd[RISE], R, C1, C2, tf[RISE]);
            // {
            NewtonRaphsonSolve0< Newton_functor, double >
            newtonRise(
                       Ceff[RISE],               // X0
                       1e-30,                 // Epsilon
                       true,                 // Twice_df
                       100,                  // max_iter
                       functorRise,              // F
                       &Newton_functor::f,   // f
                       &Newton_functor::df); // df
            
            //   cout << "running iteration with f and df defined, over all real numbers" << endl;
            newtonRise.set_check_boundary( true);
            newtonRise.set_max_x( C1 + C2 );
            newtonRise.set_min_x( C2 );
            
            Newton_functor functorFall( Rd[FALL], R, C1, C2, tf[FALL]);
            // {
            NewtonRaphsonSolve0< Newton_functor, double >
            newtonFall(
                       Ceff[FALL],               // X0
                       1e-30,                 // Epsilon
                       true,                 // Twice_df
                       100,                  // max_iter
                       functorFall,              // F
                       &Newton_functor::f,   // f
                       &Newton_functor::df); // df
            
            //   cout << "running iteration with f and df defined, over all real numbers" << endl;
            newtonFall.set_check_boundary( true);
            newtonFall.set_max_x( C1 + C2 );
            newtonFall.set_min_x( C2 );
            
            //Ceff = EdgeArray<double>(newtonRise.do_iteration( &cout), newtonFall.do_iteration( &cout)) ;
            Ceff = EdgeArray<double>(newtonRise.do_iteration( NULL/*&cout*/), newtonFall.do_iteration( NULL/*&cout*/)) ;
            //cout << "Ceff: " << Ceff << endl;
            
            if ( nearlyEqual(previousCeff[RISE], Ceff[RISE], epsilon ) &&
                nearlyEqual(previousCeff[FALL], Ceff[FALL], epsilon )	)
			{
				converged = true;
				break;
			} // end if
            
			previousCeff = Ceff;
			
		} // end for
        /*
         if ( !converged )
         cout << "[WARNING] Simulation for RC tree driven by node '" << clsNodeNames[0] << "' did not converge within " << maxIterations << " iterations.\n";
         */
		return Ceff;
	} // end method
	
	
	// -------------------------------------------------------------------------
	
	// [PAPER] Modeling the Driving-Point Characteristic of Resistive
	// Interconnect for Accurate Delay Estimation
	
	void reduceToPiModel(double &C1, double &R, double &C2) {
		if ( clsDirty )
			updateDownstreamCap();
		
		updateDrivingPoint();
		
		const Node &root = clsNodes[0];
		C1 = pow2(root.propY2) / root.propY3;
		C2 = root.propY1 - C1;
		R = -pow2(root.propY3)/pow3(root.propY2);
	} // end method
	
	// -------------------------------------------------------------------------
	
	void setNodeExtraCap(const int index, const double cap) {
		Node &node = clsNodes[index];
		node.propCap = node.propOriginalCap + cap;
		
		clsDirty = true;
	} // end method
	
	// -------------------------------------------------------------------------
    
	int getNumNodes() const { return clsNodes.size(); }
	
	const Node &getNode( const int index ) const { return clsNodes[index]; }
	const Node &getRoot() const { return clsNodes[0]; }
	
	const string &getNodeName( const int index ) const { return clsNodeNames[index]; }
	const string &getNodeTag( const int index ) const { return clsNodeTags[index]; }
	
	double getLumpedCap() const { return clsNodes[0].propDownstreamCap; }
	
}; // end class


// =============================================================================
// Flach's version of RCTree.
// =============================================================================

class RCTreeFlach : public RCTreeBase {
private:
    
	template<class RCTreeDriver>
	void stepForward( const RCTreeDriver &driver ) {
		Node &rootState = clsNodes[0];
        
		rootState.propSlew = driver.computeSlew(rootState.propEffectiveCap);
        
		const int numNodes = clsNodes.size();
		for ( int n = 1; n < numNodes; n++ ) { // 1 => skips root node
			Node &node = clsNodes[n];
			const Node &parent = clsNodes[node.propParent];
            
			const EdgeArray<double> Ceff1 = node.propEffectiveCap;
			const double R1 = node.propDrivingResistance;
            
			const EdgeArray<double> RCeff = R1*Ceff1;
			node.propDelay = parent.propDelay + RCeff;
			node.propSlew = parent.propSlew +  RCeff*1.386294361/2;
		} // end for
	} // end for
    
	// -------------------------------------------------------------------------
    
	void stepBackward() {
		const int numNodes = clsNodes.size();
        
		for ( int i = 0; i < numNodes; i++ ) {
			Node &node = clsNodes[i];
            
			clsCeffs[i] = node.propEffectiveCap;
			node.propEffectiveCap.set(node.propCap, node.propCap);
		} // end for
        
		for ( int n = numNodes - 1; n > 0; n-- ) { // n > 0 skips root node
			const Node &sink = clsNodes[n];
			Node &driver = clsNodes[sink.propParent];
            
			const EdgeArray<double> S0 = driver.propSlew;
			const EdgeArray<double> Ceff1 = clsCeffs[n];
			const double R1 = sink.propDrivingResistance;
            
			const EdgeArray<double> RCeff = R1*Ceff1;
			const EdgeArray<double> K1 = 1.0 - ( (2.0*RCeff)/S0 )*(1.0 - exp(-S0/(2.0*RCeff) ) );
            driver.propEffectiveCap += K1 * sink.propDownstreamCap;
		} // end for
	} // end for
    
    
public:
    
	template<class RCTreeDriver>
	void simulate( const RCTreeDriver &driver, const double epsilon = 1e-6, const int maxIterations = 100 ) {
		if ( clsDirty )
			updateDownstreamCap();
        
		const int numNodes = clsNodes.size();
        
		// Initially set effective capacitance equals to downstream capacitance.
		for ( int i = 0; i < numNodes; i++ ) {
			Node &node = clsNodes[i];
			node.propEffectiveCap = EdgeArray<double>(clsNodes[i].propDownstreamCap, clsNodes[i].propDownstreamCap);
		} // end for
        
		bool converged = false;
		EdgeArray<double> previousRootEffectiveCapacitance = clsNodes[0].propEffectiveCap;
		for ( int i = 0; i < maxIterations; i++ ) {
			stepForward(driver);
			stepBackward();
            
			const Node &root = clsNodes[0];
            
			if ( nearlyEqual(previousRootEffectiveCapacitance[RISE], root.propEffectiveCap[RISE], epsilon ) &&
                nearlyEqual(previousRootEffectiveCapacitance[FALL], root.propEffectiveCap[FALL], epsilon )	)
			{
				converged = true;
				break;
			} // end if
            
			previousRootEffectiveCapacitance = root.propEffectiveCap;
		} // end for
        
		if ( !converged )
			cout << "[WARNING] Simulation for RC tree driven by node '" << clsNodeNames[0] << "' did not converge within " << maxIterations << " iterations.\n";
	} // end method
    
}; // end class

// =============================================================================
// Default's version of RCTree.
// =============================================================================

class RCTreeDefault : public RCTreeBase {
private:
    
	template<class RCTreeDriver>
	void stepForward( const RCTreeDriver &driver ) {
		Node &rootState = clsNodes[0];
        
        const EdgeArray<double> downCap(rootState.propDownstreamCap,rootState.propDownstreamCap);
		rootState.propSlew = 1.05*driver.computeSlew(min(1.05*rootState.propEffectiveCap,downCap));
        
		const int numNodes = clsNodes.size();
		for ( int n = 1; n < numNodes; n++ ) { // 1 => skips root node
			Node &node = clsNodes[n];
			const Node &parent = clsNodes[node.propParent];
            
			const EdgeArray<double> S0 = parent.propSlew;
			const EdgeArray<double> Ceff1 = node.propEffectiveCap;
			const double R1 = node.propDrivingResistance;
            
            const EdgeArray<double> RCeff = R1*Ceff1;
            
			node.propDelay = parent.propDelay + RCeff;
			//node.propSlew = S0 + 2.1*RCeff;
            node.propSlew = S0 + 1.386*RCeff;
            //node.propSlew = S0 / ( 1.0 - ( (RCeff)/S0 )*(1.0 - exp(-S0/(RCeff) ) ) );
			//node.propSlew = S0*max(EdgeArray<double>(1.0,1.0),node.propDelay/S0);
		} // end for
	} // end for
    
	// -------------------------------------------------------------------------
    
	void stepBackward() {
		const int numNodes = clsNodes.size();
        
		for ( int i = 0; i < numNodes; i++ ) {
			Node &node = clsNodes[i];
            
			clsCeffs[i] = node.propEffectiveCap;
			node.propEffectiveCap.set(node.propCap, node.propCap);
		} // end for
        
		for ( int n = numNodes - 1; n > 0; n-- ) { // n > 0 skips root node
			const Node &sink = clsNodes[n];
			Node &driver = clsNodes[sink.propParent];
            
			const EdgeArray<double> S0 = driver.propSlew;
			const EdgeArray<double> Ceff1 = clsCeffs[n];
			const double R1 = sink.propDrivingResistance;
            
			const EdgeArray<double> RCeff = 2.0*R1*Ceff1;
			const EdgeArray<double> K1 = 1.0 - ( (RCeff)/S0 )*(1.0 - exp(-S0/(RCeff) ) );
            //driver.propEffectiveCap += K1 * sink.propDownstreamCap;
			//const EdgeArray<double> K1 = 1.0 - min(EdgeArray<double>(1.0,1.0),(RCeff)/S0);
            driver.propEffectiveCap += K1 * sink.propDownstreamCap;
            
		} // end for
	} // end for
    
public:
    
	template<class RCTreeDriver>
	void simulate( const RCTreeDriver &driver, const double epsilon = 1e-6, const int maxIterations = 100 ) {
		if ( clsDirty )
			updateDownstreamCap();
        
		const int numNodes = clsNodes.size();
        
		// Initially set effective capacitance equals to downstream capacitance.
		for ( int i = 0; i < numNodes; i++ ) {
			Node &node = clsNodes[i];
			node.propEffectiveCap = EdgeArray<double>(clsNodes[i].propDownstreamCap, clsNodes[i].propDownstreamCap);
		} // end for
        
		bool converged = false;
		EdgeArray<double> previousRootEffectiveCapacitance = clsNodes[0].propEffectiveCap;
		for ( int i = 0; i < maxIterations; i++ ) {
			stepForward(driver);
			stepBackward();
            
			const Node &root = clsNodes[0];
            
			if ( nearlyEqual(previousRootEffectiveCapacitance[RISE], root.propEffectiveCap[RISE], epsilon ) &&
                nearlyEqual(previousRootEffectiveCapacitance[FALL], root.propEffectiveCap[FALL], epsilon )	)
			{
				converged = true;
				break;
			} // end if
            
			previousRootEffectiveCapacitance = root.propEffectiveCap;
		} // end for
        
		if ( !converged )
			cout << "[WARNING] Simulation for RC tree driven by node '" << clsNodeNames[0] << "' did not converged within " << maxIterations << " iterations.\n";
	} // end method
    
}; // end class

// =============================================================================
// Reimann's version of RCTree.
// =============================================================================

class RCTreeReimann : public RCTreeBase {
private:
    
	template<class RCTreeDriver>
	void stepForward( const RCTreeDriver &driver ) {
		Node &rootState = clsNodes[0];
		
		const EdgeArray<double> downCap(rootState.propDownstreamCap,rootState.propDownstreamCap);
		/*if (rootState.propEffectiveCap[RISE] < rootState.propDownstreamCap*0.9)
         rootState.propEffectiveCap[RISE] *= 0.75;
         if (rootState.propEffectiveCap[FALL] < rootState.propDownstreamCap*0.9)
         rootState.propEffectiveCap[FALL] *= 0.75;*/
		rootState.propSlew = 0.99*driver.computeSlew(rootState.propEffectiveCap);
        
		const EdgeArray<double> rootSlew = rootState.propSlew;
		const int numNodes = clsNodes.size();
		for ( int n = 1; n < numNodes; n++ ) { // 1 => skips root node
			Node &node = clsNodes[n];
			const Node &parent = clsNodes[node.propParent];
            
			const EdgeArray<double> S0 = parent.propSlew;
			const EdgeArray<double> Ceff1 = node.propEffectiveCap;
			const double R1 = node.propDrivingResistance;
            
            const EdgeArray<double> RCeff = R1*Ceff1;
			const EdgeArray<double> ratio = (RCeff)/S0;
			
            //cout << S0 << Ceff1 << R1 << endl;
            
			node.propDelay = parent.propDelay + RCeff;
			//node.propSlew = S0 + 0.5*RCeff;
            node.propSlew = sqrt(S0*S0 + 1.93*RCeff*RCeff);
            //node.propSlew = ( ( (node.propDelay[RISE] > S0[RISE]*0.2) && (node.propDelay[FALL] > S0[FALL]*0.2) )?(sqrt(S0*S0 + RCeff*RCeff)):S0 );
            //node.propSlew[RISE] =  ( (node.propDelay[RISE] > S0[RISE]*2.0) )?(S0[RISE]/(1.0 - (RCeff[RISE]/S0[RISE])*(1.0 - fmath::expd(-S0[RISE]/RCeff[RISE]) ) ) ):S0[RISE];
			//node.propSlew[FALL] =  ( (node.propDelay[FALL] > S0[FALL]*2.0) )?(S0[FALL]/(1.0 - (RCeff[FALL]/S0[FALL])*(1.0 - fmath::expd(-S0[FALL]/RCeff[FALL]) ) ) ):S0[FALL];
			//node.propSlew = (node.propSlew + 19.0*S0)/20.0;
			//node.propSlew = S0 / ( 1.0 - ( ratio )*(1.0 - exp(-1.0/ratio) ) );
			//node.propSlew = sqrt(rootSlew*rootSlew + node.propDelay*node.propDelay);
			//node.propSlew = S0*max(EdgeArray<double>(1.0,1.0),node.propDelay/S0);
            
		} // end for
	} // end for
    
	// -------------------------------------------------------------------------
    
	void stepBackward() {
		const int numNodes = clsNodes.size();
        
		for ( int i = 0; i < numNodes; i++ ) {
			Node &node = clsNodes[i];
            
			clsCeffs[i] = node.propEffectiveCap;
			node.propEffectiveCap.set(node.propCap, node.propCap);
		} // end for
        
		for ( int n = numNodes - 1; n > 0; n-- ) { // n > 0 skips root node
			const Node &sink = clsNodes[n];
			Node &driver = clsNodes[sink.propParent];
            
			const EdgeArray<double> one(0.999,0.999);
			const EdgeArray<double> S0 = driver.propSlew;
			const EdgeArray<double> Ceff1 = clsCeffs[n];
			const double R1 = sink.propDrivingResistance;
            
			const EdgeArray<double> RCeff = 2.0*R1*Ceff1;
			const EdgeArray<double> ratio = (RCeff)/S0;
			const EdgeArray<double> K1 = 1.0 - ( ratio )*(1.0 - exp(-1.0/ratio) );
            //const EdgeArray<double> K1 = pow(S0/sink.propSlew,3);
            
            //driver.propEffectiveCap += K1 * sink.propDownstreamCap;
			//const EdgeArray<double> K1 = 1.0 - min(EdgeArray<double>(1.0,1.0),(RCeff)/S0);
            driver.propEffectiveCap += sink.propDownstreamCap;
            
		} // end for
	} // end for
    
public:
    
	template<class RCTreeDriver>
	void simulate( const RCTreeDriver &driver, const double epsilon = 1e-6, const int maxIterations = 100 ) {
		if ( clsDirty )
			updateDownstreamCap();
        
		const int numNodes = clsNodes.size();
        
		// Initially set effective capacitance equals to downstream capacitance.
		for ( int i = 0; i < numNodes; i++ ) {
			Node &node = clsNodes[i];
			node.propEffectiveCap = EdgeArray<double>(clsNodes[i].propDownstreamCap, clsNodes[i].propDownstreamCap);
		} // end for
        
		bool converged = false;
		clsNodes[0].propEffectiveCap = computeEffectiveCapacitanceBasedOnJessica(driver, epsilon, maxIterations);
		EdgeArray<double> previousRootEffectiveCapacitance = clsNodes[0].propEffectiveCap;
		for ( int i = 0; i < maxIterations; i++ ) {
			stepForward(driver);
			//stepBackward();
			break;
			const Node &root = clsNodes[0];
            
			if ( nearlyEqual(previousRootEffectiveCapacitance[RISE], root.propEffectiveCap[RISE], epsilon ) &&
                nearlyEqual(previousRootEffectiveCapacitance[FALL], root.propEffectiveCap[FALL], epsilon )	)
			{
				converged = true;
				break;
			} // end if
            
			previousRootEffectiveCapacitance = root.propEffectiveCap;
		} // end for
		/*if (clsNodes[0].propEffectiveCap[RISE] < clsNodes[0].propDownstreamCap*0.9)
         clsNodes[0].propEffectiveCap[RISE] *= 0.75;
         if (clsNodes[0].propEffectiveCap[FALL] < clsNodes[0].propDownstreamCap*0.9)
         clsNodes[0].propEffectiveCap[FALL] *= 0.75;*/
		//if ( !converged )
		//	cout << "[WARNING] Simulation for RC tree driven by node '" << clsNodeNames[0] << "' did not converge within " << maxIterations << " iterations.\n";
	} // end method
    
}; // end class

// =============================================================================
// DEFAULT TREE SIMULATION STRATEGY
// =============================================================================

//typedef RCTreeBase RCTree;
typedef RCTreeReimann RCTree;
//typedef RCTreeFlach RCTree;
//typedef RCTreeDefault RCTree;

#endif
