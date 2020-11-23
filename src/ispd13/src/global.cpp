/*
 *  global.cpp
 *  sizer
 *
 *  Created by Tiago Reimann on 11/01/12.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 */

#include "global.h"

LibParserCellInfo* OrgCells::findCellInst(string instType){
	
	for(int i = 0; i<this->oCells.size(); ++i)
		if ( this->oCells[i].footprint.compare(0,this->oCells[i].footprint.size(),instType,0,this->oCells[i].footprint.size())==0 )
			for(int j = 0; j<this->oCells[i].cells.size(); ++j)
				if (instType == this->oCells[i].cells[j].name)
					return &this->oCells[i].cells[j];
	
    return NULL;
}

