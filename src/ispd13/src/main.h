/*
 *  main_cl.h
 *  gr
 *
 *  Created by Tiago Reimann on 27/06/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <time.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <iomanip>
#include <float.h>
#include <limits>
#include <cassert>
#include <algorithm>

using namespace std;

void flowTiago(string benchmarkName, string dirRoot);
void flowGraci(string benchmarkName, string dirRoot);
void flowFlach(string benchmarkName, string dirRoot);
void flowJohann(string benchmarkName, string dirRoot);

void flowDefault(string benchmarkName, string dirRoot);
void flowDefaultFast(string benchmarkName, string dirRoot);
void flowDefaultLoad(string benchmarkName, string dirRoot);

