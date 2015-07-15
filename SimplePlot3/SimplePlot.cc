//	SimplePlot.cc - simple plotting functions
//	if TEST is defined, also generates a simple main program to test it.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <complex>
#include <ccomplex>
#include "complex.h"

#include "TCanvas.h"
#include "TGraph.h"
#include "TApplication.h"
#include "TROOT.h"

#include "SimplePlot.hh"






TApplication *SimplePlot::app = 0;
int SimplePlot::nPlot = 0;







void SimplePlot::init()
{
	if(app != 0) return;

	int argc=1;
	char *argv[2]; argv[0] = (char *)"test"; argv[1] = 0;
	app = new TApplication("App",&argc,argv);
}

void SimplePlot::oneD(const char *title, const double y[], int nValues)
{
	if(nValues < 2) {
		fprintf(stderr,"ERROR SimplePlot::oneD: nValues < 2\n");
		return;
	}

	// create and fill x[] with indices into y[]
	double *x = new double[nValues];
	for(int i=0; i<nValues; ++i) x[i] = (double)i;

	graph(title,x,y,nValues);

	delete[] x;
}

void SimplePlot::oneD(const char *title, const float y[], int nValues)
{
	if(nValues < 2) {
		fprintf(stderr,"ERROR SimplePlot::oneD: nValues < 2\n");
		return;
	}

	// create and fill x[] with indices into y[]
	float *x = new float[nValues];
	for(int i=0; i<nValues; ++i) x[i] = (float)i;

	graph(title,x,y,nValues);

	delete[] x;
}

void SimplePlot::graph(const char *title, const double x[], const double y[],
								int nValues)
{
	init();

	// create a canvas with a unique name; it is drawn in its own window
	char name[64];
	sprintf(name,"Plot %d",++nPlot);
	TCanvas *c = new TCanvas(name,name, 400, 400);

	// create the graph (automatically goes into *c); this plots x vs y
	TGraph *g = new TGraph(nValues,x,y);
	g->SetTitle(title);
	g->Draw("APL"); // draw Axes, markers at Points, and Lines
	c->Update();

	// Run the Root event loop, until Quit ROOT is clicked
	fprintf(stderr,"Click on 'File/Quit ROOT' to proceed\n");
	app->Run(true);

	// close the canvas and clean up
	c->Close();
	delete g;
	delete c;
}

void SimplePlot::graph(const char *title, const float x[], const float y[],
								int nValues)
{
	init();

	// create a canvas with a unique name; it is drawn in its own window
	char name[64];
	sprintf(name,"Plot %d",++nPlot);
	TCanvas *c = new TCanvas(name,name, 400, 400);

	// create the graph (automatically goes into *c); this plots x vs y
	TGraph *g = new TGraph(nValues,x,y);
	g->SetTitle(title);
	g->Draw("APL"); // draw Axes, markers at Points, and Lines
	c->Update();

	// Run the Root event loop, until Quit ROOT is clicked
	fprintf(stderr,"Click on 'File/Quit ROOT' to proceed\n");
	app->Run(true);

	// close the canvas and clean up
	c->Close();
	delete g;
	delete c;
}


#ifdef TEST


/**	main() is a simple test program for SimplePlot
 **/
int main(int argc, char *argv[])
{
	printf("oneD(double)\n");
	
   const int NVALUES=200;

    
	double values[NVALUES];
	for(int i=0; i<NVALUES; ++i)
    values[i] = sin(2.0*M_PI*i/50.0);
      
	SimplePlot::oneD("Test sin()",values,NVALUES);
    

	printf("oneD(float)\n");
	float floats[NVALUES];
	for(int i=0; i<NVALUES; ++i)
		floats[i] = cos(2.0*M_PI*i/50.0);
	SimplePlot::oneD("Test cos()",floats,NVALUES);

	printf("graph\n");
	float x[NVALUES], y[NVALUES];
	for(int i=0; i<NVALUES; ++i) {
		x[i] = cos(2.0*M_PI*i/31.0);
		y[i] = cos(2.0*M_PI*i/100.0);
	}
	SimplePlot::graph("Graph",x,y,NVALUES);

	printf("exiting\n");
	return 0;
}

#endif // TEST
