//	SimplePlot.hh - simple plotting functions

#ifndef SIMPLEPLOT_HH
#define SIMPLEPLOT_HH

#include "TApplication.h"

/**	class SimplePlot contains functions to make Root plots.
 *
 *	Each function creates a Root plot in a new window.
 *	Each window blocks until "File/Quit ROOT" is clicked.
 *	The plots scale automatically.
 **/
class SimplePlot {
	static TApplication *app;
	static int nPlot;
	static void init();
public:
	/// oneD() plots a 1-d array.
	/// The plot is displayed in a new window.
	/// These function block until "File/Quit ROOT" is clicked.
	/// The x axis is simply the integers [0,nValues).
	static void oneD(const char *title, const double y[], int nValues);
	static void oneD(const char *title, const float y[], int nValues);

	/// graph() plots a graph of (x[i],y{i]), with lines between adjacent
	/// points. x[] and y[] must be the same type (float or double), and
	/// have the same number of elements.
	/// The plot is displayed in a new window.
	/// These function block until "File/Quit ROOT" is clicked.
	static void graph(const char *title, const double x[], const double y[],
								int nValues);
	static void graph(const char *title, const float x[], const float y[],
								int nValues);
};

#endif // SIMPLEPLOT_HH
