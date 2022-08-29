#ifndef G_PLOTTING_H
#define G_PLOTTING_H

#include "PlottingEngine.hh"

namespace Garfield {

extern PlottingEngine plottingEngine;

inline void SetDefaultStyle() { 
  plottingEngine.SetDefaultStyle();
}

}

#endif
