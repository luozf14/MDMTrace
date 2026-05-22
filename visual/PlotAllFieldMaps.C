#include "PlotFieldSlices.C"
#include "PlotFieldProfiles.C"
#include "PlotDipoleTopDown.C"

void PlotAllFieldMaps() {
  PlotFieldSlices("Multipole.bin");
  PlotFieldProfiles("Multipole.bin");

  PlotDipoleTopDown("DipoleEntrance.bin", "DipoleSector.bin",
                    "DipoleExit.bin");

  PlotFieldProfiles("DipoleEntrance.bin");
  PlotFieldProfiles("DipoleSector.bin");
  PlotFieldProfiles("DipoleExit.bin");
}
