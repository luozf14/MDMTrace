#ifndef MDMFIELDMAPINTEROP_H
#define MDMFIELDMAPINTEROP_H

extern "C" {
void mdmfm_init();
void mdmfm_set_probes(double dipole_probe, double multipole_probe);
void mdmfm_eval_entry_multipole_field(double x, double y, double z,
                                      double* bx, double* by, double* bz);
void mdmfm_eval_dipole_field(double x, double y, double z,
                             double* bx, double* by, double* bz);

extern struct {
  double DATA[200][75];
  double ITITLE[200];
} blck0_;
}

#endif
