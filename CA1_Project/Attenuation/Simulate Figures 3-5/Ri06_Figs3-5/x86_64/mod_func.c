#include <stdio.h>
#include "hocdec.h"
extern int nrnmpi_myid;
extern int nrn_nobanner_;

extern void _dv_reg(void);
extern void _ih_reg(void);
extern void _spines_reg(void);
extern void _synampa_reg(void);

void modl_reg(){
  if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
    fprintf(stderr, "Additional mechanisms from files\n");

    fprintf(stderr," dv.mod");
    fprintf(stderr," ih.mod");
    fprintf(stderr," spines.mod");
    fprintf(stderr," synampa.mod");
    fprintf(stderr, "\n");
  }
  _dv_reg();
  _ih_reg();
  _spines_reg();
  _synampa_reg();
}
