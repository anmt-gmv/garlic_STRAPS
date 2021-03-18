#ifndef RESIDUALS_RECORDER_H_
#define RESIDUALS_RECORDER_H_

#include "dbtypedef.h"

void open_residuals_file(const char * filename, int index);
void write_residual(double tow, obsdata_t * residual, int index);
void close_residuals_file(int index);


#endif /* RESIDUALS_RECORDER_H_ */
