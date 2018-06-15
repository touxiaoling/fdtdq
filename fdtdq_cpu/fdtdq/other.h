#pragma once
#ifndef _OTHER_H_
#define _OTHER_H_

double *** generateZeroData(int nx, int ny, int nz);
int setnsize(int*n, double with, double lambda_25);
double electric(int N, double dt);
#endif
