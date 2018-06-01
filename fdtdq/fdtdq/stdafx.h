// stdafx.h : include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently
//

#pragma once
#ifndef _STDAFX_H_
#define _STDAFX_H_

#include "targetver.h"

#include <stdio.h>
#include <tchar.h>  
#include <float.h>
#include <malloc.h>
#include <iostream>
#include <cmath>//math.h
#include <ctime>
//程序控制常量
//#define MAXIMUM_ITERATION 2000
//要计算的总步数
#define PLOT_MODULUS 20
//程序将每PLOT_MODULUS时间步输出3D数据，除了最后一次迭代计算，总是输出。
//因此，如果MAXIMUM_ITERATION不是PLOT_MODULUS的整数倍，则最后一个时间步输出将比隔开先前输出的时间间隔更短。
#define FREQUENCY 50.0   
//计算频率（Hz）
#define EARTH_RD 6370.0e3
//（m）
#define EARTH_HD 100.0e3
//（m）

#define RE 25
#define WE 60
#define JE 120
#define TOTALSIMULATERTIME 2.0
//物理常数
#define LIGHT_SPEED 299792458.0       
// 真空中的光速（m/s）
#define LIGHT_SPEED_SQUARED 89875517873681764.0      
// m^2/s^2
#define MU_0 1.2566370614359172953850573533118011536788677597500423283899778369231265625144835994512139301368468271e-6
// 真空中磁导率（H/m）
#define EPSILON_0 8.8541878176203898505365630317107502606083701665994498081024171524053950954599821142852891607182008932e-12
// 真空中介电常数(F/m)
#define M_PI 3.141159265358979323846
//用于pi的值
//这里使用的其他常量是DBL_EPSILON和FLT_MAX

#define EPSILON_R 1.0
#define SIGMA 0.0
#define MU_R 1.0
#define SIGMA_M 0.0
#endif


// TODO: reference additional headers your program requires here
