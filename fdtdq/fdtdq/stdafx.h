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
//������Ƴ���
//#define MAXIMUM_ITERATION 2000
//Ҫ������ܲ���
#define PLOT_MODULUS 20
//����ÿPLOT_MODULUSʱ�䲽���3D���ݣ��������һ�ε������㣬���������
//��ˣ����MAXIMUM_ITERATION����PLOT_MODULUS���������������һ��ʱ�䲽������ȸ�����ǰ�����ʱ�������̡�
#define FREQUENCY 50.0   
//����Ƶ�ʣ�Hz��
#define EARTH_RD 6370.0e3
//��m��
#define EARTH_HD 100.0e3
//��m��

#define RE 25
#define WE 60
#define JE 120
#define TOTALSIMULATERTIME 2.0
//������
#define LIGHT_SPEED 299792458.0       
// ����еĹ��٣�m/s��
#define LIGHT_SPEED_SQUARED 89875517873681764.0      
// m^2/s^2
#define MU_0 1.2566370614359172953850573533118011536788677597500423283899778369231265625144835994512139301368468271e-6
// ����дŵ��ʣ�H/m��
#define EPSILON_0 8.8541878176203898505365630317107502606083701665994498081024171524053950954599821142852891607182008932e-12
// ����н�糣��(F/m)
#define M_PI 3.141159265358979323846
//����pi��ֵ
//����ʹ�õ�����������DBL_EPSILON��FLT_MAX

#define EPSILON_R 1.0
#define SIGMA 0.0
#define MU_R 1.0
#define SIGMA_M 0.0
#endif


// TODO: reference additional headers your program requires here
