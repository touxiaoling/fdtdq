#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "device_functions.h"
#include <stdio.h>
#include <iostream>
#include <tchar.h>
#include <float.h>
#include <math.h>
#include <fstream>
#include <string>
#include <ctime>
#include "sm_20_atomic_functions.h"
using namespace std;
//程序控制常量//需要程序通过变量引入，方便更改。
#define RE 25
#define WE 60
#define JE 120
//元胞划分
#define TX 16
#define TY 6
#define TZ 4
//block划分
#define FREQUENCYRESOLUTION 0.25//频率分辨率
#define TOTALSIMULATERTIME (1.0/FREQUENCYRESOLUTION)//总模拟仿真时间

#define FREQUENCY 50.0
//计算频率（Hz）

#define NOTHING 1231123
//电导率可选参数 KNEE AVG NOTHING

//物理常量//可直接使用。
#define EARTH_RD 6370.0e3
//（m）
#define EARTH_HD 100.0e3
//（m）
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

int divUp(int a, int b)
{
	return (a + b - 1) / b;
}
float electric(int N, float dt)
{
	float rtau = 0.3 / 20;//频率
	float tau = rtau / dt;
	return (1000 * exp(-(N - 3 * tau)*(N - 3 * tau) / tau / tau));
}

__global__ void EKernel1(const float*hj, const int rb, const int wb, const int jb, float * sumnum)
{
	const int x = blockIdx.x*blockDim.x + threadIdx.x;
	const int y = blockIdx.y*blockDim.y + threadIdx.y;
	const int z = blockIdx.z*blockDim.z + threadIdx.z;
	const int i = x + y*rb + z*rb*wb;
	const int i0 = x + y*rb + 0 * rb*wb;

	if ((x < rb - 1) && ((y == 0) || (y == wb - 2)) && (z < jb - 1))
	{
		atomicAdd(&(sumnum[x + (y / (wb - 2))*rb]), hj[i]);
	}
}

__global__ void EKernel2(float *er, float*ew, float*ej, const float*hr, const float*hw, const float*hj, const float* ca, const float* cb, const int rb, const int wb, const int jb, const float dr, const float dw, const float dj, const float *sumnum)
{
	const int x = blockIdx.x*blockDim.x + threadIdx.x;
	const int y = blockIdx.y*blockDim.y + threadIdx.y;
	const int z = blockIdx.z*blockDim.z + threadIdx.z;
	const int i = x + y*rb + z*rb*wb;
	const int ix_1 = (x - 1) + y *rb + z*rb*wb;
	const int iy_1 = x + (y - 1)*rb + z*rb*wb;
	const int iz_1 = x + y*rb + (z - 1)*rb*wb;
	const int izje = x + y*rb + (jb - 2)*rb*wb;

	if ((x < rb - 1) && (y == 0) &&  (z < jb))
	{
		er[i] = ca[x] * er[i] + (sin(dw / 2)*dj / (2 * M_PI*(1 - cos(dw / 2))*(EARTH_RD + (x + 0.5)*dr)))*cb[x] * sumnum[x];
	}
	if ((x < rb - 1) && (y == wb - 1) &&  (z < jb))
	{
		er[i] = ca[x] * er[i] - (sin(dw / 2)*dj / (2 * M_PI*(1 - cos(dw / 2))*(EARTH_RD + (x + 0.5)*dr)))*cb[x] * sumnum[x + rb];
	}

	if ((x < rb - 1) && (y > 0) && (y < wb - 1) && (z == 0))
		er[i] = ca[x] * er[i] + cb[x] / (((x + 0.5)*dr + EARTH_RD)*sin(y * dw))*(sin((y + 0.5)*dw)*hj[i] / dw - sin((y - 0.5)*dw)*hj[iy_1] / dw - (hw[i] - hw[izje]) / dj);
	if ((x < rb - 1) && (y > 0) && (y < wb - 1) && (z > 0) && (z < jb - 1))
		er[i] = ca[x] * er[i] + cb[x] / (((x + 0.5)*dr + EARTH_RD)*sin(y * dw))*(sin((y + 0.5)*dw)*hj[i] / dw - sin((y - 0.5)*dw)*hj[iy_1] / dw - (hw[i] - hw[iz_1]) / dj);

	if ((x > 0) && (x < rb - 1) && (y < wb - 1) && (z == 0))
		ew[i] = ca[x + rb] * ew[i] + cb[x + rb] / ((x*dr + EARTH_RD))*((hr[i] - hr[izje]) / (sin((y + 0.5)*dw)*dj) - (x + 0.5 + EARTH_RD / dr)*hj[i] + (x - 0.5 + EARTH_RD / dr)*hj[ix_1]);
	if ((x > 0) && (x < rb - 1) && (y < wb - 1) && (z > 0) && (z < jb - 1))
		ew[i] = ca[x + rb] * ew[i] + cb[x + rb] / ((x*dr + EARTH_RD))*((hr[i] - hr[iz_1]) / (sin((y + 0.5)*dw)*dj) - (x + 0.5 + EARTH_RD / dr)*hj[i] + (x - 0.5 + EARTH_RD / dr)*hj[ix_1]);

	if ((x > 0) && (x < rb - 1) && (y > 0) && (y < wb - 1) && (z < jb - 1))
		ej[i] = ca[x + 2 * rb] * ej[i] + cb[x + 2 * rb] / (((x)*dr + EARTH_RD))*(-((hr[i] - hr[iy_1]) / (dw)-(x + 0.5 + EARTH_RD / dr)*hw[i] + (x - 0.5 + EARTH_RD / dr)*hw[ix_1]));
	//__syncthreads();




}


__global__ void EKernel3(float *er, float*ew, float*ej, const int sr, const int sw, const int sj, const float S, float *sumnum, const int rb, const int wb, const int jb)
{
	const int x = blockIdx.x*blockDim.x + threadIdx.x;
	const int y = blockIdx.y*blockDim.y + threadIdx.y;
	const int z = blockIdx.z*blockDim.z + threadIdx.z;
	const int i = x + y*rb + z*rb*wb;
	const int ix_1 = (x - 1) + y *rb + z*rb*wb;
	const int iy_1 = x + (y - 1)*rb + z*rb*wb;
	const int iz_1 = x + y*rb + (z - 1)*rb*wb;
	const int iz0 = x + y*rb + 0 * rb*wb;


	if ((x < rb) && (y < wb) && (z == jb - 1))
	{
		er[i] = er[iz0];
		ew[i] = ew[iz0];
	}
	if (x == sr&&y == sw&&z == sj)
		er[i] += S;
	if ((x < rb - 1) && ((y == 0) || (y == wb - 1)) && (z == 0))
	{
		sumnum[x + (y / (wb - 1))*rb] = 0;
	}
}
__global__ void HKernel(const float *er, const float*ew, const float*ej, float*hr, float*hw, float*hj, const float da, const float db, const int rb, const int wb, const int jb, const float dr, const float dw, const float dj)
{
	const int x = blockIdx.x*blockDim.x + threadIdx.x;
	const int y = blockIdx.y*blockDim.y + threadIdx.y;
	const int z = blockIdx.z*blockDim.z + threadIdx.z;
	const int i = x + y*rb + z*rb*wb;
	const int i0 = x + y*rb + 0 * rb*wb;
	const int ix1 = (x + 1) + y *rb + z*rb*wb;
	const int iy1 = x + (y + 1)*rb + z*rb*wb;
	const int iz1 = x + y*rb + (z + 1)*rb*wb;

	if ((x > 0) && (x < rb - 1) && (y < wb - 1) && (z < jb - 1))
		hr[i] = da*hr[i] - db / ((x * dr + EARTH_RD)*sin((y + 0.5)*dw))*((sin((y + 1)*dw)*ej[iy1] - sin(y * dw)*ej[i]) / dw - (ew[iz1] - ew[i]) / dj);
	if ((x < rb - 1) && (y > 0) && (y < wb - 1) && (z < jb - 1))
		hw[i] = da*hw[i] - db / (((x + 0.5) * dr + EARTH_RD))*((er[iz1] - er[i]) / (sin(y * dw)*dj) - ((x + 1 + EARTH_RD / dr)*ej[ix1] - (x + EARTH_RD / dr)*ej[i]));
	if ((x < rb - 1) && (y < wb - 1) && (z < jb))
		hj[i] = da*hj[i] - db / (((x + 0.5)*dr + EARTH_RD))*(((x + 1 + EARTH_RD / dr)*ew[ix1] - (x + EARTH_RD / dr)*ew[i]) - (er[iy1] - er[i]) / dw);

}


int main()
{
	/////////////////////////////////////////////////////////////////////////////
	//变量声明
	const int re = RE, we = WE, je = JE;
	const int rb = re + 1, wb = we + 1, jb = je + 1;
	//沿着x，y和z轴的单元总数
	const double dr = EARTH_HD / re, dw = M_PI / we, dj = 2 * M_PI / je;
	//元胞尺寸（m）
	const double dt = dr / (2 * LIGHT_SPEED);
	double fs = 1.0 / dt;
	const int sr = 0, sw = we / 2, sj = je / 2;
	const int on = 5;
	const int or [] = { 2 }, ow[] = { sw / 2 ,sw / 2,sw / 2,sw / 2,sw / 2 }, oj[] = { sj / 2 + 10,sj / 2 + 20,sj / 2 + 30,sj / 2 + 40,sj / 2 + 50 };
	//观察点
	const double totalSimulatedTime = TOTALSIMULATERTIME;
	const int maximumIteration = (int)(totalSimulatedTime / dt) + 1;
	//const int samplingPoint = maximumIteration/(10*(FREQUENCY / FREQUENCYRESOLUTION));
	const int samplingPoint = LIGHT_SPEED / (5 * FREQUENCY *dr);
	int allocatedBytes = 0;
	//一个跟踪分配字节数的计数器
	int iteration = 0;
	//计数器来跟踪已经计算了多少次时步
	float stimulus = 0.0;
	//给定时间步的激励值
	float currentSimulatedTime = 0.0;
	//以秒为单位的时间将被程序模拟
	time_t startTime, nowTime, lastTime = 0;
	const float tau = (float)(0.5 / FREQUENCY / dt);


	float *ca;// = (2 * EPSILON_0*EPSILON_R - SIGMA*dt) / (2 * EPSILON_0*EPSILON_R + SIGMA*dt);
	float *cb;//= (2 * dt) / (2 * EPSILON_0*EPSILON_R + SIGMA*dt);
	float da = (2 * MU_0*MU_R - SIGMA_M * dt) / (2 * MU_0*MU_R + SIGMA_M*dt);
	float db = 2 * dt / (2 * MU_0*MU_R + SIGMA_M*dt);
	//三维数组指针
	float *er, *ew, *ej;
	float *hr, *hw, *hj;
	float *sumnum;
	float sigma[100];
#ifdef KNEE
	string filesite = "./knee";
#elif AVG
	string filesite = "./avg";
	float sigma_avg[100] = { -13.82032, -13.66546, -13.40335, -13.17205, -12.99382, -12.84368, -12.70918, -12.58078, -12.4639, -12.34632,
		-12.23973, -12.13118, -12.03225, -11.93083, -11.83736, -11.74204, -11.65429, -11.5665,0 - 11.48464, -11.40081,
		-11.32333, -11.24440, -11.17181, -11.09722, -11.02856, -10.95733, -10.89176, -10.82312, -10.75983, -10.69290,
		-10.62880, -10.56880, -10.51152, -10.45234, -10.38773, -10.31917, -10.25097, -10.17755, -10.11117, -10.04006,
		-9.969597, -9.882136, -9.819363, -9.746188, -9.689247, -9.628824, -9.587276, -9.555107, -9.534937, -9.507883,
		-9.481722, -9.464875, -9.441254, -9.408667, -9.377145, -9.292050, -9.220743, -9.104745, -9.013201, -8.864476,
		-8.753869, -8.573321, -8.446165, -8.238668, -8.098838, -7.873779, -7.726281, -7.496615, -7.34718, -7.170000,
		-7.023970, -6.847879, -6.722940, -6.546849, -6.370757, -6.245819, -6.120880, -6.02397, -5.92706, -5.833638,
		-5.756798, -5.662242, -5.584637, -5.486151, -5.405922, -5.287574, -5.194666, -5.050842, -4.941185, -4.766291,
		-4.641053, -4.433463, -4.293188, -4.043739, -3.886141, -3.576347, -3.397373, -3.013649, -2.813157, -2.613600
	};
#else
	string filesite = "./nothing";
#endif
	string filename = "/matlab";
	string fileform = ".txt";
	ofstream FileStream[on];
	//指向txt文件的指针
	for (int rr = 0; rr < on; rr++)
	{
		FileStream[rr].open(filesite + filename + to_string(rr) + fileform);
	}

	const dim3 blockSize(TX, TY, TZ);
	const dim3 gridSize(divUp(rb, TX), divUp(wb, TY), divUp(jb, TZ));

	//udaSetDevice(0);//当有多块GPU时默认使用第0块
	cudaMallocManaged(&er, rb*wb*jb * sizeof(float));
	cudaMallocManaged(&ew, rb*wb*jb * sizeof(float));
	cudaMallocManaged(&ej, rb*wb*jb * sizeof(float));

	cudaMallocManaged(&hr, rb*wb*jb * sizeof(float));
	cudaMallocManaged(&hw, rb*wb*jb * sizeof(float));
	cudaMallocManaged(&hj, rb*wb*jb * sizeof(float));

	cudaMallocManaged(&sumnum, rb * 2 * sizeof(float));
	cudaMallocManaged(&ca, rb * 3 * sizeof(float));
	cudaMallocManaged(&cb, rb * 3 * sizeof(float));

	for (int i = 0; i < rb*wb*jb; i++)
	{
		er[i] = 0;
		ew[i] = 0;
		ej[i] = 0;

		hr[i] = 0;
		hw[i] = 0;
		hj[i] = 0;

	}
	for (int i = 0; i < 2 * rb; i++)
	{
		sumnum[i] = 0;
	}
	//knee
#ifdef KNEE
	for (int h = 0; h < 100; h++)
	{
		if(h<=55)
			sigma[h] = 5.6e-10*exp((h - 55) / 8.3) ;
		else
			igma[h] = 5.6e-10*exp((h - 55) / 2.9);
	}

#elif AVG
	for (int i = 0; i < 100; i++)
	{
		sigma[i] = pow(10, sigma_avg[i]);
	}
#else
	for (int i = 0; i < 100; i++)
	{
		sigma[i] = 0;
	}
#endif
	for (int i = 0; i < re; i++)
	{
		//(2 * EPSILON_0*EPSILON_R - SIGMA*dt) / (2 * EPSILON_0*EPSILON_R + SIGMA*dt)
		ca[i] = (2 * EPSILON_0*EPSILON_R - sigma[2 + 4 * i] * dt) / (2 * EPSILON_0*EPSILON_R + sigma[2 + 4 * i] * dt);

		cb[i] = (2 * dt) / (2 * EPSILON_0*EPSILON_R + sigma[2 + 4 * i] * dt);
	}
	for (int i = 1; i < re; i++)
	{
		ca[i + rb] = (2 * EPSILON_0*EPSILON_R - sigma[4 * i] * dt) / (2 * EPSILON_0*EPSILON_R + sigma[4 * i] * dt);
		ca[i + 2 * rb] = (2 * EPSILON_0*EPSILON_R - sigma[4 * i] * dt) / (2 * EPSILON_0*EPSILON_R + sigma[4 * i] * dt);

		cb[i + rb] = (2 * dt) / (2 * EPSILON_0*EPSILON_R + sigma[4 * i] * dt);
		cb[i + 2 * rb] = (2 * dt) / (2 * EPSILON_0*EPSILON_R + sigma[4 * i] * dt);

	}

	time(&startTime);

	int usesec, allusesec, usedsec;
	for (iteration = 0; iteration < maximumIteration; iteration++)
	{



		EKernel1 << <gridSize, blockSize >> > (hj, rb, wb, jb, sumnum);
		//cudaDeviceSynchronize();

		EKernel2 << <gridSize, blockSize >> > (er, ew, ej, hr, hw, hj, ca, cb, rb, wb, jb, dr, dw, dj, sumnum);

		stimulus = -(electric(iteration + 1, dt) - electric(iteration, dt))*dt / (EPSILON_0*EPSILON_R) / (sin(sw*dw)*dw*(EARTH_RD + sr*dr)*dj*(EARTH_RD + sr*dr)*dw*dr);
		//cudaDeviceSynchronize();

		EKernel3 << <gridSize, blockSize >> > (er, ew, ej, sr, sw, sj, stimulus, sumnum, rb, wb, jb);

		//cudaDeviceSynchronize();

		HKernel << <gridSize, blockSize >> > (er, ew, ej, hr, hw, hj, da, db, rb, wb, jb, dr, dw, dj);

		//cudaDeviceSynchronize();
		//此处不加会在840m上有输出问题，但1080ti上没有
		time(&nowTime);
		if (nowTime>lastTime)
		{
			currentSimulatedTime = dt*(double)iteration;
			//仿真时间模拟已经进行：
			system("cls");
			//time(&nowTime);
			usedsec = nowTime - startTime;
			usesec = (int)(((double)(nowTime - startTime))*(maximumIteration - iteration) / (iteration + 1));
			allusesec = (int)(((double)(nowTime - startTime))*(maximumIteration) / (iteration + 1));
			//打印到标准输出迭代号码和当前模拟时间：
			cout << iteration << " / " << maximumIteration << " " << currentSimulatedTime << "sec " << endl;
			cout << "use:" << usedsec / 3600 << "h " << (usedsec % 3600) / 60 << "m " << usedsec % 60 << "s   ";

			cout << "need:  " << usesec / 3600 << "h " << (usesec % 3600) / 60 << "m " << usesec % 60 << "s   ";
			cout << "all:  " << allusesec / 3600 << "h " << (allusesec % 3600) / 60 << "m " << allusesec % 60 << "s " << endl;


			cout << "speed:" << ((double)rb*wb*jb*iteration) / (nowTime - startTime) / 1.0e6 << " Mceil/s" << endl;
			lastTime = nowTime;
			/*	for (int k = 0; k < jb; k += 4)
			{
			for (int j = 0; j < wb; j += 4)
			{
			cout << er[or +j*rb + k*rb*wb] << " ";
			}
			cout << endl;
			}*/
		}


		if (!(iteration%samplingPoint))
		{
			cudaDeviceSynchronize();
			for (int oo = 0; oo < on; oo++)
			{
				for (int r = 0; r < re; r++)
				{
					FileStream[oo] << er[r + ow[oo] * rb + oj[oo] * rb*wb] << " " <<
						ew[r + ow[oo] * rb + oj[oo] * rb*wb] << " " <<
						ej[r + ow[oo] * rb + oj[oo] * rb*wb] << " " <<
						hr[r + ow[oo] * rb + oj[oo] * rb*wb] << " " <<
						hw[r + ow[oo] * rb + oj[oo] * rb*wb] << " " <<
						hj[r + ow[oo] * rb + oj[oo] * rb*wb] << " ";
				}
				FileStream[oo] << endl;
			}

		}
	}
	for (int rr = 0; rr < on; rr++)
		FileStream[rr].close();
	cout << "all ok" << endl;
	cudaFree(er);
	cudaFree(ew);
	cudaFree(ej);
	cudaFree(hr);
	cudaFree(hw);
	cudaFree(hj);
	cudaFree(sumnum);
	cudaFree(ca);
	cudaFree(cb);

	return 0;
}
