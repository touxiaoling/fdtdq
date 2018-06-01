#include "stdafx.h"
#include "other.h"
#include "time.h"
using namespace std;
//程序控制常量
#define MAXIMUM_ITERATION 2000
//要计算的总步数
#define PLOT_MODULUS 20
//程序将每PLOT_MODULUS时间步输出3D数据，除了最后一次迭代计算，总是输出。
//因此，如果MAXIMUM_ITERATION不是PLOT_MODULUS的整数倍，则最后一个时间步输出将比隔开先前输出的时间间隔更短。
#define FREQUENCY 10.0e9     
//计算频率（Hz）
#define GUIDE_X 0.025
//（m）
#define GUIDE_Y 0.0125
//（m）
#define GUIDE_Z 0.06
//（m）

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

int main()
{
	/////////////////////////////////////////////////////////////////////////////
	//变量声明
	int i, j, k;
	const int nx=40, ny=20, nz=96;
	//沿着x，y和z轴的单元总数
	const double dx=GUIDE_X/nx, dy=GUIDE_Y/ny, dz=GUIDE_Z/nz;
	//元胞尺寸（m）
	double dt = dx / (2 * LIGHT_SPEED);
	//时间差(s)
	double fs = 1.0 / dt;
	const int ox = 20, oy = 10, oz = 48;
	const int sx = 25, sy = 10, sz = 40;
	//观察点

	int allocatedBytes = 0;
	//一个跟踪分配字节数的计数器
	int iteration = 0;
	//计数器来跟踪已经计算了多少次时步
	double stimulus = 0.0;
	//给定时间步的激励值
	double currentSimulatedTime = 0.0;
	//模拟时间内仿真进行的时间
	const double totalSimulatedTime = MAXIMUM_ITERATION*dt;
	//以秒为单位的时间将被程序模拟
	double omega= 2.0*M_PI*FREQUENCY;
	//角频率
	double lambda= LIGHT_SPEED / FREQUENCY;
	//波长(m)
	double tau = 0.5 / FREQUENCY;

	//场更新方程中使用的物理常数
	double dtmudx, dtepsdx;
	double dtmudy, dtepsdy; 
	double dtmudz, dtepsdz;
	
	//三维数组指针
	double ***ex, ***ey, ***ez;
	double ***hx, ***hy, ***hz;
	
	double simulationMin = FLT_MAX;
	//跟踪整个模拟输出的最小值
	double simulationMax = -FLT_MAX;
	
	FILE *FilePointer;
	//指向txt文件的指针
	time_t startTime, nowTime, lastTime = 0;
	//在场更新方程中使用的常量：
	dtmudx = dt / (MU_0*dx);
	dtepsdx = dt / (EPSILON_0*dx);
	dtmudy = dt / (MU_0*dy);
	dtepsdy = dt / (EPSILON_0*dy);
	dtmudz = dt / (MU_0*dz);
	dtepsdz = dt / (EPSILON_0*dz);

	// 为E面数组分配内存：

	ex=generateZeroData(nx,ny+1,nz+1);
	allocatedBytes += ((nx)*(ny + 1)*(nz + 1) * sizeof(double));

	ey = generateZeroData(nx+1, ny, nz + 1);
	allocatedBytes += ((nx + 1)*(ny)*(nz + 1) * sizeof(double));

	ez = generateZeroData(nx + 1, ny + 1, nz);
	allocatedBytes += ((nx + 1)*(ny + 1)*(nz) * sizeof(double));

	//分配H字段数组：

	hx = generateZeroData(nx -1, ny, nz);
	allocatedBytes += ((nx - 1)*(ny)*(nz) * sizeof(double));

	hy = generateZeroData(nx, ny -1, nz);
	allocatedBytes += ((nx)*(ny - 1)*(nz) * sizeof(double));

	hz = generateZeroData(nx , ny , nz-1);
	allocatedBytes += ((nx)*(ny)*(nz - 1) * sizeof(double));

	/////////////////////////////////////////////////////////////////////////////
	//向标准输出写入一些进度说明

	//打印输出文件的尺寸：
	fprintf(stdout, "bob -cmap chengGbry.cmap -s %dx%dx%d *.bob\n", nx + 1, ny + 1, nz);
	fprintf(stdout, "\n");
	//打印出分配了多少内存：
	fprintf(stdout, "Dynamically allocated %d bytes\n", allocatedBytes);
	fprintf(stdout, "\n");
	//打印出一些模拟参数：
	fprintf(stdout, "Stimulus = %lg Hertz ", FREQUENCY);
	fprintf(stdout, "%dx%dx%d cells\n", nx, ny, nz);
	fprintf(stdout, "dx=%lg, dy=%lg, dz=%lg meters\n", dx, dy, dz);
	fprintf(stdout, "%lg x %lg x %lg meter^3 simulation region\n",GUIDE_X, GUIDE_Y, GUIDE_Z);
	fprintf(stdout, "\n");
	fprintf(stdout, "Time simulated will be %lg seconds, %d timesteps\n",totalSimulatedTime, MAXIMUM_ITERATION);
	fprintf(stdout, "\n");
	cin.get();
	/////////////////////////////////////////////////////////////////////////////
	//打开并开始编写.txt文件
	//如果需要的话，这个文件将会很方便地将参数提供给matlab

	while (fopen_s(&FilePointer, "matlab.txt", "w"))
	{
		fprintf(stderr, "Difficulty opening ToyFDTD1c.viz");
		perror(" ");
	}


	/////////////////////////////////////////////////////////////////////////////
	// main loop:
	time(&startTime);
	int usesec, allusesec, usedsec;
	for (iteration = 0; iteration < MAXIMUM_ITERATION; iteration++)
	{// mainloop

		//循环的第一次，所有数组的数据都是零。 如果有任何东西不是零，这里就有bug。：>
		currentSimulatedTime = dt*(double)iteration;


		fprintf(FilePointer, "%f ", ey[sx][sy][sz]);

		if ((iteration % PLOT_MODULUS) == 0)
		{
			//仿真时间模拟已经进行：
			system("cls");
			time(&nowTime);
			usedsec = nowTime - startTime;
			usesec = (int)(((double)(nowTime - startTime))*(MAXIMUM_ITERATION - iteration) / (iteration + 1));
			allusesec = (int)(((double)(nowTime - startTime))*(MAXIMUM_ITERATION) / (iteration + 1));
			//打印到标准输出迭代号码和当前模拟时间：
			cout << iteration << " / " << MAXIMUM_ITERATION << " " << currentSimulatedTime << "sec " << endl;
			cout << "use:" << usedsec / 3600 << "h " << (usedsec % 3600) / 60 << "m " << usedsec % 60 << "s   ";

			cout << "need:  " << usesec / 3600 << "h " << (usesec % 3600) / 60 << "m " << usesec % 60 << "s   ";
			cout << "all:  " << allusesec / 3600 << "h " << (allusesec % 3600) / 60 << "m " << allusesec % 60 << "s " << endl;


			cout << "speed:" << ((double)nx*ny*nz*iteration) / (nowTime - startTime) / 1.0e6 << " Mceil/s";
			lastTime = nowTime;

		}

		 /////////////////////////////////////////////////////////////////////////
		 //计算激励
		stimulus = 8 * exp(-(iteration*dt-tau*3)*(iteration*dt-tau * 3)/ tau/tau);
		for (j = 0; j < ny ; j++)
		{
			ey[ox][j][oz] += stimulus;
			
		}
		/////////////////////////////////////////////////////////////////////////
		//更新网格的内部：
		
		//更新hx值：
		for (i = 0; i<(nx - 1); i++)
		{
			for (j = 0; j<(ny); j++)
			{
				for (k = 0; k<(nz); k++)
				{
					hx[i][j][k] += (dtmudz*(ey[i + 1][j][k + 1] - ey[i + 1][j][k]) -
						dtmudy*(ez[i + 1][j + 1][k] - ez[i + 1][j][k]));
				}
			}
		}

		//更新hy值：
		for (i = 0; i<(nx); i++)
		{
			for (j = 0; j<(ny - 1); j++)
			{
				for (k = 0; k<(nz); k++)
				{
					hy[i][j][k] += (dtmudx*(ez[i + 1][j + 1][k] - ez[i][j + 1][k]) -
						dtmudz*(ex[i][j + 1][k + 1] - ex[i][j + 1][k]));
				}
			}
		}

		//更新hz值：
		for (i = 0; i<(nx); i++)
		{
			for (j = 0; j<(ny); j++)
			{
				for (k = 0; k<(nz - 1); k++)
				{
					hz[i][j][k] += (dtmudy*(ex[i][j + 1][k + 1] - ex[i][j][k + 1]) -
						dtmudx*(ey[i + 1][j][k + 1] - ey[i][j][k + 1]));
				}
			}
		}

		//更新E场矢量分量。
		//网格面上的值不会在这里更新; 它们由边界条件计算（和刺激计算）处理。

		//更新ex值：
		for (i = 0; i<(nx); i++)
		{
			for (j = 1; j<(ny); j++)
			{
				for (k = 1; k<(nz); k++)
				{
					ex[i][j][k] += (dtepsdy*(hz[i][j][k - 1] - hz[i][j - 1][k - 1]) -
						dtepsdz*(hy[i][j - 1][k] - hy[i][j - 1][k - 1]));
				}
			}
		}

		//更新ey值：
		for (i = 1; i<(nx); i++)
		{
			for (j = 0; j<(ny); j++)
			{
				for (k = 1; k<(nz); k++)
				{
					ey[i][j][k] += (dtepsdz*(hx[i - 1][j][k] - hx[i - 1][j][k - 1]) -
						dtepsdx*(hz[i][j][k - 1] - hz[i - 1][j][k - 1]));
				}
			}
		}

		//更新ez值：
		for (i = 1; i<(nx); i++)
		{
			for (j = 1; j<(ny); j++)
			{
				for (k = 0; k<(nz); k++)
				{
					ez[i][j][k] += (dtepsdx*(hy[i][j - 1][k] - hy[i - 1][j - 1][k]) -
						dtepsdy*(hx[i - 1][j][k] - hx[i - 1][j - 1][k]));
				}
			}
		}
	}//结束主循环

	 /////////////////////////////////////////////////////////////////////
	 //输出部分：
	 //最后一次重复输出例程以写出最后计算的数据。

	 //仿真时间模拟已经进行：
	currentSimulatedTime = dt*(double)iteration;
	//打印到标准输出迭代号码
	//和当前模拟时间：
	fprintf(stdout, "#%d %lgsec", iteration, currentSimulatedTime);

	//将此迭代的网格角点和最大值和最小值写入viz文件：...
	fprintf(FilePointer, "%lg ", ex[sx][sy][sz]);

	/////////////////////////////////////////////////////////////////////////////
	//关闭这个模拟的viz文件：
	fclose(FilePointer);

	//打印出分配了多少内存： 
	cout<< "Dynamically allocated %d bytes\n"<< allocatedBytes<<endl<<endl;
	//打印出一些模拟参数：
	cout << "PLOT_MODULUS = " << PLOT_MODULUS << endl;
	cout << "Stimulus = " << FREQUENCY <<" Hertz"<< endl;
	cout << nx << ny << nz << "cells" << endl;
	cout<<"dx="<< dx<<"dy="<<dy<<" dz="<<dz<<" meters"<<endl;
	cout << "%lg x %lg x %lg meter^3 simulation region\n" <<
		GUIDE_X << GUIDE_Y << GUIDE_Z << endl;
	fprintf(stdout, "\n");
	fprintf(stdout, "Time simulated was %lg seconds, %d timesteps\n",
		totalSimulatedTime, MAXIMUM_ITERATION);
	
}// end main    

