#include "stdafx.h"
#include "other.h"

using namespace std;


int main()
{
	/////////////////////////////////////////////////////////////////////////////
	//变量声明
	int i, j, k;
	const int re = RE, we = WE, je = JE;
	const int rb = re+1, wb = we+1, jb = je+1;
	//沿着x，y和z轴的单元总数
	const double dr = EARTH_HD/re , dw = M_PI/ we, dj = 2*M_PI/ je;
	//元胞尺寸（m）
	const double dt = dr / (2 * LIGHT_SPEED);
	//时间差(s)
	double fs = 1.0 / dt;
	const int sr = 1, sw = we/2+1, sj = je/2+1;
	const int or = 1, ow = sw/2+1, oj = sj/2+1-10;
	
	
	//观察点
	const double totalSimulatedTime = TOTALSIMULATERTIME;
	const int maximumIteration = (int)(totalSimulatedTime / dt)+1;
	int iteration = 0;
	//计数器来跟踪已经计算了多少次时步
	double stimulus = 0.0;
	//给定时间步的激励值
	double currentSimulatedTime = 0.0;
	//模拟时间内仿真进行的时间
	time_t startTime, nowTime,lastTime=0;
	//以秒为单位的时间将被程序模拟

	//场更新方程中使用的物理常数
	double c1, c2, c3, c4;
	//三维数组指针
	double ***er, ***ew, ***ej;
	double ***hr, ***hw, ***hj;
	int ***i1, ***i2;

	FILE *FilePointer, *FilePointew, *FilePointej;
	FILE *FilePointhr, *FilePointhw, *FilePointhj;

	//指向txt文件的指针

	// 为E面数组分配内存：

	er = generateZeroData(rb, wb, jb);
	ew = generateZeroData(rb, wb, jb);
	ej = generateZeroData(rb, wb, jb);
	//分配H字段数组：

	hr = generateZeroData(rb, wb, jb);
	hw = generateZeroData(rb, wb, jb);
	hj = generateZeroData(rb, wb, jb);

	//在场更新方程中使用的常量：
	c1 = (2 * EPSILON_0*EPSILON_R - SIGMA*dt) / (2 * EPSILON_0*EPSILON_R + SIGMA*dt);
	c2 = (2 * dt) / (2 * EPSILON_0*EPSILON_R + SIGMA*dt);
	c3 = (2 * MU_0*MU_R - SIGMA_M * dt) / (2 * MU_0*MU_R + SIGMA_M*dt);
	c4 = 2 * dt / (2 * MU_0*MU_R + SIGMA_M*dt);
	i1 = (int ***)malloc((rb) * sizeof(int **));
	for (i = 0; i<(rb); i++)
	{
		i1[i] = (int **)malloc((wb) * sizeof(int *));
		for (j = 0; j<(wb); j++)
		{
			i1[i][j] = (int *)malloc((jb) * sizeof(int));
			for (k = 0; k<(jb); k++)
			{
				i1[i][j][k] = i;
			}
		}
	}
	i2 = (int ***)malloc((rb) * sizeof(int **));
	for (i = 0; i<(rb); i++)
	{
		i2[i] = (int **)malloc((wb) * sizeof(int *));
		for (j = 0; j<(wb); j++)
		{
			i2[i][j] = (int *)malloc((jb) * sizeof(int));
			for (k = 0; k<(jb); k++)
			{
				i2[i][j][k] = j;
			}
		}
	}
	/////////////////////////////////////////////////////////////////////////////
	//向标准输出写入一些进度说明
	cout << "all step is " << maximumIteration << endl;
	cout<<"start for ENTER"<<endl;
	//cin.get();
	/////////////////////////////////////////////////////////////////////////////
	//打开并开始编写.txt文件
	//如果需要的话，这个文件将会很方便地将参数提供给matlab
	while (fopen_s(&FilePointer, "matlaber.txt", "w"))
	{
		fprintf(stderr, "Can't opening matlab.txt");
		perror(" ");
	}
	while (fopen_s(&FilePointew, "matlabew.txt", "w"))
	{
		fprintf(stderr, "Can't opening matlab.txt");
		perror(" ");
	}
	while (fopen_s(&FilePointej, "matlabej.txt", "w"))
	{
		fprintf(stderr, "Can't opening matlab.txt");
		perror(" ");
	}
	while (fopen_s(&FilePointhr, "matlabhr.txt", "w"))
	{
		fprintf(stderr, "Can't opening matlab.txt");
		perror(" ");
	}
	while (fopen_s(&FilePointhw, "matlabhw.txt", "w"))
	{
		fprintf(stderr, "Can't opening matlab.txt");
		perror(" ");
	}
	while (fopen_s(&FilePointhj, "matlabhj.txt", "w"))
	{
		fprintf(stderr, "Can't opening matlab.txt");
		perror(" ");
	}

	/////////////////////////////////////////////////////////////////////////////
	// main loop:
	time(&startTime);
	int lastcellnum=0,cellnum;
	int usesec, allusesec, usedsec;
	for (iteration = 0; iteration < maximumIteration; iteration++)
	{// mainloop

	 //循环的第一次，所有数组的数据都是零。 如果有任何东西不是零，这里就有bug。：>
		currentSimulatedTime = dt*(double)iteration;

		time(&nowTime);
		if (nowTime>lastTime)
		{
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

			
			cout << "speed:" << ((double)rb*wb*jb*iteration) / (nowTime - startTime)/1.0e6<<" Mceil/s";
			lastTime = nowTime;

		}


		//更新er值：
		for (i = 0; i<(re); i++)
		{
			for (j = 1; j<(we); j++)
			{
				for (k = 1; k<(jb); k++)
				{
					er[i][j][k] = c1*er[i][j][k] + c2 / (((i1[i][j][k] + 0.5)*dr + EARTH_RD)*sin(i2[i][j][k] * dw))*(sin((i2[i][j][k] + 0.5)*dw)*hj[i][j][k] / dw - sin((i2[i][j][k] - 0.5)*dw)*hj[i][j - 1][k] / dw - (hw[i][j][k] - hw[i][j][k - 1]) / dj);
				}
			}
		}

		//更新ew值：
		for (i = 1; i<(re); i++)
		{
			for (j = 0; j<(we); j++)
			{
				for (k = 1; k<(jb); k++)
				{
					ew[i][j][k] = c1*ew[i][j][k] + c2 / (((i1[i][j][k])*dr + EARTH_RD))*((hr[i][j][k] - hr[i][j][k-1])/ (sin((i2[i][j][k] + 0.5)*dw)*dj) - (i1[i][j][k] + 0.5 + EARTH_RD / dr)*hj[i][j][k] + (i1[i][j][k] - 0.5 + EARTH_RD / dr)*hj[i-1][j][k]);
				}
			}
		}

		//更新ej值：
		for (i = 1; i<(re); i++)
		{
			for (j = 1; j<(we); j++)
			{
				for (k = 0; k<(je); k++)
				{
					ej[i][j][k] = c1*ej[i][j][k] + c2 / (((i1[i][j][k])*dr +EARTH_RD))*(-((hr[i][j][k] - hr[i][j-1][k]) / (dw)-(i1[i][j][k] + 0.5 + EARTH_RD / dr)*hw[i][j][k] + (i1[i][j][k] - 0.5 + EARTH_RD / dr)*hw[i-1][j][k]));
				}
		}
			}
		//极点边界条件
		for (i = 0; i<(re); i++)
		{
			for (j = 0; j<wb;j+=we)
			{
				double sumnum = 0;
				for (k = 1; k < jb; k++)
				{
					sumnum += hj[i][j][k];
				}
				for (k = 1; k<jb; k++)
				{
					er[i][j][k] = c1*er[i][j][k] + (sin(dw / 2)*dw / (2*M_PI*(1 - cos(dw / 2))*(EARTH_RD + (i1[i][j][k] + 0.5)*dr)))*c2*sumnum;
				}
			}
		}
		//周期边界条件
		for (i = 0; i<rb; i++)
		{
			for (j = 0; j<wb; j ++)
			{
				er[i][j][0] = er[i][j][je];
				ew[i][j][0] = ew[i][j][je];
			}
		}

		/////////////////////////////////////////////////////////////////////////
		//计算激励
		er[sr][sw][sj] +=(electric(iteration + 1) - electric(iteration))*dt / (EPSILON_0*EPSILON_R) / (sin(sw*dw)*dw*(EARTH_RD + sr*dr)*dj*(EARTH_RD + sr*dr)*dw*dr);
		//cout<< er[sr][sw][sj]<<endl;
		//更新hr值：
		for (i = 1; i<re; i++)
		{
			for (j = 0; j<we; j++)
			{
				for (k = 0; k<je; k++)
				{

					hr[i][j][k] = c3*hr[i][j][k] - c4/ ((i1[i][j][k]*dr + EARTH_RD)*sin((i2[i][j][k] + 0.5)*dw))*((sin((i2[i][j][k] + 1)*dw)*ej[i][j+1][k] - sin(i2[i][j][k]*dw)*ej[i][j][k]) / dw - (ew[i][j][k+1] - ew[i][j][k]) / dj);
				}
			}
		}

		//更新hw值：
		for (i = 0; i<re; i++)
		{
			for (j = 1; j<we; j++)
			{
				for (k = 0; k<je; k++)
				{

					hw[i][j][k] = c3*hw[i][j][k] - c4 / (((i1[i][j][k] + 0.5)*dr + EARTH_RD))*((er[i][j][k + 1] - er[i][j][k]) / (sin(i2[i][j][k] * dw)*dj) - ((i1[i][j][k] + 1 + EARTH_RD / dr)*ej[i + 1][j][k] - (i1[i][j][k] + EARTH_RD / dr)*ej[i][j][k]));
				}
			}
		}

		//更新hj值：
		for (i = 0; i<re; i++)
		{
			for (j = 0; j<we; j++)
			{
				for (k = 1; k<je; k++)
				{

					hj[i][j][k] = c3*hj[i][j][k] - c4 / (((i1[i][j][k] + 0.5)*dr + EARTH_RD))*(((i1[i][j][k] + 1 + EARTH_RD / dr)*ew[i+1][j][k] - (i1[i][j][k] + EARTH_RD / dr)*ew[i][j][k]) - (er[i][j+1][k] - er[i][j][k]) / dw);
				}
			}
		}
		//周期边界条件
		for (i = 0; i<rb; i++)
		{
			for (j = 0; j<wb; j ++)
			{
				hr[i][j][je] = hr[i][j][0];
				hw[i][j][je] = hw[i][j][0];
			}
		}
		fprintf(FilePointer, "%E ", er[or ][ow][oj]);
		fprintf(FilePointew, "%E ", ew[or ][ow][oj]);
		fprintf(FilePointej, "%E ", ej[or ][ow][oj]);
		fprintf(FilePointhr, "%E ", hr[or ][ow][oj]);
		fprintf(FilePointhw, "%E ", hw[or ][ow][oj]);
		fprintf(FilePointhj, "%E ", hj[or ][ow][oj]);
	}//结束主循环


	/////////////////////////////////////////////////////////////////////////////
	//关闭这个模拟的viz文件：
	fclose(FilePointer);

	cout << "all simulater ok" << endl;
//	cin.get();
}// end main    

