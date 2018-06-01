#include "stdafx.h"
#include "other.h"
#include "time.h"
using namespace std;
//������Ƴ���
#define MAXIMUM_ITERATION 2000
//Ҫ������ܲ���
#define PLOT_MODULUS 20
//����ÿPLOT_MODULUSʱ�䲽���3D���ݣ��������һ�ε������㣬���������
//��ˣ����MAXIMUM_ITERATION����PLOT_MODULUS���������������һ��ʱ�䲽������ȸ�����ǰ�����ʱ�������̡�
#define FREQUENCY 10.0e9     
//����Ƶ�ʣ�Hz��
#define GUIDE_X 0.025
//��m��
#define GUIDE_Y 0.0125
//��m��
#define GUIDE_Z 0.06
//��m��

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

int main()
{
	/////////////////////////////////////////////////////////////////////////////
	//��������
	int i, j, k;
	const int nx=40, ny=20, nz=96;
	//����x��y��z��ĵ�Ԫ����
	const double dx=GUIDE_X/nx, dy=GUIDE_Y/ny, dz=GUIDE_Z/nz;
	//Ԫ���ߴ磨m��
	double dt = dx / (2 * LIGHT_SPEED);
	//ʱ���(s)
	double fs = 1.0 / dt;
	const int ox = 20, oy = 10, oz = 48;
	const int sx = 25, sy = 10, sz = 40;
	//�۲��

	int allocatedBytes = 0;
	//һ�����ٷ����ֽ����ļ�����
	int iteration = 0;
	//�������������Ѿ������˶��ٴ�ʱ��
	double stimulus = 0.0;
	//����ʱ�䲽�ļ���ֵ
	double currentSimulatedTime = 0.0;
	//ģ��ʱ���ڷ�����е�ʱ��
	const double totalSimulatedTime = MAXIMUM_ITERATION*dt;
	//����Ϊ��λ��ʱ�佫������ģ��
	double omega= 2.0*M_PI*FREQUENCY;
	//��Ƶ��
	double lambda= LIGHT_SPEED / FREQUENCY;
	//����(m)
	double tau = 0.5 / FREQUENCY;

	//�����·�����ʹ�õ�������
	double dtmudx, dtepsdx;
	double dtmudy, dtepsdy; 
	double dtmudz, dtepsdz;
	
	//��ά����ָ��
	double ***ex, ***ey, ***ez;
	double ***hx, ***hy, ***hz;
	
	double simulationMin = FLT_MAX;
	//��������ģ���������Сֵ
	double simulationMax = -FLT_MAX;
	
	FILE *FilePointer;
	//ָ��txt�ļ���ָ��
	time_t startTime, nowTime, lastTime = 0;
	//�ڳ����·�����ʹ�õĳ�����
	dtmudx = dt / (MU_0*dx);
	dtepsdx = dt / (EPSILON_0*dx);
	dtmudy = dt / (MU_0*dy);
	dtepsdy = dt / (EPSILON_0*dy);
	dtmudz = dt / (MU_0*dz);
	dtepsdz = dt / (EPSILON_0*dz);

	// ΪE����������ڴ棺

	ex=generateZeroData(nx,ny+1,nz+1);
	allocatedBytes += ((nx)*(ny + 1)*(nz + 1) * sizeof(double));

	ey = generateZeroData(nx+1, ny, nz + 1);
	allocatedBytes += ((nx + 1)*(ny)*(nz + 1) * sizeof(double));

	ez = generateZeroData(nx + 1, ny + 1, nz);
	allocatedBytes += ((nx + 1)*(ny + 1)*(nz) * sizeof(double));

	//����H�ֶ����飺

	hx = generateZeroData(nx -1, ny, nz);
	allocatedBytes += ((nx - 1)*(ny)*(nz) * sizeof(double));

	hy = generateZeroData(nx, ny -1, nz);
	allocatedBytes += ((nx)*(ny - 1)*(nz) * sizeof(double));

	hz = generateZeroData(nx , ny , nz-1);
	allocatedBytes += ((nx)*(ny)*(nz - 1) * sizeof(double));

	/////////////////////////////////////////////////////////////////////////////
	//���׼���д��һЩ����˵��

	//��ӡ����ļ��ĳߴ磺
	fprintf(stdout, "bob -cmap chengGbry.cmap -s %dx%dx%d *.bob\n", nx + 1, ny + 1, nz);
	fprintf(stdout, "\n");
	//��ӡ�������˶����ڴ棺
	fprintf(stdout, "Dynamically allocated %d bytes\n", allocatedBytes);
	fprintf(stdout, "\n");
	//��ӡ��һЩģ�������
	fprintf(stdout, "Stimulus = %lg Hertz ", FREQUENCY);
	fprintf(stdout, "%dx%dx%d cells\n", nx, ny, nz);
	fprintf(stdout, "dx=%lg, dy=%lg, dz=%lg meters\n", dx, dy, dz);
	fprintf(stdout, "%lg x %lg x %lg meter^3 simulation region\n",GUIDE_X, GUIDE_Y, GUIDE_Z);
	fprintf(stdout, "\n");
	fprintf(stdout, "Time simulated will be %lg seconds, %d timesteps\n",totalSimulatedTime, MAXIMUM_ITERATION);
	fprintf(stdout, "\n");
	cin.get();
	/////////////////////////////////////////////////////////////////////////////
	//�򿪲���ʼ��д.txt�ļ�
	//�����Ҫ�Ļ�������ļ�����ܷ���ؽ������ṩ��matlab

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

		//ѭ���ĵ�һ�Σ�������������ݶ����㡣 ������κζ��������㣬�������bug����>
		currentSimulatedTime = dt*(double)iteration;


		fprintf(FilePointer, "%f ", ey[sx][sy][sz]);

		if ((iteration % PLOT_MODULUS) == 0)
		{
			//����ʱ��ģ���Ѿ����У�
			system("cls");
			time(&nowTime);
			usedsec = nowTime - startTime;
			usesec = (int)(((double)(nowTime - startTime))*(MAXIMUM_ITERATION - iteration) / (iteration + 1));
			allusesec = (int)(((double)(nowTime - startTime))*(MAXIMUM_ITERATION) / (iteration + 1));
			//��ӡ����׼�����������͵�ǰģ��ʱ�䣺
			cout << iteration << " / " << MAXIMUM_ITERATION << " " << currentSimulatedTime << "sec " << endl;
			cout << "use:" << usedsec / 3600 << "h " << (usedsec % 3600) / 60 << "m " << usedsec % 60 << "s   ";

			cout << "need:  " << usesec / 3600 << "h " << (usesec % 3600) / 60 << "m " << usesec % 60 << "s   ";
			cout << "all:  " << allusesec / 3600 << "h " << (allusesec % 3600) / 60 << "m " << allusesec % 60 << "s " << endl;


			cout << "speed:" << ((double)nx*ny*nz*iteration) / (nowTime - startTime) / 1.0e6 << " Mceil/s";
			lastTime = nowTime;

		}

		 /////////////////////////////////////////////////////////////////////////
		 //���㼤��
		stimulus = 8 * exp(-(iteration*dt-tau*3)*(iteration*dt-tau * 3)/ tau/tau);
		for (j = 0; j < ny ; j++)
		{
			ey[ox][j][oz] += stimulus;
			
		}
		/////////////////////////////////////////////////////////////////////////
		//����������ڲ���
		
		//����hxֵ��
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

		//����hyֵ��
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

		//����hzֵ��
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

		//����E��ʸ��������
		//�������ϵ�ֵ�������������; �����ɱ߽��������㣨�ʹ̼����㣩����

		//����exֵ��
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

		//����eyֵ��
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

		//����ezֵ��
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
	}//������ѭ��

	 /////////////////////////////////////////////////////////////////////
	 //������֣�
	 //���һ���ظ����������д������������ݡ�

	 //����ʱ��ģ���Ѿ����У�
	currentSimulatedTime = dt*(double)iteration;
	//��ӡ����׼�����������
	//�͵�ǰģ��ʱ�䣺
	fprintf(stdout, "#%d %lgsec", iteration, currentSimulatedTime);

	//���˵���������ǵ�����ֵ����Сֵд��viz�ļ���...
	fprintf(FilePointer, "%lg ", ex[sx][sy][sz]);

	/////////////////////////////////////////////////////////////////////////////
	//�ر����ģ���viz�ļ���
	fclose(FilePointer);

	//��ӡ�������˶����ڴ棺 
	cout<< "Dynamically allocated %d bytes\n"<< allocatedBytes<<endl<<endl;
	//��ӡ��һЩģ�������
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

