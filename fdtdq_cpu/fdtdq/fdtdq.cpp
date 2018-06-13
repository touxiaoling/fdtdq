#include "stdafx.h"
#include "other.h"
#include <fstream>
#include <string>
using namespace std;


int main()
{
	/////////////////////////////////////////////////////////////////////////////
	//��������
	int i, j, k;
	const int re = RE, we = WE, je = JE;
	const int rb = re+1, wb = we+1, jb = je+1;
	//����x��y��z��ĵ�Ԫ����
	const double dr = EARTH_HD/re , dw = M_PI/ we, dj = 2*M_PI/ je;
	//Ԫ���ߴ磨m��
	const double dt = dr / (2 * LIGHT_SPEED);
	//ʱ���(s)
	const int samplingPoint = LIGHT_SPEED / (5 * FREQUENCY * dr);
	double fs = 1.0 / dt;
	const int sr = 0, sw = we / 2, sj = je / 2;
	const int on = 5;
	const int or [] = {2}, ow[] = {sw / 2, sw / 2, sw / 2, sw / 2, sw / 2}, oj[] = {sj / 2 + 10, sj / 2 + 20, sj / 2 + 30, sj / 2 + 40, sj / 2 + 50};

	//�۲��
	const double totalSimulatedTime = TOTALSIMULATERTIME;
	const int maximumIteration = (int)(totalSimulatedTime / dt)+1;
	int iteration = 0;
	//�������������Ѿ������˶��ٴ�ʱ��
	double stimulus = 0.0;
	//����ʱ�䲽�ļ���ֵ
	double currentSimulatedTime = 0.0;
	//ģ��ʱ���ڷ�����е�ʱ��
	time_t startTime, nowTime,lastTime=0;
	//����Ϊ��λ��ʱ�佫������ģ��

	//�����·�����ʹ�õ�������
	double c1, c2, c3, c4;
	//��ά����ָ��
	double ***er, ***ew, ***ej;
	double ***hr, ***hw, ***hj;
	int ***i1, ***i2;

	FILE *FilePointer, *FilePointew, *FilePointej;
	FILE *FilePointhr, *FilePointhw, *FilePointhj;

	//ָ��txt�ļ���ָ��

	// ΪE����������ڴ棺

	er = generateZeroData(rb, wb, jb);
	ew = generateZeroData(rb, wb, jb);
	ej = generateZeroData(rb, wb, jb);
	//����H�ֶ����飺

	hr = generateZeroData(rb, wb, jb);
	hw = generateZeroData(rb, wb, jb);
	hj = generateZeroData(rb, wb, jb);

	//�ڳ����·�����ʹ�õĳ�����
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
	//���׼���д��һЩ����˵��
	cout << "all step is " << maximumIteration << endl;
	cout<<"start for ENTER"<<endl;
	//cin.get();
	/////////////////////////////////////////////////////////////////////////////
	string filesite = "./nothing";
	string filename = "/matlab";
	string fileform = ".txt";
	ofstream FileStream[on];
	//ָ��txt�ļ���ָ��
	for (int rr = 0; rr < on; rr++)
	{
		FileStream[rr].open(filesite + filename + to_string(rr) + fileform);
	}
	/////////////////////////////////////////////////////////////////////////////
	// main loop:
	time(&startTime);
	int usesec, allusesec, usedsec;
	for (iteration = 0; iteration < maximumIteration; iteration++)
	{// mainloop

		//ѭ���ĵ�һ�Σ�������������ݶ����㡣 ������κζ��������㣬�������bug����>

		//����߽�����
		for (i = 0; i < (re); i++)
		{
			j = 0;
			{
				double sumnum = 0;
				for (k = 1; k < jb; k++)
				{
					sumnum += hj[i][j][k];
				}
				for (k = 0; k < jb; k++)
				{
					er[i][j][k] = c1 * er[i][j][k] + (sin(dw / 2) * dw / (2 * M_PI * (1 - cos(dw / 2)) * (EARTH_RD + (i1[i][j][k] + 0.5) * dr))) * c2 * sumnum;
				}
			}
		}

		for (i = 0; i < (re); i++)
		{
			j = we;
			{
				double sumnum = 0;
				for (k = 1; k < jb; k++)
				{
					sumnum += hj[i][j][k];
				}
				for (k = 0; k < jb; k++)
				{
					er[i][j][k] = c1 * er[i][j][k] - (sin(dw / 2) * dw / (2 * M_PI * (1 - cos(dw / 2)) * (EARTH_RD + (i1[i][j][k] + 0.5) * dr))) * c2 * sumnum;
				}
			}
		}
		//����erֵ��
		for (i = 0; i < (re); i++)
		{
			for (j = 1; j < (we); j++)
			{
				k = 0;
				{
					er[i][j][k] = c1 * er[i][j][k] + c2 / (((i1[i][j][k] + 0.5) * dr + EARTH_RD) * sin(i2[i][j][k] * dw)) * (sin((i2[i][j][k] + 0.5) * dw) * hj[i][j][k] / dw - sin((i2[i][j][k] - 0.5) * dw) * hj[i][j - 1][k] / dw - (hw[i][j][k] - hw[i][j][je-1]) / dj);
				}
			}
		}
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

		//����ewֵ��
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

		//����ejֵ��
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

		//���ڱ߽�����
		for (i = 0; i<rb; i++)
		{
			for (j = 0; j<wb; j ++)
			{
				er[i][j][je] = er[i][j][0];
				ew[i][j][je] = ew[i][j][0];
			}
		}

		/////////////////////////////////////////////////////////////////////////
		//���㼤��
		er[sr][sw][sj] +=(electric(iteration + 1) - electric(iteration))*dt / (EPSILON_0*EPSILON_R) / (sin(sw*dw)*dw*(EARTH_RD + sr*dr)*dj*(EARTH_RD + sr*dr)*dw*dr);
		//cout<< er[sr][sw][sj]<<endl;
		//����hrֵ��
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

		//����hwֵ��
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

		//����hjֵ��
		for (i = 0; i<re; i++)
		{
			for (j = 0; j<we; j++)
			{
				for (k = 0; k<je; k++)
				{

					hj[i][j][k] = c3*hj[i][j][k] - c4 / (((i1[i][j][k] + 0.5)*dr + EARTH_RD))*(((i1[i][j][k] + 1 + EARTH_RD / dr)*ew[i+1][j][k] - (i1[i][j][k] + EARTH_RD / dr)*ew[i][j][k]) - (er[i][j+1][k] - er[i][j][k]) / dw);
				}
			}
		}
		//���ڱ߽�����
		for (i = 0; i<rb; i++)
		{
			for (j = 0; j<wb; j ++)
			{
				hj[i][j][je] = hr[i][j][1];
			}
		}

		time(&nowTime);
		if (nowTime > lastTime)
		{
			currentSimulatedTime = dt * (double)iteration;
			//����ʱ��ģ���Ѿ����У�
			system("cls");
			//time(&nowTime);
			usedsec = nowTime - startTime;
			usesec = (int)(((double)(nowTime - startTime)) * (maximumIteration - iteration) / (iteration + 1));
			allusesec = (int)(((double)(nowTime - startTime)) * (maximumIteration) / (iteration + 1));
			//��ӡ����׼�����������͵�ǰģ��ʱ�䣺
			cout << iteration << " / " << maximumIteration << " " << currentSimulatedTime << "sec " << endl;
			cout << "use:" << usedsec / 3600 << "h " << (usedsec % 3600) / 60 << "m " << usedsec % 60 << "s   ";

			cout << "need:  " << usesec / 3600 << "h " << (usesec % 3600) / 60 << "m " << usesec % 60 << "s   ";
			cout << "all:  " << allusesec / 3600 << "h " << (allusesec % 3600) / 60 << "m " << allusesec % 60 << "s " << endl;

			cout << "speed:" << ((double)rb * wb * jb * iteration) / (nowTime - startTime) / 1.0e6 << " Mceil/s" << endl;
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

		if (!(iteration % samplingPoint))
		{
			for (int oo = 0; oo < on; oo++)
			{
				for (int r = 0; r < re; r++)
				{
					FileStream[oo] << er[r][ow[oo]][oj[oo]] << " " << ew[r][ow[oo]][oj[oo]] << " " << ej[r][ow[oo]][oj[oo]] << " " << hr[r][ow[oo]][oj[oo]] << " " << hw[r][ow[oo]][oj[oo]] << " " << hj[r][ow[oo]][oj[oo]] << " ";
				}
				FileStream[oo] << endl;
			}
		}
	}//������ѭ��


	/////////////////////////////////////////////////////////////////////////////
	//�ر����ģ���viz�ļ���
	for (int rr = 0; rr < on; rr++)
		FileStream[rr].close();
	cout << "all ok" << endl;
	//	cin.get();
}// end main

