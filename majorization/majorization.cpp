/************************************************************************/
/* Author:      weihan                     */
/* Email:       weihan@qq.com                                            */
/* FileName:    PSO.cpp                                                  */
/* LastChange:  2023-2-28                                                 */
/************************************************************************/


#ifndef _PSO_H
#endif
#include <iostream>
#include <istream>
#include <streambuf>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include"PSO.H"
using namespace std;


double temp[1000] = { 0 };//定义的一个中间数组 用来存储读入数据调用操作



//微粒的无参构造函数
PARTICLE::PARTICLE()
{
	X = 0;
	V = 0;
	XBest = 0;
	Dim = 0;
}

//微粒的有参构造函数
PARTICLE::PARTICLE(int n)
{
	Dim = n;
	X = new double[Dim];
	V = new double[Dim];
	XBest = new double[Dim];
}

//微粒的析构函数
PARTICLE::~PARTICLE()
{
	if (Dim)
	{
		delete[] X;
		delete[] V;
		delete[] XBest;
	}
}

//设置微粒的维数
void PARTICLE::SetDim(int d)
{
	if (X)
	{
		delete[] X;
	}
	if (V)
	{
		delete[] V;
	}
	if (XBest)
	{
		delete[] XBest;
	}
	Dim = d;
	X = new double[Dim];
	V = new double[Dim];
	XBest = new double[Dim];
}

//PSO的无参构造函数
PSO::PSO()
{
	Particle = 0;
	PNum = 0;
	GBestIndex = 0;
	Xup = 0;
	Xdown = 0;
	//	W=1;
	W_max = 1;
	W_min = 0.6;

	C1 = 2;
	C2 = 2;
	Com = 0;
}

//PSO的有参构造函数
PSO::PSO(int dim, int n)
{
	Particle = new PARTICLE[n];
	for (int i = 0; i < n; i++)
	{
		Particle[i].SetDim(dim);
	}
	PNum = n;
	GBestIndex = 0;
	Xup = new double[dim];
	Xdown = new double[dim];
	Vmax = new double[dim];
	//	W=1;
	W_max = 1;
	W_min = 0.6;

	C1 = 2;
	C2 = 2;
	Com = 0;
}

//PSO的析构函数
PSO::~PSO()
{
	if (Particle)
	{
		delete[] Particle;
	}
	if (Xup)
	{
		delete[] Xup;
	}
	if (Xdown)
	{
		delete[] Xdown;
	}
	if (Vmax)
	{
		delete[] Vmax;
	}
}
//设置坐标上界
void PSO::SetXup(double* up)
{
	if (!Particle)
	{
		return;
	}
	for (int i = 0; i < Particle[0].Dim; i++)
	{
		Xup[i] = up[i];
	}
}

//设置坐标下界
void PSO::SetXdown(double* d)
{
	if (!Particle)
	{
		return;
	}
	for (int i = 0; i < Particle[0].Dim; i++)
	{
		Xdown[i] = d[i];
	}
}

//设置最大速度
void PSO::SetVmax(double* max)
{
	if (!Particle)
	{
		return;
	}
	for (int i = 0; i < Particle[0].Dim; i++)
	{
		Vmax[i] = max[i];
	}
}

void PSO::SetVmax(double p)
{
	if (!Particle)
	{
		return;
	}
	for (int i = 0; i < Particle[0].Dim; i++)
	{
		Vmax[i] = (Xup[i] - Xdown[i]) * p;
	}
}

//初始化群体
void PSO::Initialize()
{
	if (!Particle)
	{
		return;
	}
	static int kk = (unsigned)time(NULL);
	srand((unsigned)time(NULL) + kk++);

	GBestIndex = 0;

	//初始化所有粒子的个体
	for (int i = 0; i < PNum; i++)
	{
		for (int j = 0; j < Particle[i].Dim; j++)
		{
			Particle[i].X[j] = rand() / (double)RAND_MAX * (Xup[j] - Xdown[j]) + Xdown[j];//随机初始化坐标
			Particle[i].XBest[j] = Particle[i].X[j];
			Particle[i].V[j] = rand() / (double)RAND_MAX * Vmax[j] - Vmax[j] / 2;//随机初始化速度
		}
		Particle[i].Fit = GetFit(Particle[i]);//计算每个微粒适合度
		Particle[i].FitBest = Particle[i].Fit;//计算最优适合度值
		if (Particle[i].Fit > Particle[GBestIndex].Fit)
		{
			//如果这个鸟的适合度大于群体的最大适合度的话，记录下查找群体的最优微粒
			GBestIndex = i;
		}
	}
}
//计算群体各个微粒的适合度

void PSO::CalFit()
{
	if (!Particle)
	{
		return;
	}
	for (int i = 0; i < PNum; i++)
	{
		Particle[i].Fit = GetFit(Particle[i]);
	}
}
//微粒飞翔，产生新一代微粒

void PSO::ParticleFly()
{
	static double FitBak[100];//用来存放备份的适合度值
	if (!Particle)
	{
		return;
	}
	static int tt = (unsigned)time(NULL);
	srand((unsigned)time(NULL) * tt++);

	static int kk = 2;//迭代次数
	double W;
	W = W_max - kk * (W_max - W_min) / IteorMax;
	kk++;

	//整个群体飞向新的位置
	for (int i = 0; i < PNum; i++)
	{
		for (int j = 0; j < Particle[i].Dim; j++)
		{
			Particle[i].V[j] = W * Particle[i].V[j] +
				rand() / (double)RAND_MAX * C1 * (Particle[i].XBest[j] - Particle[i].X[j]) +
				rand() / (double)RAND_MAX * C2 * (Particle[GBestIndex].XBest[j] - Particle[i].X[j]);

		}
		for (int j = 0; j < Particle[i].Dim; j++)
		{
			if (Particle[i].V[j] > Vmax[j])
			{
				Particle[i].V[j] = Vmax[j];
			}
			if (Particle[i].V[j] < -Vmax[j])
			{
				Particle[i].V[j] = -Vmax[j];
			}
		}
		for (int j = 0; j < Particle[i].Dim; j++)
		{
			Particle[i].X[j] += Particle[i].V[j];//修改坐标
			if (Particle[i].X[j] > Xup[j])
			{
				Particle[i].X[j] = Xup[j];
			}
			if (Particle[i].X[j] < Xdown[j])
			{
				Particle[i].X[j] = Xdown[j];
			}
		}
	}

	//计算各微粒的适合度
	CalFit();
	for (int i = 0; i < PNum; i++)
	{
		FitBak[i] = Particle[i].Fit;
	}
	//设置个体的最好位置
	for (int i = 0; i < PNum; i++)
	{
		if (Particle[i].Fit >= Particle[i].FitBest)
		{
			Particle[i].FitBest = Particle[i].Fit;
			for (int j = 0; j < Particle[i].Dim; j++)
			{
				Particle[i].XBest[j] = Particle[i].X[j];
			}
		}
	}

	//设置群体中新的最优个体
	GBestIndex = 0;
	for (int i = 0; i < PNum; i++)
	{
		if ((Particle[i].FitBest >= Particle[GBestIndex].FitBest) && i != GBestIndex)
		{
			GBestIndex = i;
		}
	}
}

//按最多运行次数运行群粒算法，返回最优粒子
PARTICLE& PSO::Run(int n)
{
	Initialize();
	double* opt_p = new double[Particle[0].Dim];//通讯用数组，最优点坐标
	double** opt_a = new double* [PNum];		  //通讯用数组，所有点坐标

	for (int i = 0; i < n; i++)
	{
		ParticleFly();
		if (Com)			//通讯函数存在，完成通讯
		{
			for (int k = 0; k < Particle[0].Dim; k++)
			{
				opt_p[k] = Particle[GBestIndex].XBest[k];//拷贝当前最优点坐标
			}
			for (int k = 0; k < PNum; k++)
			{
				opt_a[k] = Particle[k].X;//指向所有点坐标
			}
			if (!Com(Particle[GBestIndex].FitBest, opt_p, opt_a, GBestIndex))
			{
				break;
			}
		}
	}
	delete[] opt_p;
	delete[] opt_a;
	return Particle[GBestIndex];
}

//按最佳适合度运行群粒算法
PARTICLE& PSO::Run(double fit)
{
	Initialize();
	double* opt_p = new double[Particle[0].Dim];//通讯用数组，最优点坐标
	double** opt_a = new double* [PNum];		  //通讯用数组，所有点坐标

	do
	{
		ParticleFly();
		if (Com)			//通讯函数存在，完成通讯
		{
			for (int k = 0; k < Particle[0].Dim; k++)
			{
				opt_p[k] = Particle[GBestIndex].XBest[k];//拷贝最优点坐标
			}
			for (int k = 0; k < PNum; k++)
			{
				opt_a[k] = Particle[k].X;//指向所有点坐标
			}
			if (!Com(Particle[GBestIndex].FitBest, opt_p, opt_a, GBestIndex))
			{
				break;
			}
		}
	} while (Particle[GBestIndex].FitBest < fit);
	delete[] opt_p;
	delete[] opt_a;
	return Particle[GBestIndex];
}

//返回最佳个体
double PSO::GetBest(double* r)
{
	for (int i = 0; i < Particle[GBestIndex].Dim; i++)
	{
		r[i] = Particle[GBestIndex].XBest[i];
	}
	return Particle[GBestIndex].FitBest;
}





bool read()
{
	int data_num;
	int demension_count = 0;//计数维数
	int count1 = 0;//计数数组光标
	int count2 = 1;//用来计数数据块数
	int count3 = 0;//用来存储维数位的光标
	int demension = 0;//用来传输维数

	////******读数据操作  将数据读入一维数组里

	 cout << "请输入你要读取的数据块数：" <<   endl;
	 cin >> data_num;

	 count3 = count1;
	 count1++;
	 temp[count1] = count2;
	 count1++;
	 for (int iii = 1; iii <= data_num; iii++)
	 {
		 string DataFile_name;
		 cout << "请输入第"<<iii<<"个数据文件名（加后缀）:" << endl;
		 cin >> DataFile_name;
		 ifstream csv_data(DataFile_name, ios::in);//数据文件为test
		 string line;//字符串line用来存储一行的数据

		 if (!csv_data.is_open())//打开成功检验
		 {
			 cout << "Error: opening file fail" << endl;
			 exit(1);
			 return false;
		 }
		 else {
			 cout << "the data"<<iii<<"is opened successfully!!!!" << endl;
		 }

		 istringstream sin;         //将整行字符串line读入到字符串istringstream中
		 vector< string> words; //声明一个字符串向量
		 string word;
		 // 读取标题行
		 // getline(csv_data, line);

		 // 读取数据
		 while (getline(csv_data, line))
		 {

			 sin.clear();
			 sin.str(line);
			 words.clear();

			 demension_count = 0;
			 while (getline(sin, word, ',')) //将字符串流sin中的字符读到field字符串中，以逗号为分隔符
			 {
				 demension_count++;
				 double n = atof(word.c_str());
				 temp[count1] = n;
				 words.push_back(word); //将每一格中的数据逐个push
				 //outFile << word <<  endl;

				 count1++;
			 }
			 temp[count1] = -3.1415926;
			 count1++;
			 temp[count3] = demension_count;
		 }
		 

		 temp[count1] = demension;
		 csv_data.close();//读数据关闭

	 }
	
	return true;
}







int main()//用来读取数据和操作
{
	read();
	PARTICLE A;
	A.SetDim(3);
	for (int i = 0; i < 1000; i++)
	{
		cout << temp[i] << endl;
	}
	

	
	return 0;
	
}





//读写数据
//bool read()
//{
//	int count1 = 0;//计数行
//	int count2 = 0;//计数列
//	////******写数据操作
//	 ofstream outFile;
//	outFile.open("test1.csv",  ios::out |  ios::trunc);//创建文件test1
//	outFile << "函数x值" <<  endl;// 写入标题行
//
//	////******读数据操作
//	 ifstream csv_data("test.csv",  ios::in);//数据文件为test
//	 string line;//字符串line用来存储一行的数据
//
//	if (!csv_data.is_open())//打开成功检验
//	{
//		 cout << "Error: opening file fail" <<  endl;
//		 exit(1);
//		return false;
//	}
//	else {
//		 cout << "opening file successfully!!!!" <<  endl;
//	}
//
//	 istringstream sin;         //将整行字符串line读入到字符串istringstream中
//	 vector< string> words; //声明一个字符串向量
//	 string word;
//	// 读取标题行
//	 getline(csv_data, line);
//	// 读取数据
//	while ( getline(csv_data, line))
//	{
//
//		sin.clear();
//		sin.str(line);
//		words.clear();
//		while ( getline(sin, word, ',')) //将字符串流sin中的字符读到field字符串中，以逗号为分隔符
//		{
//			double n = atof(word.c_str());
//			temp[count1][count2] = n;
//			words.push_back(word); //将每一格中的数据逐个push
//
//
//			//outFile << word <<  endl;
//
//			count2++;
//		}
//		count2 = 0;
//
//		 string temmp;
//		words.push_back(temmp);
//
//		 cout <<  endl;
//		
//		count1++;
//	}
//	csv_data.close();//读数据关闭
//
//	outFile.close();//写数据关闭
//	
//	return true;
//}