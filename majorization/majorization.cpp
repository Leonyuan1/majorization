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

int main()//用来读取数据和操作
{
	//int row=100;
	//int columns=100;
	//std::cout << "please enter the number of rowsz:" << std::endl;
	//std::cin >> row;
	//std::cout << "please enter the number of columns:" << std::endl;
	//std::cin >> columns;
	double temp[20][10];
	
	////******写数据操作
	std::ofstream outFile;
	outFile.open("test1.csv", std::ios::out | std::ios::trunc);//创建文件test1
	outFile << "函数x值" << std::endl;// 写入标题行

	////******读数据操作
	std::ifstream csv_data("test.csv", std::ios::in);//数据文件为test
	std::string line;//字符串line用来存储一行的数据

	if (!csv_data.is_open())//打开成功检验
	{
		std::cout << "Error: opening file fail" << std::endl;
		std::exit(1);
	}
	else {
		std::cout << "opening file successfully!!!!" << std::endl;
	}

	std::istringstream sin;         //将整行字符串line读入到字符串istringstream中
	std::vector<std::string> words; //声明一个字符串向量
	std::string word;

	// 读取标题行
	std::getline(csv_data, line);
	// 读取数据
	while (std::getline(csv_data, line))
	{
		sin.clear();
		sin.str(line);
		words.clear();
		while (std::getline(sin, word, ',')) //将字符串流sin中的字符读到field字符串中，以逗号为分隔符
		{
			words.push_back(word); //将每一格中的数据逐个push
			//调用粒子群算法进行优化 然后将得到的结果写进test1.csv文件中 
			//模拟测试 构造一个二次函数 y=-x^2+3,给定x一个数列 求解其在数列中的最大值
			//word = word + ",";
			std::cout << word;//用于测试读入的数据是否正确，其中第一行为标题
			outFile << word << std::endl;
		}
		
			std:: string temmp;
			words.push_back(temmp);
		
		std::cout << std::endl;
		// do something
	}

	csv_data.close();//读数据关闭
	outFile.close();//写数据关闭
	{
		std::cout << "writing file successfully!!!" << std::endl;
	}

	//std::vector<std::vector<int>> user_arr;
	//std::ifstream fp("test.csv"); //定义声明一个ifstream对象，指定文件路径
	//std::string line;
	//std::getline(fp, line); //跳过列名，第一行不做处理
	//while (getline(fp, line)) { //循环读取每行数据
	//	std::vector<int> data_line;
	//	std::string number;
	//	std::istringstream readstr(line); //string数据流化
	//	//将一行数据按'，'分割
	//	for (int j = 0; j < 11; j++) { //可根据数据的实际情况取循环获取
	//		getline(readstr, number, ','); //循环读取数据
	//		data_line.push_back(atoi(number.c_str())); //字符串传int
	//	}
	//	user_arr.push_back(data_line); //插入到vector中
	//}
	return 0;
	
}