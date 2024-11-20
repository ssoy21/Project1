#pragma once
#include "matplotlibcpp.h"
#include <vector>
#include <string>

using namespace std;


class Multifractal
{
public:

	vector<vector<float>> mask_data;  //存储掩码数据

	vector<string> geo_name;  //存储元素名称
	vector<vector<float>> geo_data;  //存储元素含量数据

	vector<string> Ngeo_name;  //存储非元素列名称
	vector <vector <float>> Ngeo_data; //存储非元素列数据

	vector<float> ele_data;  //存储为一维数组的元素含量顺序数据

public:

	float ele_mean;  //均值
	float ele_median;  //中位数
	float ele_variance;  //方差
	float ele_stdDev;  //标准差
	float ele_skewness;  //偏度
	float ele_kurtosis;  //峰度
	float ele_maxVal;  //最大值
	float ele_minVal;  //最小值

public:

	int GroupNum;  //分割阈值数目
	int LineNum;  //分割线段数目

public:

	float w_min;  //窗口最小半径
	float w_max;  //窗口最大半径
	float w_step;  //窗口移动步长
	vector<vector<float>> alpha_data; //存储奇异性指数

public:

	// 清空屏幕
	void clearScreen();
	// 延迟特定时间
	void delay(int milliseconds);
	// 输入文件路径
	void getInputFromFileManager();
	// 输入客户选择
	int getUserChoice();
	// 根据用户选择进行流程实现
	void processChoice();
	// 读取数据
	bool readTifFile(const string& filePath, const string& filePath_mask);
	// 多元统计分析
	void calculateStatistics();
	// 含量-面积模型（C-A）
	void CAmodel();
	// 含量-距离模型（C-D）
	void CDmodel();
	// 含量-频率模型（C-Q)
	void CQmodel();
	// 频谱-面积模型（S-A）
	void SAmodel();
	// 局部奇异性分析
	void calculateHeterogeneity();
	//根据用户输入进行对应数据处理
	void anomaly_detection(int choice);
};