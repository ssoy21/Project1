#pragma once
#include "matplotlibcpp.h"
#include <vector>
#include <string>

using namespace std;


class Multifractal
{
public:

	vector<vector<float>> mask_data;  //�洢��������

	vector<string> geo_name;  //�洢Ԫ������
	vector<vector<float>> geo_data;  //�洢Ԫ�غ�������

	vector<string> Ngeo_name;  //�洢��Ԫ��������
	vector <vector <float>> Ngeo_data; //�洢��Ԫ��������

	vector<float> ele_data;  //�洢Ϊһά�����Ԫ�غ���˳������

public:

	float ele_mean;  //��ֵ
	float ele_median;  //��λ��
	float ele_variance;  //����
	float ele_stdDev;  //��׼��
	float ele_skewness;  //ƫ��
	float ele_kurtosis;  //���
	float ele_maxVal;  //���ֵ
	float ele_minVal;  //��Сֵ

public:

	int GroupNum;  //�ָ���ֵ��Ŀ
	int LineNum;  //�ָ��߶���Ŀ

public:

	float w_min;  //������С�뾶
	float w_max;  //�������뾶
	float w_step;  //�����ƶ�����
	vector<vector<float>> alpha_data; //�洢������ָ��

public:

	// �����Ļ
	void clearScreen();
	// �ӳ��ض�ʱ��
	void delay(int milliseconds);
	// �����ļ�·��
	void getInputFromFileManager();
	// ����ͻ�ѡ��
	int getUserChoice();
	// �����û�ѡ���������ʵ��
	void processChoice();
	// ��ȡ����
	bool readTifFile(const string& filePath, const string& filePath_mask);
	// ��Ԫͳ�Ʒ���
	void calculateStatistics();
	// ����-���ģ�ͣ�C-A��
	void CAmodel();
	// ����-����ģ�ͣ�C-D��
	void CDmodel();
	// ����-Ƶ��ģ�ͣ�C-Q)
	void CQmodel();
	// Ƶ��-���ģ�ͣ�S-A��
	void SAmodel();
	// �ֲ������Է���
	void calculateHeterogeneity();
	//�����û�������ж�Ӧ���ݴ���
	void anomaly_detection(int choice);
};