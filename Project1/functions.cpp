#include <algorithm>
#include <gdal.h>
#include <gdal_priv.h>
#include <Eigen/Dense>
#include <fftw3.h>
#include <limits>
#include <fstream>
#include <ctime>
#include <iostream>
#include <cmath>
#include <numeric>
#include <fstream>
#include <cstdlib>
#include <complex>
#include "functions.h"


using namespace Eigen;
using namespace std;
using namespace matplotlibcpp;

const double _E = 2.71828;

// �����Ļ
void Multifractal::clearScreen() {
    cout << "\033[2J\033[1;1H"; // �����Ļ��������֧�� ANSI ת�����е��ն�
}

// �ӳ�ʱ��
void Multifractal::delay(int milliseconds) {
    clock_t end_time = clock() + milliseconds * CLOCKS_PER_SEC / 1000; // �������ʱ��
    while (clock() < end_time) {
        // �������ѡ��ִ������������߱��ֿ�ѭ��
    }
}

// ģ���ļ���������������
void Multifractal::getInputFromFileManager() {
    string filePath = "D:\\test_data\\cu_ok.tif";
    string filePath_mask = "D:\\test_data\\cu_ok_raster1.tif";
    while (true) { // ѭ��ֱ���ɹ������ļ�·��
        cout << "�����������ļ�·����" << filePath << endl;
        //getline(cin, filePath); // ��ȡ�û�������ļ�·��

        ifstream file(filePath); // ���Դ��ļ�

        // ����ļ��Ƿ�ɹ���
        if (!file) {
            cout << "�޷��������ļ��������ļ�·���Ƿ���ȷ!!" << endl;
            delay(1000);
            clearScreen(); // ��ս���
            // �������ѡ������ӳ٣������Ҫ
            continue; // ����ѭ����������ʾ�����ļ�·��
        }
        else {
            cout << "�����ļ���ȡ�ɹ�!!" << endl;
        }

        file.close(); // �ر��ļ�

        cout << "�����������ļ�·����" << filePath_mask << endl;
        //getline(cin, filePath); // ��ȡ�û�������ļ�·��

        ifstream file_mask(filePath_mask); // ���Դ��ļ�

        // ����ļ��Ƿ�ɹ���
        if (!file_mask) {
            cout << "�޷��������ļ��������ļ�·���Ƿ���ȷ!!" << endl;
            delay(1000);
            clearScreen(); // ��ս���
            // �������ѡ������ӳ٣������Ҫ
            continue; // ����ѭ����������ʾ�����ļ�·��
        }
        else {
            cout << "�����ļ���ȡ�ɹ�!!" << endl;
        }

        file_mask.close(); // �ر��ļ�

        //��ȡ�ļ��Ĳ���
        readTifFile(filePath, filePath_mask);
        break; // �ɹ���ȡ�ļ�������ѭ��
    }
}

// ��������ȡ .tif �ļ�����
bool Multifractal::readTifFile(const string& filePath, const string& filePath_mask) {
    GDALAllRegister(); // ע�����и�ʽ

    // ���ļ�
    GDALDataset* dataset = (GDALDataset*)GDALOpen(filePath.c_str(), GA_ReadOnly);
    GDALDataset* dataset_mask = (GDALDataset*)GDALOpen(filePath_mask.c_str(), GA_ReadOnly);

    if (!dataset) {
        cerr << "�޷��������ļ���" << filePath << endl;
        return false;
    }
    if (!dataset_mask) {
        cerr << "�޷��������ļ���" << filePath_mask << endl;
        return false;
    }

    int xSize = dataset->GetRasterXSize(); // ��ȡ�������С
    int ySize = dataset->GetRasterYSize(); // ��ȡ�������С
    float* rasterData = new float[xSize * ySize]; // ��ʱ����

    // ���� geo_data �Ĵ�С
    geo_data.resize(ySize, vector<float>(xSize));

    // ��ȡ����
    GDALRasterBand* band = dataset->GetRasterBand(1); // ��ȡ��һ������
    CPLErr err = band->RasterIO(GF_Read, 0, 0, xSize, ySize, rasterData, xSize, ySize, GDT_Float32, 0, 0);


    if (err != CE_None) {
        cerr << "��ȡ����ʧ�ܡ�" << endl;
        delete[] rasterData; // �ͷ���ʱ����
        GDALClose(dataset); // �ر��ļ�
        return false;
    }

    // �����ݴ洢�� geo_data�����޳���ֵ
    for (int y = 0; y < ySize; ++y) {
        for (int x = 0; x < xSize; ++x) {
            float value = rasterData[y * xSize + x];
            geo_data[y][x] = value; // �洢��Ч����

        }
    }

    // ���� geo_data �Ĵ�С
    mask_data.resize(ySize, vector<float>(xSize));

    // ��ȡ����

    band = dataset_mask->GetRasterBand(1); // ��ȡ��һ������
    err = band->RasterIO(GF_Read, 0, 0, xSize, ySize, rasterData, xSize, ySize, GDT_Float32, 0, 0);


    if (err != CE_None) {
        cerr << "��ȡ����ʧ�ܡ�" << endl;
        delete[] rasterData; // �ͷ���ʱ����
        GDALClose(dataset_mask); // �ر��ļ�
        return false;
    }
    // �����ݴ洢�� geo_data�����޳���ֵ
    for (int y = 0; y < ySize; ++y) {
        for (int x = 0; x < xSize; ++x) {
            float value = rasterData[y * xSize + x];
            mask_data[y][x] = value; // �洢��Ч����
        }
    }

    int nan_data = 0;
    //�޳���Чֵ
    for (int y = 0; y < ySize; ++y) {
        for (int x = 0; x < xSize; ++x) {
            if (mask_data[y][x] == 0) {
                ele_data.push_back(geo_data[y][x]);
            }
            else
            {
                nan_data++;
                geo_data[y][x] = NAN;
            }
        }
    }
    //cout << ele_data[0];

    //�Ժ������ݽ��д�С��������
    sort(ele_data.begin(), ele_data.end());
    cout << "��ֵ��" << nan_data << "��!!" << endl;

    delete[] rasterData; // �ͷ���ʱ����

    GDALClose(dataset); // �ر��ļ�
    GDALClose(dataset_mask);

    return true;
}

// ����Ƿ��� NaN
bool is_nan(float val) {
    return isnan(val);
}

// ��ֵ�����������ƽ��ֵ
void interpolate_once(vector<vector<float>>& data) {
    int rows = data.size();
    int cols = data[0].size();

    // ����һ���������洢��ֵ���
    vector<vector<float>> copy = data;

    // ������������
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            if (is_nan(data[i][j])) {
                float sum = 0.0f;
                int count = 0;

                // ���8����
                for (int di = -1; di <= 1; ++di) {
                    for (int dj = -1; dj <= 1; ++dj) {
                        int ni = i + di, nj = j + dj;
                        if (ni >= 0 && ni < rows && nj >= 0 && nj < cols && !is_nan(data[ni][nj])) {
                            sum += data[ni][nj];
                            count++;
                        }
                    }
                }

                // �������������Чֵ������ƽ��ֵ�
                if (count > 0) {
                    copy[i][j] = sum / count;
                }
            }
        }
    }

    // ����ԭ����
    data = copy;
}

// ������ֵ����
void iterative_interpolation(vector<vector<float>>& data, int max_iter = 100, float tolerance = 1e-6) {
    int rows = data.size();
    int cols = data[0].size();

    for (int iter = 0; iter < max_iter; ++iter) {
        vector<vector<float>> previous = data;

        // ִ��һ�β�ֵ
        interpolate_once(data);

        // ����Ƿ�����
        float max_diff = 0.0f;
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                max_diff = max(max_diff, abs(data[i][j] - previous[i][j]));
            }
        }

        // ����仯С�����̶ȣ���Ϊ��ֵ�Ѿ�����
        if (max_diff < tolerance) {
            cout << "��ֵ����������������" << iter + 1 << endl;
            break;
        }
    }
}

// ��Ԫͳ�Ʒ���
void Multifractal::calculateStatistics() {
    // ��ֵ
    ele_mean = accumulate(ele_data.begin(), ele_data.end(), 0.0) / ele_data.size();
    cout << "��ֵ��" << ele_mean << endl;
    delay(1000);

    // ��λ��
    vector<float> sortedData = ele_data;
    sort(sortedData.begin(), sortedData.end());
    size_t n = sortedData.size();
    ele_median = (n % 2 == 0) ? (sortedData[n / 2 - 1] + sortedData[n / 2]) / 2.0 : sortedData[n / 2];
    cout << "��λ����" << ele_median << endl;
    delay(1000);

    // ����
    ele_variance = 0.0;
    for (float x : ele_data) {
        ele_variance += (x - ele_mean) * (x - ele_mean);
    }
    ele_variance /= ele_data.size();
    cout << "���" << ele_variance << endl;
    delay(1000);

    // ��׼��
    ele_stdDev = sqrt(ele_variance);
    cout << "��׼�" << ele_stdDev << endl;
    delay(1000);

    // ƫ��
    ele_skewness = 0.0;
    for (float x : ele_data) {
        ele_skewness += pow((x - ele_mean) / ele_stdDev, 3);
    }
    ele_skewness /= ele_data.size();
    cout << "ƫ�ȣ�" << ele_skewness << endl;
    delay(1000);

    // ��ȣ�����ȣ�
    ele_kurtosis = 0.0;
    for (float x : ele_data) {
        ele_kurtosis += pow((x - ele_mean) / ele_stdDev, 4);
    }
    ele_kurtosis = ele_kurtosis / ele_data.size() - 3.0;
    cout << "��ȣ�" << ele_kurtosis << endl;
    delay(1000);

    // ���ֵ����Сֵ
    ele_maxVal = *max_element(ele_data.begin(), ele_data.end());
    ele_minVal = *min_element(ele_data.begin(), ele_data.end());
    cout << "���ֵ��" << ele_maxVal << endl;
    cout << "��Сֵ��" << ele_minVal << endl;
    delay(1000);
}

// �Ե���ѧ���ݽ�����Ӧ�����̼���
void Multifractal::processChoice() {
    cout << "���ݶ�Ԫͳ�Ʒ�����" << endl;
    cout << "================================" << endl;
    calculateStatistics();
    cout << "================================" << endl;

    cout << "";


}

// ������Ϻ���������б�ʺͽؾ�
pair<float, float> linearFit(const vector<float>& x, const vector<float>& y) {
    int n = x.size();
    MatrixXd A(n, 2);
    VectorXd Y(n);
    for (int i = 0; i < n; ++i) {
        A(i, 0) = x[i];
        A(i, 1) = 1.0;  // �ؾ�
        Y(i) = y[i];
    }

    // ʹ��QR�ֽ���� Ax = Y
    Vector2d coeffs = A.colPivHouseholderQr().solve(Y);
    float slope = coeffs(0);
    float intercept = coeffs(1);

    return { slope, intercept };
}

// ��С���˷��ֶ���Ϻ���
vector<float> segmentFitting(const vector<float>& refs, const vector<float>& nc, int k) {
    int n = refs.size();
    int segment_len = n / k;

    vector<float>log_refs, log_nc;
    for (size_t i = 0; i < refs.size(); ++i) {
        log_refs.push_back(log(refs[i]));
        log_nc.push_back(log(nc[i]));
    }

    scatter(log_refs, log_nc);
    xlabel("log(refs)");
    ylabel("log(nc)");
    title("log-log");
    vector<string> colors = { "r", "g", "b", "m", "c", "y", "k" }; // ���Ը�����Ҫ���Ӹ�����ɫ

    vector<float> slopes(k), intercepts(k);
    vector<float> intersection_x(k - 1);
    vector<pair<float, float>> intersections;  // �洢���� (x, y)

    // �ֶ���С�������
    for (int i = 0; i < k; ++i) {
        vector< float> x_segment, y_segment;
        if (i == k - 1) {
            x_segment = vector< float>(log_refs.begin() + i * segment_len, log_refs.end());
            y_segment = vector< float>(log_nc.begin() + i * segment_len, log_nc.end());
        }
        else {
            x_segment = vector< float>(log_refs.begin() + i * segment_len, log_refs.begin() + (i + 1) * segment_len);
            y_segment = vector< float>(log_nc.begin() + i * segment_len, log_nc.begin() + (i + 1) * segment_len);
        }

        // ��ȡÿ�ε�б�ʺͽؾ�
        auto [slope, intercept] = linearFit(x_segment, y_segment);
        slopes[i] = slope;
        intercepts[i] = intercept;

        // �����߶ε������յ�
        float x_start = x_segment.front();
        float x_end = x_segment.back();
        float y_start = slope * x_start + intercept;
        float y_end = slope * x_end + intercept;

        // ����ÿ�������
        plot({ x_start, x_end }, { y_start, y_end }, colors[i % colors.size()] + "-"); // ʹ�ò�ͬ��ɫ�������
        
        // ����������һ�Σ��������ڶε��ص���
        if (i < k - 1) {
            float x_overlap = x_segment.back();  // ��ǰ�ε�ĩβ��
            float y1 = slopes[i] * x_overlap + intercepts[i];       // ��ǰ�ε� y ֵ
            float y2 = slopes[i + 1] * x_overlap + intercepts[i + 1];  // ��һ�ε� y ֵ

            intersections.push_back({ x_overlap, y1 });  // �����¼ (x_overlap, y2)
        }
    }


    // ���ƽ���
    vector<float> section_x, section_y;
    // ������������
    cout << "Intersection X-coordinates:\n";
    for (const auto& intersection : intersections) {
        section_x.push_back(intersection.first);
        section_y.push_back(intersection.second);
        cout << intersection.first << "  ";

    }
    scatter(section_x, section_y, 10.0, { {"color", "red"} });  // ��������
    show();

    return section_x;
}

// ѭ�����
vector<float> residual(const vector<float>& refs, const vector<float>& nc) {
    // x,y����
    vector<float> x, y;
    vector<string> colors = { "r", "g", "b", "m", "c", "y", "k" }; // ���Ը�����Ҫ���Ӹ�����ɫ

    // x,y�����ת��
    for (int i = 0; i < refs.size(); i++) {
        x.push_back(log(refs[i]));
        y.push_back(log(nc[i]));
    }

    // ������ϲ���
    float residual_threshold = 0.2;  // ��ϵ�����в�
    int min_points_per_segment = 6;   // ÿ����ϵ���С����

    vector<tuple<float, float, float, float>> segments;  // �洢ÿ�ε���ϲ���
    vector<float> current_x = { x[0] };  // ��ǰ�߶ε� x ����
    vector<float> current_y = { y[0] };  // ��ǰ�߶ε� y ����

    // �����������ݵ�������
    for (size_t i = 1; i < x.size(); ++i) {
        current_x.push_back(x[i]);
        current_y.push_back(y[i]);

        // ʹ�� Eigen �������Իع����
        VectorXd X(current_x.size());
        VectorXd Y(current_y.size());
        for (size_t j = 0; j < current_x.size(); ++j) {
            X(j) = current_x[j];
            Y(j) = current_y[j];
        }

        // ���Իع飺y = mx + b��ʹ����С���˷����
        MatrixXd A(current_x.size(), 2);
        for (size_t j = 0; j < current_x.size(); ++j) {
            A(j, 0) = current_x[j];  // x
            A(j, 1) = 1;              // 1 (������)
        }
        VectorXd coeff = A.colPivHouseholderQr().solve(Y);

        float slope = coeff(0);  // б��
        float intercept = coeff(1);  // �ؾ�

        // ���㵱ǰ���Ԥ��ֵ
        float predicted_y = slope * x[i] + intercept;

        // ����в�
        float residual = abs(y[i] - predicted_y);

        // ����в������ֵ������Ϊ��ǰ�߶���Ͻ���
        if (residual > residual_threshold) {
            // ��¼��ǰ�ε���Ͻ��
            segments.push_back({ slope, intercept, current_x[0], current_x.back() });

            // ���õ�ǰ�߶�
            current_x = { x[i] };
            current_y = { y[i] };
        }

    }

    // ���һ�����
    double slope = 0.0;
    double intercept = 0.0;
    slope = (slope == 0.0) ? 0.0 : slope;
    intercept = (intercept == 0.0) ? 0.0 : intercept;
    segments.push_back({ slope, intercept, current_x[0], current_x.back() });

    // ������Ͻ��
    vector<float> xd(x.size());
    for (size_t i = 0; i < xd.size(); ++i) {
        xd[i] = x[i];
    }

    scatter(x, y, 10, { {"color", "gray"} });  // ����ԭʼ����

    vector<float> x_segment;
    vector<float> y_segment;

    // ����ÿһ����ϵ�ֱ��
    // ������������
    cout << "Intersection X-coordinates:\n";
    for (size_t i = 0; i < segments.size(); ++i) {
        float slope = get<0>(segments[i]);
        float intercept = get<1>(segments[i]);
        float x_start = get<2>(segments[i]);
        float x_end = get<3>(segments[i]);

        if (i != segments.size() - 1)
            cout << x_end;
        for (size_t j = 0; j < xd.size(); ++j) {
            if (xd[j] >= x_start && xd[j] <= x_end) {
                x_segment.push_back(xd[j]);
                y_segment.push_back(slope * xd[j] + intercept);
            }
        }
        
        plot(x_segment, y_segment, { {"label", "Segment " + to_string(i + 1)} });
    }

    legend();
    show();


    return x_segment;
}

// ����-���ģ�ͣ�C-A��
void Multifractal::CAmodel() {
    //��鴫��Ĳ�����Ч��
    if (ele_data.empty() || GroupNum < 3)
    {
        cout << "��������Ϊ�ջ���ֵ�ָ����̫�٣�����������!!" << endl;
    }

    float dCmin = *min_element(ele_data.begin(), ele_data.end());
    float dCmax = *max_element(ele_data.begin(), ele_data.end());

    float delta = (log(dCmax / dCmin)) / GroupNum;

    vector<float> refs, nc;
    //����ÿһ������ĺ���֮��
    for (int lev = 0; lev < GroupNum; lev++)
    {
        float Ci = dCmin * pow(_E, lev * delta);
        float Ni = 0;
        for (int cIdx = 0; cIdx < ele_data.size(); cIdx++)
        {
            if (ele_data[cIdx] >= Ci)
            {
                Ni += 1;
            }
        }
        refs.push_back(Ci);
        nc.push_back(Ni);
    }
    vector<float> intersections = segmentFitting(refs, nc, LineNum);
    //residual(refs, nc);

    const int rows = geo_data.size();
    const int cols = geo_data[0].size();
    vector<float> classified_data(rows * cols);

    // ����geo_data�����е�ÿ��Ԫ��
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            float value = geo_data[i][j];
            int category = 0;

            // ����intersection��ֵ�����ݽ��з���
            for (size_t k = 0; k < intersections.size(); ++k) {
                if (log(value) < intersections[k]) {
                    category = k + 1; // ���������Ϊ1, 2, ..., k
                    break;
                }
            }

            // ���ֵ���ڵ������һ���߽�ֵ����������һ��
            if (log(value) >= intersections.back()) {
                category = intersections.size() + 1;
            }

            // ���������洢��classified_data������
            classified_data[i * cols + j] = category;
        }
    }
    imshow(classified_data.data(), rows, cols, 1, {});
    show();
}

// ����-����ģ�ͣ�C-D��
void Multifractal::CDmodel() {
    
}

// ����-����ģ�ͣ�C-Q)
void Multifractal::CQmodel() {
    //��鴫��Ĳ�����Ч��
    if (ele_data.empty() || GroupNum < 3)
    {
        cout << "��������Ϊ�ջ���ֵ�ָ����̫�٣�����������!!" << endl;
    }

    float dCmin = *min_element(ele_data.begin(), ele_data.end());
    float dCmax = *max_element(ele_data.begin(), ele_data.end());

    float delta = (log(dCmax / dCmin)) / GroupNum;

    vector< float> refs, qc;
    //����ÿһ������ĺ���֮��
    for (int lev = 0; lev < GroupNum; lev++)
    {
        float Ci = dCmin * pow(_E, lev * delta);
        float Qi = 0;
        for (int cIdx = 0; cIdx < ele_data.size(); cIdx++)
        {
            if (ele_data[cIdx] >= Ci)
            {
                Qi += ele_data[cIdx];
            }
        }
        refs.push_back(Ci);
        qc.push_back(Qi / Ci);
    }
    vector<float> intersections = residual(refs, qc);

    const int rows = geo_data.size();
    const int cols = geo_data[0].size();
    vector<float> classified_data(rows * cols);

    // ����geo_data�����е�ÿ��Ԫ��
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            float value = geo_data[i][j];
            int category = 0;

            // ����intersection��ֵ�����ݽ��з���
            for (size_t k = 0; k < intersections.size(); ++k) {
                if (log(value) < intersections[k]) {
                    category = k + 1; // ���������Ϊ1, 2, ..., k
                    break;
                }
            }

            // ���ֵ���ڵ������һ���߽�ֵ����������һ��
            if (log(value) >= intersections.back()) {
                category = intersections.size() + 1;
            }

            // ���������洢��classified_data������
            classified_data[i * cols + j] = category;
        }
    }
    imshow(classified_data.data(), rows, cols, 1, { });
    show();
}

// Ƶ��-���ģ�ͣ�S-A��
void Multifractal::SAmodel() {
    // ���е�����ֵ
    iterative_interpolation(geo_data);

    // �����ά��������
    using Matrix = vector<vector<float>>;

    // ����Ҷ�任
    int rows = geo_data.size();
    int cols = geo_data[0].size();
    const int size = rows * cols; // ��Ԫ����

    fftw_complex* in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * rows * cols);
    fftw_complex* out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * rows * cols);
    fftw_plan plan = fftw_plan_dft_2d(rows, cols, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

    // ��ʼ����������
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            in[i * cols + j][0] = geo_data[i][j]; // ʵ��
            in[i * cols + j][1] = 0.0;            // �鲿
        }
    }

    // ִ�и���Ҷ�任
    fftw_execute(plan);

    // ��ȡ�������
    vector<vector<complex<float>>> W(rows, vector<complex<float>>(cols));//�洢����Ҷ�任���
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            W[i][j] = complex<float>(out[i * cols + j][0], out[i * cols + j][1]);
        }
    }

    // ����
    fftw_destroy_plan(plan);
    fftw_free(in);
    fftw_free(out);

    // ���㸴����ƽ��
    vector<vector<complex<float>>> S(rows, vector<complex<float>>(cols)); //�洢�����ܶ�
    vector<vector<complex<float>>> S1(rows, vector<complex<float>>(cols, complex<float>(0.0f, 0.0f))); //�洢����ֵ
    vector<vector<complex<float>>> S2(rows, vector<complex<float>>(cols, complex<float>(0.0f, 0.0f))); //�洢�쳣ֵ
    vector<vector<complex<float>>> S3(rows, vector<complex<float>>(cols, complex<float>(0.0f, 0.0f))); //�洢����ֵ

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            S[i][j] = W[i][j] * W[i][j]; // ���㸴����ƽ��
        }
    }


    // ����ÿ��������ģ���洢��һ��һά������
    vector<float> mods;  // �洢ÿ��������ģ����С��

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            // ʹ�� std::abs ���㸴����ģ
            mods.push_back(abs(S[i][j]));
        }
    }

    // ��ģ������������
    sort(mods.begin(), mods.end());


    float dCmin = *min_element(mods.begin(), mods.end());
    float dCmax = *max_element(mods.begin(), mods.end());

    float delta = (log(dCmax / dCmin)) / GroupNum;

    vector<float> refs, sc;
    //����ÿһ������ĺ���֮��
    for (int lev = 0; lev < GroupNum; lev++)
    {
        float Ci = dCmin * pow(_E, lev * delta);
        float Ni = 0;
        for (int cIdx = 0; cIdx < mods.size(); cIdx++)
        {
            if (mods[cIdx] >= Ci)
            {
                Ni += 1;
            }
        }
        refs.push_back(Ci);
        sc.push_back(Ni);
    }

    vector<float> intersections = segmentFitting(refs, sc, LineNum);

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            if (log(abs(S[i][j])) > intersections[LineNum - 2])
                S1[i][j] = S[i][j];
            else if (log(abs(S[i][j])) <= intersections[0])
                S3[i][j] = S[i][j];
            else
                S2[i][j] = S[i][j];
        }
    }

    // ��ʼ�� FFTW �������
    fftw_complex* input = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * size);
    fftw_complex* output = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * size);

    // ��������ά�渵��Ҷ�任
    auto inverse_fft = [&](vector<vector<complex<float>>>& src, vector<vector<complex<float>>>& dest) {
        // ��� FFTW �������飬���� geo_data ��Ӱ��
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                int idx = i * cols + j;
                input[idx][0] = src[i][j].real() * geo_data[i][j]; // ʵ��
                input[idx][1] = src[i][j].imag() * geo_data[i][j]; // �鲿
            }
        }

        // �����渵��Ҷ�任�ƻ�
        fftw_plan plan = fftw_plan_dft_2d(rows, cols, input, output, FFTW_BACKWARD, FFTW_ESTIMATE);

        // ִ�б任
        fftw_execute(plan);

        // ������洢��Ŀ������ dest ��
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                int idx = i * cols + j;
                dest[i][j] = complex<float>(output[idx][0], output[idx][1]);
            }
        }

        // ���ټƻ�
        fftw_destroy_plan(plan);
        };

    // ���ζ� S1��S2 �� S3 �����渵��Ҷ�任
    vector<vector<complex<float>>> B(rows, vector<complex<float>>(cols)); // ���� S1 ����任���
    vector<vector<complex<float>>> A(rows, vector<complex<float>>(cols)); // ���� S2 ����任���
    vector<vector<complex<float>>> Q(rows, vector<complex<float>>(cols)); // ���� S3 ����任���

    inverse_fft(S1, B);
    inverse_fft(S2, A);
    inverse_fft(S3, Q);

    // ��� B1 �Ľ��
    for (int i = 0; i < 5; ++i) {
        for (int j = 0; j < 5; ++j) {
            cout << B[i][j] << " ";
        }
        cout << endl;
    }

    // �ͷ� FFTW �ڴ�
    fftw_free(input);
    fftw_free(output);

    vector<float> classified_data(rows * cols);

    // ����geo_data�����е�ÿ��Ԫ��
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            float value = B[i][j].real();
            classified_data[i * cols + j] = value;
        }
    }
    imshow(classified_data.data(), rows, cols, 1, {});
    show();

}

// �ֲ������Է���
void Multifractal::calculateHeterogeneity() {
    int M = geo_data.size(); // ��ȡ����
    int N = geo_data[0].size(); //��ȡ����
    vector<float> classified_data; //��ά����ת�д洢����
    alpha_data.resize(M, vector<float>(N)); //�洢������ָ��

    MatrixXd geo_matrix(M, N);
    for (int i = 0; i < M; ++i) {
        for (int j = 0; j < N; ++j) {
            geo_matrix(i, j) = geo_data[i][j]; //����������ת�浽����������
        }
    }

    // ����դ�����ݵ�ÿһ����
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {

            // �洢��ǰ���ڲ�ͬ r �µ� r �� c ֵ
            vector<pair<float, float>> rc_values;

            // ��С����ѭ�����ڰ뾶 r
            for (float r = w_min; r <= w_max; r += w_step) {

                // ȷ����ǰ���ڵı߽�
                int row_min = max(0, i - static_cast<int>(r));
                int row_max = min(M - 1, i + static_cast<int>(r));
                int col_min = max(0, j - static_cast<int>(r));
                int col_max = min(N - 1, j + static_cast<int>(r));

                // ��ȡ�����ڵ�Ԫ�غ���������ƽ��ֵ
                MatrixXd window_data = geo_matrix.block(row_min, col_min, row_max - row_min + 1, col_max - col_min + 1);
                float c = window_data.mean(); // ������Ԫ�ص�ƽ������

                // ����ǰ�뾶 r �Ͷ�Ӧ�� c ֵ���� rc_values ����
                rc_values.push_back(make_pair(r, c));
            }

            // ���� log(r) �� log(c)
            vector<float> log_r, log_c;
            for (const auto& pair : rc_values) {
                log_r.push_back(log(pair.first));
                log_c.push_back(log(pair.second));
            }

            // ʹ����С���˷������������
            int n = log_r.size();
            if (n > 1) {
                VectorXd x(n), y(n);
                for (int k = 0; k < n; k++) {
                    x[k] = log_r[k];
                    y[k] = log_c[k];
                }

                // ������� X �� Y
                MatrixXd X(n, 2);
                X.col(0) = VectorXd::Ones(n); // ������
                X.col(1) = x;

                // ��С���˷���⣺ (X'X)^(-1) X'Y
                VectorXd coeffs = (X.transpose() * X).ldlt().solve(X.transpose() * y);

                // ��ȡб��
                double slope = coeffs[1];

                // �����б����Ϊ�õ��������ָ��
                alpha_data[i][j] = slope + 2;
            }
        }
    }

    // ����geo_data�����е�ÿ��Ԫ��
    for (int i = 0; i < M; ++i) {
        for (int j = 0; j < N; ++j) {
            float value = alpha_data[i][j];
            // ���������洢��classified_data������
            classified_data.push_back(value);
        }
    }
    imshow(classified_data.data(), M, N, 1, {});
    show();

}

// ��ʾ�˵�����ȡ�û�ѡ��
int Multifractal::getUserChoice() {
    int choice;
    while (true) {
        cout << "��ѡ���Ӧ���쳣ʶ���㷨��" << endl;
        cout << "================================" << endl;
        cout << "1. ����-���ģ�ͣ�C-A��\n";
        cout << "2. ����-����ģ�ͣ�C-D��\n";
        cout << "3. ����-����ģ�ͣ�C-Q��\n";
        cout << "4. Ƶ��-���ģ�ͣ�S-A��\n";
        cout << "5. �ֲ������Է���\n";
        cout << "================================" << endl;
        cout << "���������ѡ��: ";
        cin >> choice;

        if (cin.fail() || choice < 1 || choice > 5) { // �����Ч����
            cin.clear(); // ��������־
            cin.ignore(numeric_limits<streamsize>::max(), '\n'); // ����ʣ������
            cout << "��Ч���룬������������ȷ�ĵ����֡�\n";
            delay(1000); // �ӳ� 1000 ���루1 �룩
            clearScreen(); // ��ս���
        }
        else {
            cin.ignore(numeric_limits<streamsize>::max(), '\n'); // ������뻺��
            return choice; // ������Ч����
        }
    }
}

// �����û�������ж�Ӧ���ݴ���
void Multifractal::anomaly_detection(int choice) {
    if (choice > 0 && choice < 5) {
        cout << "��ѡ��ʹ�ö��ط��εķ����������쳣��ȡ!!" << endl;
        cout << "������ָ���ֵ��Ŀ GroupNum :";
        cin >> GroupNum;
        cout << "������ָ��߶���Ŀ LineNum :";
        cin >> LineNum;

        delay(1000); // �ӳ� 1000 ���루1 �룩
        clearScreen(); // ��ս���
    }
    else {
        cout << "��ѡ��ʹ�þֲ������Է��������������쳣��ȡ!!" << endl;
        cout << "�����봰����С�뾶 w_min :";
        cin >> w_min;
        cout << "�����봰�����뾶 w_max :";
        cin >> w_max;
        cout << "�����봰���ƶ����� w_step :";
        cin >> w_step;

        delay(1000); // �ӳ� 1000 ���루1 �룩
        clearScreen(); // ��ս���
    }
    switch (choice)
    {
    case 1:
        CAmodel();
        break;

    case 2:
        CDmodel();
        break;

    case 3:
        CQmodel();
        break;

    case 4:
        SAmodel();
        break;

    case 5:
        calculateHeterogeneity();
    default:
        return; // No fallthrough here, explicitly breaking
    }
    return;
}