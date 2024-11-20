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

// 清空屏幕
void Multifractal::clearScreen() {
    cout << "\033[2J\033[1;1H"; // 清空屏幕，适用于支持 ANSI 转义序列的终端
}

// 延迟时间
void Multifractal::delay(int milliseconds) {
    clock_t end_time = clock() + milliseconds * CLOCKS_PER_SEC / 1000; // 计算结束时间
    while (clock() < end_time) {
        // 这里可以选择执行其他代码或者保持空循环
    }
}

// 模拟文件管理器数据输入
void Multifractal::getInputFromFileManager() {
    string filePath = "D:\\test_data\\cu_ok.tif";
    string filePath_mask = "D:\\test_data\\cu_ok_raster1.tif";
    while (true) { // 循环直到成功输入文件路径
        cout << "请输入数据文件路径：" << filePath << endl;
        //getline(cin, filePath); // 获取用户输入的文件路径

        ifstream file(filePath); // 尝试打开文件

        // 检查文件是否成功打开
        if (!file) {
            cout << "无法打开数据文件，请检查文件路径是否正确!!" << endl;
            delay(1000);
            clearScreen(); // 清空界面
            // 这里可以选择添加延迟，如果需要
            continue; // 继续循环，重新提示输入文件路径
        }
        else {
            cout << "数据文件读取成功!!" << endl;
        }

        file.close(); // 关闭文件

        cout << "请输入掩码文件路径：" << filePath_mask << endl;
        //getline(cin, filePath); // 获取用户输入的文件路径

        ifstream file_mask(filePath_mask); // 尝试打开文件

        // 检查文件是否成功打开
        if (!file_mask) {
            cout << "无法打开掩码文件，请检查文件路径是否正确!!" << endl;
            delay(1000);
            clearScreen(); // 清空界面
            // 这里可以选择添加延迟，如果需要
            continue; // 继续循环，重新提示输入文件路径
        }
        else {
            cout << "掩码文件读取成功!!" << endl;
        }

        file_mask.close(); // 关闭文件

        //读取文件的部分
        readTifFile(filePath, filePath_mask);
        break; // 成功读取文件后跳出循环
    }
}

// 函数：读取 .tif 文件数据
bool Multifractal::readTifFile(const string& filePath, const string& filePath_mask) {
    GDALAllRegister(); // 注册所有格式

    // 打开文件
    GDALDataset* dataset = (GDALDataset*)GDALOpen(filePath.c_str(), GA_ReadOnly);
    GDALDataset* dataset_mask = (GDALDataset*)GDALOpen(filePath_mask.c_str(), GA_ReadOnly);

    if (!dataset) {
        cerr << "无法打开数据文件：" << filePath << endl;
        return false;
    }
    if (!dataset_mask) {
        cerr << "无法打开掩码文件：" << filePath_mask << endl;
        return false;
    }

    int xSize = dataset->GetRasterXSize(); // 获取横坐标大小
    int ySize = dataset->GetRasterYSize(); // 获取纵坐标大小
    float* rasterData = new float[xSize * ySize]; // 临时数组

    // 调整 geo_data 的大小
    geo_data.resize(ySize, vector<float>(xSize));

    // 读取数据
    GDALRasterBand* band = dataset->GetRasterBand(1); // 获取第一个波段
    CPLErr err = band->RasterIO(GF_Read, 0, 0, xSize, ySize, rasterData, xSize, ySize, GDT_Float32, 0, 0);


    if (err != CE_None) {
        cerr << "读取数据失败。" << endl;
        delete[] rasterData; // 释放临时数组
        GDALClose(dataset); // 关闭文件
        return false;
    }

    // 将数据存储到 geo_data，并剔除空值
    for (int y = 0; y < ySize; ++y) {
        for (int x = 0; x < xSize; ++x) {
            float value = rasterData[y * xSize + x];
            geo_data[y][x] = value; // 存储有效数据

        }
    }

    // 调整 geo_data 的大小
    mask_data.resize(ySize, vector<float>(xSize));

    // 读取数据

    band = dataset_mask->GetRasterBand(1); // 获取第一个波段
    err = band->RasterIO(GF_Read, 0, 0, xSize, ySize, rasterData, xSize, ySize, GDT_Float32, 0, 0);


    if (err != CE_None) {
        cerr << "读取数据失败。" << endl;
        delete[] rasterData; // 释放临时数组
        GDALClose(dataset_mask); // 关闭文件
        return false;
    }
    // 将数据存储到 geo_data，并剔除空值
    for (int y = 0; y < ySize; ++y) {
        for (int x = 0; x < xSize; ++x) {
            float value = rasterData[y * xSize + x];
            mask_data[y][x] = value; // 存储有效数据
        }
    }

    int nan_data = 0;
    //剔除无效值
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

    //对含量数据进行从小到大排序
    sort(ele_data.begin(), ele_data.end());
    cout << "空值有" << nan_data << "个!!" << endl;

    delete[] rasterData; // 释放临时数组

    GDALClose(dataset); // 关闭文件
    GDALClose(dataset_mask);

    return true;
}

// 检查是否是 NaN
bool is_nan(float val) {
    return isnan(val);
}

// 插值：计算邻域的平均值
void interpolate_once(vector<vector<float>>& data) {
    int rows = data.size();
    int cols = data[0].size();

    // 创建一个副本来存储插值结果
    vector<vector<float>> copy = data;

    // 遍历整个数组
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            if (is_nan(data[i][j])) {
                float sum = 0.0f;
                int count = 0;

                // 检查8邻域
                for (int di = -1; di <= 1; ++di) {
                    for (int dj = -1; dj <= 1; ++dj) {
                        int ni = i + di, nj = j + dj;
                        if (ni >= 0 && ni < rows && nj >= 0 && nj < cols && !is_nan(data[ni][nj])) {
                            sum += data[ni][nj];
                            count++;
                        }
                    }
                }

                // 如果邻域中有有效值，计算平均值填补
                if (count > 0) {
                    copy[i][j] = sum / count;
                }
            }
        }
    }

    // 更新原数组
    data = copy;
}

// 迭代插值方法
void iterative_interpolation(vector<vector<float>>& data, int max_iter = 100, float tolerance = 1e-6) {
    int rows = data.size();
    int cols = data[0].size();

    for (int iter = 0; iter < max_iter; ++iter) {
        vector<vector<float>> previous = data;

        // 执行一次插值
        interpolate_once(data);

        // 检查是否收敛
        float max_diff = 0.0f;
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                max_diff = max(max_diff, abs(data[i][j] - previous[i][j]));
            }
        }

        // 如果变化小于容忍度，认为插值已经收敛
        if (max_diff < tolerance) {
            cout << "插值收敛，迭代次数：" << iter + 1 << endl;
            break;
        }
    }
}

// 多元统计分析
void Multifractal::calculateStatistics() {
    // 均值
    ele_mean = accumulate(ele_data.begin(), ele_data.end(), 0.0) / ele_data.size();
    cout << "均值：" << ele_mean << endl;
    delay(1000);

    // 中位数
    vector<float> sortedData = ele_data;
    sort(sortedData.begin(), sortedData.end());
    size_t n = sortedData.size();
    ele_median = (n % 2 == 0) ? (sortedData[n / 2 - 1] + sortedData[n / 2]) / 2.0 : sortedData[n / 2];
    cout << "中位数：" << ele_median << endl;
    delay(1000);

    // 方差
    ele_variance = 0.0;
    for (float x : ele_data) {
        ele_variance += (x - ele_mean) * (x - ele_mean);
    }
    ele_variance /= ele_data.size();
    cout << "方差：" << ele_variance << endl;
    delay(1000);

    // 标准差
    ele_stdDev = sqrt(ele_variance);
    cout << "标准差：" << ele_stdDev << endl;
    delay(1000);

    // 偏度
    ele_skewness = 0.0;
    for (float x : ele_data) {
        ele_skewness += pow((x - ele_mean) / ele_stdDev, 3);
    }
    ele_skewness /= ele_data.size();
    cout << "偏度：" << ele_skewness << endl;
    delay(1000);

    // 峰度（超峰度）
    ele_kurtosis = 0.0;
    for (float x : ele_data) {
        ele_kurtosis += pow((x - ele_mean) / ele_stdDev, 4);
    }
    ele_kurtosis = ele_kurtosis / ele_data.size() - 3.0;
    cout << "峰度：" << ele_kurtosis << endl;
    delay(1000);

    // 最大值和最小值
    ele_maxVal = *max_element(ele_data.begin(), ele_data.end());
    ele_minVal = *min_element(ele_data.begin(), ele_data.end());
    cout << "最大值：" << ele_maxVal << endl;
    cout << "最小值：" << ele_minVal << endl;
    delay(1000);
}

// 对地球化学数据进行相应的流程计算
void Multifractal::processChoice() {
    cout << "数据多元统计分析：" << endl;
    cout << "================================" << endl;
    calculateStatistics();
    cout << "================================" << endl;

    cout << "";


}

// 线性拟合函数，返回斜率和截距
pair<float, float> linearFit(const vector<float>& x, const vector<float>& y) {
    int n = x.size();
    MatrixXd A(n, 2);
    VectorXd Y(n);
    for (int i = 0; i < n; ++i) {
        A(i, 0) = x[i];
        A(i, 1) = 1.0;  // 截距
        Y(i) = y[i];
    }

    // 使用QR分解求解 Ax = Y
    Vector2d coeffs = A.colPivHouseholderQr().solve(Y);
    float slope = coeffs(0);
    float intercept = coeffs(1);

    return { slope, intercept };
}

// 最小二乘法分段拟合函数
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
    vector<string> colors = { "r", "g", "b", "m", "c", "y", "k" }; // 可以根据需要增加更多颜色

    vector<float> slopes(k), intercepts(k);
    vector<float> intersection_x(k - 1);
    vector<pair<float, float>> intersections;  // 存储交点 (x, y)

    // 分段最小二乘拟合
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

        // 获取每段的斜率和截距
        auto [slope, intercept] = linearFit(x_segment, y_segment);
        slopes[i] = slope;
        intercepts[i] = intercept;

        // 计算线段的起点和终点
        float x_start = x_segment.front();
        float x_end = x_segment.back();
        float y_start = slope * x_start + intercept;
        float y_end = slope * x_end + intercept;

        // 绘制每段拟合线
        plot({ x_start, x_end }, { y_start, y_end }, colors[i % colors.size()] + "-"); // 使用不同颜色的拟合线
        
        // 如果不是最后一段，计算相邻段的重叠点
        if (i < k - 1) {
            float x_overlap = x_segment.back();  // 当前段的末尾点
            float y1 = slopes[i] * x_overlap + intercepts[i];       // 当前段的 y 值
            float y2 = slopes[i + 1] * x_overlap + intercepts[i + 1];  // 下一段的 y 值

            intersections.push_back({ x_overlap, y1 });  // 或仅记录 (x_overlap, y2)
        }
    }


    // 绘制交点
    vector<float> section_x, section_y;
    // 输出交点横坐标
    cout << "Intersection X-coordinates:\n";
    for (const auto& intersection : intersections) {
        section_x.push_back(intersection.first);
        section_y.push_back(intersection.second);
        cout << intersection.first << "  ";

    }
    scatter(section_x, section_y, 10.0, { {"color", "red"} });  // 参数调整
    show();

    return section_x;
}

// 循环拟合
vector<float> residual(const vector<float>& refs, const vector<float>& nc) {
    // x,y数组
    vector<float> x, y;
    vector<string> colors = { "r", "g", "b", "m", "c", "y", "k" }; // 可以根据需要增加更多颜色

    // x,y坐标的转换
    for (int i = 0; i < refs.size(); i++) {
        x.push_back(log(refs[i]));
        y.push_back(log(nc[i]));
    }

    // 设置拟合参数
    float residual_threshold = 0.2;  // 拟合点的最大残差
    int min_points_per_segment = 6;   // 每段拟合的最小点数

    vector<tuple<float, float, float, float>> segments;  // 存储每段的拟合参数
    vector<float> current_x = { x[0] };  // 当前线段的 x 数据
    vector<float> current_y = { y[0] };  // 当前线段的 y 数据

    // 遍历所有数据点进行拟合
    for (size_t i = 1; i < x.size(); ++i) {
        current_x.push_back(x[i]);
        current_y.push_back(y[i]);

        // 使用 Eigen 进行线性回归拟合
        VectorXd X(current_x.size());
        VectorXd Y(current_y.size());
        for (size_t j = 0; j < current_x.size(); ++j) {
            X(j) = current_x[j];
            Y(j) = current_y[j];
        }

        // 线性回归：y = mx + b，使用最小二乘法拟合
        MatrixXd A(current_x.size(), 2);
        for (size_t j = 0; j < current_x.size(); ++j) {
            A(j, 0) = current_x[j];  // x
            A(j, 1) = 1;              // 1 (常数项)
        }
        VectorXd coeff = A.colPivHouseholderQr().solve(Y);

        float slope = coeff(0);  // 斜率
        float intercept = coeff(1);  // 截距

        // 计算当前点的预测值
        float predicted_y = slope * x[i] + intercept;

        // 计算残差
        float residual = abs(y[i] - predicted_y);

        // 如果残差大于阈值，则认为当前线段拟合结束
        if (residual > residual_threshold) {
            // 记录当前段的拟合结果
            segments.push_back({ slope, intercept, current_x[0], current_x.back() });

            // 重置当前线段
            current_x = { x[i] };
            current_y = { y[i] };
        }

    }

    // 最后一段拟合
    double slope = 0.0;
    double intercept = 0.0;
    slope = (slope == 0.0) ? 0.0 : slope;
    intercept = (intercept == 0.0) ? 0.0 : intercept;
    segments.push_back({ slope, intercept, current_x[0], current_x.back() });

    // 绘制拟合结果
    vector<float> xd(x.size());
    for (size_t i = 0; i < xd.size(); ++i) {
        xd[i] = x[i];
    }

    scatter(x, y, 10, { {"color", "gray"} });  // 绘制原始数据

    vector<float> x_segment;
    vector<float> y_segment;

    // 绘制每一段拟合的直线
    // 输出交点横坐标
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

// 含量-面积模型（C-A）
void Multifractal::CAmodel() {
    //检查传入的参数有效性
    if (ele_data.empty() || GroupNum < 3)
    {
        cout << "输入数据为空或阈值分割个数太少，请重新输入!!" << endl;
    }

    float dCmin = *min_element(ele_data.begin(), ele_data.end());
    float dCmax = *max_element(ele_data.begin(), ele_data.end());

    float delta = (log(dCmax / dCmin)) / GroupNum;

    vector<float> refs, nc;
    //计算每一个级别的含量之和
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

    // 遍历geo_data数组中的每个元素
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            float value = geo_data[i][j];
            int category = 0;

            // 根据intersection的值对数据进行分类
            for (size_t k = 0; k < intersections.size(); ++k) {
                if (log(value) < intersections[k]) {
                    category = k + 1; // 将类别设置为1, 2, ..., k
                    break;
                }
            }

            // 如果值大于等于最后一个边界值，则归入最后一类
            if (log(value) >= intersections.back()) {
                category = intersections.size() + 1;
            }

            // 将分类结果存储在classified_data数组中
            classified_data[i * cols + j] = category;
        }
    }
    imshow(classified_data.data(), rows, cols, 1, {});
    show();
}

// 含量-距离模型（C-D）
void Multifractal::CDmodel() {
    
}

// 含量-总量模型（C-Q)
void Multifractal::CQmodel() {
    //检查传入的参数有效性
    if (ele_data.empty() || GroupNum < 3)
    {
        cout << "输入数据为空或阈值分割个数太少，请重新输入!!" << endl;
    }

    float dCmin = *min_element(ele_data.begin(), ele_data.end());
    float dCmax = *max_element(ele_data.begin(), ele_data.end());

    float delta = (log(dCmax / dCmin)) / GroupNum;

    vector< float> refs, qc;
    //计算每一个级别的含量之和
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

    // 遍历geo_data数组中的每个元素
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            float value = geo_data[i][j];
            int category = 0;

            // 根据intersection的值对数据进行分类
            for (size_t k = 0; k < intersections.size(); ++k) {
                if (log(value) < intersections[k]) {
                    category = k + 1; // 将类别设置为1, 2, ..., k
                    break;
                }
            }

            // 如果值大于等于最后一个边界值，则归入最后一类
            if (log(value) >= intersections.back()) {
                category = intersections.size() + 1;
            }

            // 将分类结果存储在classified_data数组中
            classified_data[i * cols + j] = category;
        }
    }
    imshow(classified_data.data(), rows, cols, 1, { });
    show();
}

// 频谱-面积模型（S-A）
void Multifractal::SAmodel() {
    // 进行迭代插值
    iterative_interpolation(geo_data);

    // 定义二维数据类型
    using Matrix = vector<vector<float>>;

    // 傅里叶变换
    int rows = geo_data.size();
    int cols = geo_data[0].size();
    const int size = rows * cols; // 总元素数

    fftw_complex* in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * rows * cols);
    fftw_complex* out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * rows * cols);
    fftw_plan plan = fftw_plan_dft_2d(rows, cols, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

    // 初始化输入数据
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            in[i * cols + j][0] = geo_data[i][j]; // 实部
            in[i * cols + j][1] = 0.0;            // 虚部
        }
    }

    // 执行傅里叶变换
    fftw_execute(plan);

    // 提取输出数据
    vector<vector<complex<float>>> W(rows, vector<complex<float>>(cols));//存储傅里叶变换结果
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            W[i][j] = complex<float>(out[i * cols + j][0], out[i * cols + j][1]);
        }
    }

    // 清理
    fftw_destroy_plan(plan);
    fftw_free(in);
    fftw_free(out);

    // 计算复数的平方
    vector<vector<complex<float>>> S(rows, vector<complex<float>>(cols)); //存储能谱密度
    vector<vector<complex<float>>> S1(rows, vector<complex<float>>(cols, complex<float>(0.0f, 0.0f))); //存储背景值
    vector<vector<complex<float>>> S2(rows, vector<complex<float>>(cols, complex<float>(0.0f, 0.0f))); //存储异常值
    vector<vector<complex<float>>> S3(rows, vector<complex<float>>(cols, complex<float>(0.0f, 0.0f))); //存储噪声值

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            S[i][j] = W[i][j] * W[i][j]; // 计算复数的平方
        }
    }


    // 计算每个复数的模并存储在一个一维向量中
    vector<float> mods;  // 存储每个复数的模（大小）

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            // 使用 std::abs 计算复数的模
            mods.push_back(abs(S[i][j]));
        }
    }

    // 对模进行升序排序
    sort(mods.begin(), mods.end());


    float dCmin = *min_element(mods.begin(), mods.end());
    float dCmax = *max_element(mods.begin(), mods.end());

    float delta = (log(dCmax / dCmin)) / GroupNum;

    vector<float> refs, sc;
    //计算每一个级别的含量之和
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

    // 初始化 FFTW 输入输出
    fftw_complex* input = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * size);
    fftw_complex* output = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * size);

    // 函数：二维逆傅里叶变换
    auto inverse_fft = [&](vector<vector<complex<float>>>& src, vector<vector<complex<float>>>& dest) {
        // 填充 FFTW 输入数组，考虑 geo_data 的影响
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                int idx = i * cols + j;
                input[idx][0] = src[i][j].real() * geo_data[i][j]; // 实部
                input[idx][1] = src[i][j].imag() * geo_data[i][j]; // 虚部
            }
        }

        // 创建逆傅里叶变换计划
        fftw_plan plan = fftw_plan_dft_2d(rows, cols, input, output, FFTW_BACKWARD, FFTW_ESTIMATE);

        // 执行变换
        fftw_execute(plan);

        // 将结果存储到目标数组 dest 中
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                int idx = i * cols + j;
                dest[i][j] = complex<float>(output[idx][0], output[idx][1]);
            }
        }

        // 销毁计划
        fftw_destroy_plan(plan);
        };

    // 依次对 S1、S2 和 S3 进行逆傅里叶变换
    vector<vector<complex<float>>> B(rows, vector<complex<float>>(cols)); // 储存 S1 的逆变换结果
    vector<vector<complex<float>>> A(rows, vector<complex<float>>(cols)); // 储存 S2 的逆变换结果
    vector<vector<complex<float>>> Q(rows, vector<complex<float>>(cols)); // 储存 S3 的逆变换结果

    inverse_fft(S1, B);
    inverse_fft(S2, A);
    inverse_fft(S3, Q);

    // 输出 B1 的结果
    for (int i = 0; i < 5; ++i) {
        for (int j = 0; j < 5; ++j) {
            cout << B[i][j] << " ";
        }
        cout << endl;
    }

    // 释放 FFTW 内存
    fftw_free(input);
    fftw_free(output);

    vector<float> classified_data(rows * cols);

    // 遍历geo_data数组中的每个元素
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            float value = B[i][j].real();
            classified_data[i * cols + j] = value;
        }
    }
    imshow(classified_data.data(), rows, cols, 1, {});
    show();

}

// 局部奇异性分析
void Multifractal::calculateHeterogeneity() {
    int M = geo_data.size(); // 获取行数
    int N = geo_data[0].size(); //获取列数
    vector<float> classified_data; //二维数组转列存储数据
    alpha_data.resize(M, vector<float>(N)); //存储奇异性指数

    MatrixXd geo_matrix(M, N);
    for (int i = 0; i < M; ++i) {
        for (int j = 0; j < N; ++j) {
            geo_matrix(i, j) = geo_data[i][j]; //将含量数据转存到矩阵数组中
        }
    }

    // 遍历栅格数据的每一个点
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {

            // 存储当前点在不同 r 下的 r 和 c 值
            vector<pair<float, float>> rc_values;

            // 从小到大循环窗口半径 r
            for (float r = w_min; r <= w_max; r += w_step) {

                // 确定当前窗口的边界
                int row_min = max(0, i - static_cast<int>(r));
                int row_max = min(M - 1, i + static_cast<int>(r));
                int col_min = max(0, j - static_cast<int>(r));
                int col_max = min(N - 1, j + static_cast<int>(r));

                // 提取窗口内的元素含量并计算平均值
                MatrixXd window_data = geo_matrix.block(row_min, col_min, row_max - row_min + 1, col_max - col_min + 1);
                float c = window_data.mean(); // 窗口内元素的平均含量

                // 将当前半径 r 和对应的 c 值存入 rc_values 数组
                rc_values.push_back(make_pair(r, c));
            }

            // 计算 log(r) 和 log(c)
            vector<float> log_r, log_c;
            for (const auto& pair : rc_values) {
                log_r.push_back(log(pair.first));
                log_c.push_back(log(pair.second));
            }

            // 使用最小二乘法进行线性拟合
            int n = log_r.size();
            if (n > 1) {
                VectorXd x(n), y(n);
                for (int k = 0; k < n; k++) {
                    x[k] = log_r[k];
                    y[k] = log_c[k];
                }

                // 构造矩阵 X 和 Y
                MatrixXd X(n, 2);
                X.col(0) = VectorXd::Ones(n); // 常数项
                X.col(1) = x;

                // 最小二乘法求解： (X'X)^(-1) X'Y
                VectorXd coeffs = (X.transpose() * X).ldlt().solve(X.transpose() * y);

                // 获取斜率
                double slope = coeffs[1];

                // 将拟合斜率作为该点的异质性指数
                alpha_data[i][j] = slope + 2;
            }
        }
    }

    // 遍历geo_data数组中的每个元素
    for (int i = 0; i < M; ++i) {
        for (int j = 0; j < N; ++j) {
            float value = alpha_data[i][j];
            // 将分类结果存储在classified_data数组中
            classified_data.push_back(value);
        }
    }
    imshow(classified_data.data(), M, N, 1, {});
    show();

}

// 显示菜单并获取用户选择
int Multifractal::getUserChoice() {
    int choice;
    while (true) {
        cout << "请选择对应的异常识别算法：" << endl;
        cout << "================================" << endl;
        cout << "1. 含量-面积模型（C-A）\n";
        cout << "2. 含量-距离模型（C-D）\n";
        cout << "3. 含量-总量模型（C-Q）\n";
        cout << "4. 频谱-面积模型（S-A）\n";
        cout << "5. 局部奇异性分析\n";
        cout << "================================" << endl;
        cout << "请输入你的选择: ";
        cin >> choice;

        if (cin.fail() || choice < 1 || choice > 5) { // 检查无效输入
            cin.clear(); // 清除错误标志
            cin.ignore(numeric_limits<streamsize>::max(), '\n'); // 忽略剩余输入
            cout << "无效输入，请重新输入正确的的数字。\n";
            delay(1000); // 延迟 1000 毫秒（1 秒）
            clearScreen(); // 清空界面
        }
        else {
            cin.ignore(numeric_limits<streamsize>::max(), '\n'); // 清除输入缓存
            return choice; // 返回有效输入
        }
    }
}

// 根据用户输入进行对应数据处理
void Multifractal::anomaly_detection(int choice) {
    if (choice > 0 && choice < 5) {
        cout << "您选择使用多重分形的方法来进行异常提取!!" << endl;
        cout << "请输入分割阈值数目 GroupNum :";
        cin >> GroupNum;
        cout << "请输入分割线段数目 LineNum :";
        cin >> LineNum;

        delay(1000); // 延迟 1000 毫秒（1 秒）
        clearScreen(); // 清空界面
    }
    else {
        cout << "您选择使用局部奇异性分析方法来进行异常提取!!" << endl;
        cout << "请输入窗口最小半径 w_min :";
        cin >> w_min;
        cout << "请输入窗口最大半径 w_max :";
        cin >> w_max;
        cout << "请输入窗口移动步长 w_step :";
        cin >> w_step;

        delay(1000); // 延迟 1000 毫秒（1 秒）
        clearScreen(); // 清空界面
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