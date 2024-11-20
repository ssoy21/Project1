#include "functions.h"

using namespace matplotlibcpp;
using namespace std;

int main() {
    Multifractal Mf;

    Mf.clearScreen();

    // 提示输入数据
    Mf.getInputFromFileManager();

    Mf.delay(2000); // 延迟 1000 毫秒（1 秒）
    Mf.clearScreen(); // 清空界面

    // 进行流程计算
    Mf.processChoice();
    Mf.delay(1000); // 延迟 1000 毫秒（1 秒）
    Mf.clearScreen(); // 清空界面

    // 选择异常识别算法
    int choice = Mf.getUserChoice();
    Mf.anomaly_detection(choice);

    cout << "计算完成。\n";

    return 0;
}

