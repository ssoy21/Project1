#include "functions.h"

using namespace matplotlibcpp;
using namespace std;

int main() {
    Multifractal Mf;

    Mf.clearScreen();

    // ��ʾ��������
    Mf.getInputFromFileManager();

    Mf.delay(2000); // �ӳ� 1000 ���루1 �룩
    Mf.clearScreen(); // ��ս���

    // �������̼���
    Mf.processChoice();
    Mf.delay(1000); // �ӳ� 1000 ���루1 �룩
    Mf.clearScreen(); // ��ս���

    // ѡ���쳣ʶ���㷨
    int choice = Mf.getUserChoice();
    Mf.anomaly_detection(choice);

    cout << "������ɡ�\n";

    return 0;
}

