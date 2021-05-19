#include <iostream>
#include <random>
#include <Eigen/Dense>
using namespace Eigen;
using namespace std;
const int D = 10;
const int population = 100;
const int generation = 500;
const double c1 = 0.5, c2 = 0.5;
const double rmax = 32, rmin = 32;
double Sphere(MatrixXd M){
    double re = (M.transpose()*M).sum();
    return re;
}
int main()
{
    srand(31);

    MatrixXd X = MatrixXd::Random(population,D);
    MatrixXd V = MatrixXd::Random(population,D);
    X = (X.array() - 0.5) * (rmax - rmin) + rmin;
    MatrixXd pid = X;
    VectorXd pid_v(population);
    double min = 1e8;
    int index = 0;


    for(int i=1;i<=generation;i++){

        for(int j=0;j<population;j++) {
            double tmp = Sphere(X.row(j));

            if(i == 1){
                pid_v(j) = tmp;
                pid.row(j) = X.row(j);
            }else if(tmp < pid_v(j)){ // 記錄個人最佳
                pid_v(j) = tmp;
                pid.row(j) = X.row(j);
            }
            if (tmp < min) {
                min = tmp;
                index = j;
            }
        }
        for(int j=0;j<population;j++) {
            double r1 = (double) rand() / (RAND_MAX + 1.0);
            double r2 = (double) rand() / (RAND_MAX + 1.0);

            V.row(j) += c1 * r1 * (pid.row(j) - X.row(j)) + c2 * r2 * (pid.row(index) - X.row(j));
            X.row(j) += V.row(j);
        }
        cout<<"generation:"<<i<<" "<<min<<endl;
    }

}
