#include "main.h"
#include <cmath>
#include <iostream>
#include "gnuplot-iostream-master/gnuplot-iostream.h"
#include "Eigen/Dense"

using namespace Eigen;

#define GP_DSAC1 "plot \"../DataSets/dsetA_class1\"\n" //GnuPlot_DataSetAClass1
#define GP_DSAC2 "plot \"../DataSets/dsetA_class2\"\n"
#define GP_DSBC1 "plot \"../DataSets/dsetB_class1\"\n"
#define GP_DSBC2 "plot \"../DataSets/dsetB_class2\"\n"
#define GP_DSA "plot for[i = 1:2] \"../DataSets/dsetA_class\".i.\"\" \n"
#define GP_DSB "plot for[i = 1:2] \"../DataSets/dsetB_class\".i.\"\" \n"
#define GP_DSA_DBOUND "plot  \"../DataSets/disc1\", for[i = 1:2] \"../DataSets/dsetA_class\".i.\"\" \n"
#define GP_DSA_CLASS_DBOUND " plot  \"../DataSets/dboundA\", for[i = 1:2] \"../DataSets/dsetA_classesfor\".i.\"\" \n"
#define GP_DSB_CLASS_DBOUND " plot  \"../DataSets/dboundB\", for[i = 1:2] \"../DataSets/dsetB_classesfor\".i.\"\" \n"
#define GP_DSA_SAMPLES_DBOUND " plot  \"../DataSets/dboundA\", for[i = 1:2] \"../DataSets/dsetA_class\".i.\"\" \n"
#define GP_DSB_SAMPLES_DBOUND " plot  \"../DataSets/dboundB\", for[i = 1:2] \"../DataSets/dsetB_class\".i.\"\" \n"

std::vector<std::vector<group>> prepare_datasets(bool writeToFile=false); //Write datasets to files, and return both datasets
void write_group_to_file(const group &A, const std::string &name);
void write_points_to_file(const std::vector<std::pair<double, double>> & points, const std::string & name);
double disc_case1(const group & A, const sample& observation );
double disc_case3(const group & B, const sample& observation );
double disc_euclidean(const group & A, const sample & observation);
std::vector<std::vector<std::pair<double, double>>> classify(const std::vector<group> & dset, int discase = 1);
std::vector<std::pair<double, double>> create_line_from_disc_case1(const std::vector<group> & dset);
std::vector<std::pair<double, double>> create_line_from_disc_case3(const std::vector<group> & dset);
std::vector<std::pair<double, double>> create_line_from_disc_case4(const std::vector<group> & dset);
double misclasses(const std::vector<std::vector<std::pair<double, double>>> & classifications, const std::vector<group> & dset);

int main() {
    std::vector<std::vector<group>> dsets = prepare_datasets(true);

    std::cout << "Hello! Would you like to use dataset A or B? (Note: Program may take up to 5 mins to compelte running due to extremely "
                 "inefficient search function): ";
    char input;
    std::cin >> input;

    std::cout << "Would you like to use a Euclidean distance classifier(y/n): ";
    char useE;
    std::cin >> useE;

    if(input == 'A') {
        auto discPoints = create_line_from_disc_case1(dsets[0]);

        //Currently supports discriminant cases 1 and 3. Case 4 is Euclidean distance
        auto classifications = classify(dsets[0], useE == 'y' ? 4 : 1);
        misclasses(classifications, dsets[0]);
        write_points_to_file(discPoints, "../DataSets/dboundA");
        write_points_to_file(classifications[0], "../DataSets/dsetA_classesfor1");
        write_points_to_file(classifications[1], "../DataSets/dsetA_classesfor2");

        Gnuplot gp;
        Gnuplot gp2;
        gp << GP_DSA_CLASS_DBOUND;
        gp2 << GP_DSA_SAMPLES_DBOUND;
//    std::vector<int> classes = classify_case1(dsets[0]);
//    Gnuplot gp2;
//    Gnuplot gp3;
//    gp << GP_DSA;
//    gp2 << GP_DSB;
    }

    else if (input == 'B'){
          auto discPoints = create_line_from_disc_case3(dsets[1]);
          auto classifications = classify(dsets[1], useE == 'y' ? 4 : 3);
          misclasses(classifications, dsets[1]);
          write_points_to_file(discPoints, "../DataSets/dboundB");
          write_points_to_file(classifications[0], "../DataSets/dsetB_classesfor1");
          write_points_to_file(classifications[1], "../DataSets/dsetB_classesfor2");

        Gnuplot gp;
        Gnuplot gp2;
        gp << GP_DSB_SAMPLES_DBOUND;
        gp2 << GP_DSB_CLASS_DBOUND;

    }
    else {
        fprintf(stderr, "Please input either A or B");
        exit(69);
    }

    return 0;
}
std::vector<std::vector<group>> prepare_datasets(bool writeToFile){
    std::vector<group> A;
    A.emplace_back(group(1, 1, 1, 1,60000, 200000));
    A.emplace_back(group(4,4,1,1,140000, 200000));



    std::vector<group> B;
    B.emplace_back(group(1,1,1,1,40000,200000));
    B.emplace_back(group(4,4,4,8,160000,200000));

    if(writeToFile) {
        write_group_to_file(A[0], "../DataSets/dsetA_class1");
        write_group_to_file(A[1], "../DataSets/dsetA_class2");
        write_group_to_file(B[0], "../DataSets/dsetB_class1");
        write_group_to_file(B[1], "../DataSets/dsetB_class2");
    }
    std::vector<std::vector<group>> dsets = {A, B};

    return dsets;
}

void write_group_to_file(const group & A, const std::string & name){
    std::ofstream outdata;
    outdata.open(name.c_str(), std::ios::trunc);

    std::vector<std::vector<double>> samples = A.samples;

    for (auto & sample : samples){
        outdata << sample[0] << " " << sample[1] << "\n";
    }
    outdata.close();

}

std::vector<std::vector<std::pair<double, double>>> classify(const std::vector<group> & dset, int discase){

    group group1 = dset[0];
    group group2 = dset[1];

    std::vector<sample> samples;
    samples.insert(samples.end(), group1.samples.begin(), group1.samples.end());
    samples.insert(samples.end(), group2.samples.begin(), group2.samples.end());

    std::vector<std::vector<std::pair<double, double>>> classifications(2);
        for (auto &sample: samples) {
            double g1, g2;
            if (discase == 1) {
                 g1 = disc_case1(group1, sample);
                 g2 = disc_case1(group2, sample);
            }
            else if (discase == 3){
                g1 = disc_case3(group1, sample);
                g2 = disc_case3(group2, sample);
            }
            else if (discase == 4){
                g1 = disc_euclidean(group1, sample);
                g2 = disc_euclidean(group2, sample);
            }
            double gx = g1 - g2;

            if (gx > 0){
                classifications[0].emplace_back(std::make_pair(sample[0], sample[1]));
            }
            else {
                classifications[1].emplace_back(std::make_pair(sample[0], sample[1]));
            }
        }


    return classifications ;
}
double disc_euclidean(const group & A, const sample & observation){
    VectorXd obs(2,1);
    VectorXd means(2,1);

    obs(0,0) = observation[0];
    obs(1,0) = observation[1];
    means(0,0) = A.meanx;
    means(1,0) = A.meany;

    auto result = obs -  means;
    return (- (result.norm() * result.norm()));
}
double disc_case1(const group & A, const sample& observation ){
    //Get wi and wi0
    VectorXd means(2,1);
    MatrixXd cov(2,2);

    means[0] = A.meanx;
    means[1] = A.meany;

    cov(0,0) = A.stdevx;
    cov(1,1) = A.stdevy;
    cov(0,1) = 0;
    cov(1,0) = 0;

    auto wi_vals = (1/(A.stdevx * A.stdevx)) * means;


    auto wi0 =  ( (-1.0/(2.0* A.stdevx * A.stdevx)) * means.transpose() * means) + std::log(A.prior);
    VectorXd wi(2,1);
    wi(0,0) = wi_vals(0,0);
    wi(1,0) = wi_vals(1,0);




    //Now, run our observation through the discriminant
    VectorXd x(2,1);
    x(0,0) = observation[0];
    x(1,0) = observation[1];
    double value = ((wi.transpose()).dot(x)) + wi0;

    return value;
}

std::vector<std::pair<double, double>> create_line_from_disc_case1(const std::vector<group> & dset) {

    std::vector<std::pair<VectorXd, double>> discrims;
    for (auto & A : dset) {
        //Get wi and wi0
        VectorXd means(2, 1);
        MatrixXd cov(2, 2);

        means[0] = A.meanx;
        means[1] = A.meany;

        cov(0, 0) = A.stdevx;
        cov(1, 1) = A.stdevy;
        cov(0, 1) = 0;
        cov(1, 0) = 0;
//
//    std::cout << "Means: " << means[0] << ", " << means[1] << "\n";
//    std::cout << "Cov:\n"
//                << cov(0,0) << " " << cov(0,1) << "\n"
//                << cov(1,0) << " " << cov(1,1) << "\n";
        auto wi_vals = (1.0 / (A.stdevx * A.stdevx)) * means;
        auto wi0 = (-1 * (1.0 / (2.0 * A.stdevx * A.stdevx)) * means.transpose() * means) + std::log(A.prior);
        VectorXd wi(2, 1);
        wi(0, 0) = wi_vals(0, 0);
        wi(1, 0) = wi_vals(1, 0);

        discrims.emplace_back(std::make_pair(wi, wi0));
    }

    //Math Time - first, set discriminants = to each other and isolate wi0 and wi
    double c = std::get<1>(discrims[1]) - std::get<1>(discrims[0]);
    auto wi = std::get<0>(discrims[0]) - std::get<0>(discrims[1]);

    //Then, set one side of equation to solve for x1 and x2. Get sample points
    std::vector<std::pair<double, double>> points;
    for (double x = -10; x < 10; x += 0.05f){
        double y = (c - (x * wi(0, 0))) / wi(1, 0);

        points.emplace_back(std::make_pair(x, y));
    }
    return points;
}
std::vector<std::pair<double, double>> create_line_from_disc_case3(const std::vector<group> & dset){

    std::vector<std::vector<double>> coefficients(2); //Represents coefficients of discriminant. Use ordered coefficients to solve for x2
    short count = 0;
    for (auto & g : dset) {
        VectorXd means(2,1);
        MatrixXd cov(2,2);

        means[0] = g.meanx;
        means[1] = g.meany;

        cov(0,0) = g.stdevx;
        cov(1,1) = g.stdevy;
        cov(0,1) = 0;
        cov(1,0) = 0;

        auto Wi = (-1.0/2) * cov.inverse();
        auto wi = cov.inverse()*means;
        //Convert matrix to array for ln, then back to cov
        double magcov = cov(0,0) * cov(1,1);

        auto term1 = (-1.0/2) * means.transpose() * cov.inverse() * means;
        auto term2 = (1.0/2) * std::log(magcov);
        auto term3 = std::log(g.prior);

        auto wi0 =term1 - term2 + term3;

        coefficients[count].emplace_back(Wi(0,0)); //Last coefficient x1^2
        coefficients[count].emplace_back(Wi(1,1)); //4th to last, x2^2
        coefficients[count].emplace_back(wi.transpose()(0,0)); //3rd to last, x1
        coefficients[count].emplace_back(wi.transpose()(0,1)); //2nd to last, x2, is the 2nd term of wit
        coefficients[count].emplace_back(wi0); //Last coefficient is constant, wi0

        count += 1;
    }

    //Now solve for x1
    std::vector<std::pair<double,double>> points;
    for (double x = -10; x < 10; x += 0.05){

        double c1, c2, c4, w10, k1,k2,k4,w20;
        c1 = coefficients[0][0];
        c2 = coefficients[0][1];
        c4 = coefficients[0][3];
        w10 = coefficients[0][4];
        k1 = coefficients[1][0];
        k2 = coefficients[1][1];
        k4 = coefficients[1][3];
        w20 = coefficients[1][4];
        double numerator =  ((k2-c2)*x*x) + ((k4 - c4) *x) +  w20 - w10;
        double denominator = c1 - k1;
        if ( (numerator/denominator) < 0 || denominator == 0){
            continue;
        }
        double y = sqrt(numerator/denominator);

        points.emplace_back(std::make_pair(x,y));
    }

    return points;

}

void write_points_to_file(const std::vector<std::pair<double, double>> & points, const std::string & name){
    std::ofstream f(name, std::ios::trunc);
    for (auto & point : points){
        f << std::get<0>(point) << " " << std::get<1>(point) << "\n";
    }
    f.close();
}

//I know it's O(n^2) but it's still computable in a kinda reasonable amount of time and the day is short and my legs have been sitting for the past 5 hours forgive me
double misclasses(const std::vector<std::vector<std::pair<double, double>>> & classifications, const std::vector<group> & dset){

    //Print misclassification for each class

    short count = 0;
    double totalNotFound = 0;
    double totalSamples = 0;

    for (auto & g : dset){
        double notFoundNum = 0;
        totalSamples += g.sampleNum;
        for (auto & sample : g.samples){
            std::pair s(sample[0], sample[1]);

            bool found = false;
            for (auto & val : classifications[count]){
                if (std::get<0>(s) == std::get<0>(val) && std::get<1>(s) == std::get<1>(val)){
                    found = true;
                    break;
                }
            }

            if(!found) notFoundNum++;
        }

        std::cout << "Misclass rate for class " << count << ": " << notFoundNum/g.sampleNum << "\n";
        totalNotFound += notFoundNum;
        count += 1;
    }

    std::cout << "Total misclass rate: " << totalNotFound/totalSamples;
    return 0;
}
double disc_case3(const group & B, const sample& observation ){
    VectorXd means(2,1);
    MatrixXd cov(2,2);

    means[0] = B.meanx;
    means[1] = B.meany;

    cov(0,0) = B.stdevx;
    cov(1,1) = B.stdevy;
    cov(0,1) = 0;
    cov(1,0) = 0;

    auto Wi = (-1.0/2) * cov.inverse();
    auto wi = cov.inverse()*means;
    //Convert matrix to array for ln, then back to cov
    double magcov = cov(0,0) * cov(1,1);

    auto term1 = (-1.0/2) * means.transpose() * cov.inverse() * means;
    auto term2 = (1.0/2) * std::log(magcov);
    auto term3 = std::log(B.prior);

    auto wi0 =term1 - term2 + term3;


    //Run our observation through:
    VectorXd obs(2,1);
    obs(0,0) = observation[0];
    obs(1,0) = observation[1];

    double val1 = obs.transpose() * Wi * obs;
    double val2 = wi.transpose() * obs;
    double val3 = wi0;


    double value = val1 + val2 + val3;
    return value;
}