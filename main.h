#ifndef CS479_PROJ1_MAIN_H
#define CS479_PROJ1_MAIN_H
#include "boxmuller.h"
#include <vector>

typedef std::vector<double> sample;

typedef struct group{

    int meanx = 0;
    int meany = 0;
    int stdevx = 0;
    int stdevy = 0;

    int sampleNum = 0;
    float prior = 0;
    std::vector<sample> samples;

    group(int mx, int my, int sx, int sy, float snum, float totnum) : meanx(mx), meany(my), stdevx(sx), stdevy(sy), sampleNum(snum){
        prior = snum/totnum;
        for (int i = 0; i < sampleNum; i++){
            double samplex = box_muller(meanx, stdevx);
            double sampley = box_muller(meany, stdevy);

            std::vector<double> sample;
            sample.push_back(samplex);
            sample.push_back(sampley);

            samples.push_back(sample);

        }
    };

} group;






#endif //CS479_PROJ1_MAIN_H
