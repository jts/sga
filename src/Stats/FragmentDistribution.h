#ifndef FRAGMENTDISTRIBUTION_H
#define FRAGMENTDISTRIBUTION_H

#include "IntDist.h"
#include <map>

typedef std::map<int, int> IntIntMap;

class FragmentDistribution
{
    public:

        //
        FragmentDistribution() {}
        void readFromFile(std::string filename);
        IntDist convertToIntDist(double trim = 1.0f);

        //
        void increment(int index, int by = 1);
        void findPRange(double p, int& min, int& max);

        //
        int getCount(int index) const;
        double getFreq(int index) const;

    private:
        size_t m_totalCount;
        IntIntMap m_counts;
};
#endif

