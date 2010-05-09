#ifndef IntDist_H
#define IntDist_H

#include "Util.h"

//
// Probabilty distribution over a range of integer indices
//
class IntDist
{
    public:
        IntDist() {}
        IntDist(int start, int end);

        // get/set
        void addWeight(int pos, double w);
        void setP(int pos, double p);
        double getP(int index) const;
        double getWeight(int index) const;
        int getStart() const { return m_start; }
        int getEnd() const { return m_start + m_values.size(); }
        double expectedValue() const;

        //
        // normalize the distribution by calculating probabilities from weights
        //
        void normalize();

        // utility funcs
        size_t pos2Idx(int pos) const;
        friend std::ostream& operator<<(std::ostream& out, const IntDist& pdf);

    private:
        int m_start;
        bool m_normalized;
        DoubleVec m_values;
        static double ERROR_LIMIT;
};

#endif
