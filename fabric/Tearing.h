#ifndef TEARING_H
#define TEARING_H

#include <FronTier.h>

#include <cassert>
#include <utility>
#include <unordered_set>
#include <vector>


class FabricEdge
{
    public:

        POINT* beg;
        POINT* end;

        double tension {0.0};

        //FabricEdge() = default;

        FabricEdge(POINT* p1, POINT* p2)
            : beg{p1}, end{p2}
        {}

        void setRestLength(double l);
        void setSpringConstant(double k);
        void setTearingThreshold(double T);

        double getLength();
        double getTension();

        void print();

        bool checkForTear();

        //Testing Functions
        bool checkForTearTest(int i);

    private:

        double ks {0.0};
        double tear_threshold {HUGE};
        
        double length0 {-1.0};
        double length {-1.0};
        double disp[3];

        void computeTension();
        bool underTension();
        void computeLength();
        void computeDisplacement();
};


class FabricTearer
{
    public:

        ~FabricTearer();

        void collectFabricEdges(const INTERFACE* intfc);
        void setSpringData(double ks, double max_tension);
        
        std::vector<std::pair<long int, long int>> recordGindexPairs();
        void readGindexPairs(POINT** gpoints,
                const std::vector<std::pair<long int, long int>>& gindexpairs);
        
        std::vector<double> recordRestingEdgeLengths();
        void readRestingEdgeLengths(const std::vector<double>& restlengths);

        void tearFabric();

        void printEdges();

        //Testing Functions
            //void setEdgeTension(int index, double T);
        void tearFabricTest();

    private:
        
        std::vector<long int> tear_idx;
        std::vector<FabricEdge*> edges;

        void clearEdges();
        
        void checkForTearingEvents();
        void processTearingEvents();

        //Testing Functions
        void checkForTearingEventsTest();
};


#endif
