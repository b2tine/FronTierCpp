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

        //TRI* left;
        //TRI* right;

            //double tension {0.0};

        //FabricEdge() = default;

        FabricEdge(POINT* p1, POINT* p2);

        void setRestLength(double l);
        void setSpringConstant(double k);
        void setTearingThreshold(double T);

        double getLength() const;
        double getTension() const;

        bool checkForTear();

        void setWeakPointFlag(bool flag);

        void print() const;

        //Testing Functions
        bool checkForTearTest(int i);

    private:

        double ks {0.0};
        double tear_threshold {HUGE};
        double length0 {-1.0};
        
        double disp[3];
        double length {-1.0};
        double tension {0.0};

        void computeTension();
        bool underTension() const;
        void computeLength();
        void computeDisplacement();

        bool has_weakpt {false};
        double weakpt_factor {1.0}; //TODO: read from AF_PARAMS

        bool hasWeakPoint() const;
};


class FabricTearer
{
    public:

        ~FabricTearer();

        void collectFabricEdges(const INTERFACE* intfc);
        void setSpringData(double ks, double max_tension);
        
        std::vector<std::pair<long int, long int>>
            recordGindexPointPairs() const;
        
        std::vector<double> recordRestingEdgeLengths();

        std::vector<long int> recordGindexWeakPoints() const;
        
        void readGindexPointPairs(POINT** gpoints,
                const std::vector<std::pair<long int, long int>>& gindex_pairs);
        
        void readRestingEdgeLengths(const std::vector<double>& restlengths);

        void readGindexWeakPoints(
                const std::vector<long int>& gindex_weakpts);


        void tearFabric();
        //void tearFabric(Front* front);

        void printEdges() const;

        //Testing Functions
            //void setEdgeTension(int index, double T);
        void tearFabricTest();
        //void tearFabricTest(Front* front);

    private:
        
        std::vector<FabricEdge*> edges;
        std::vector<long int> tear_idx;
        std::unordered_set<long int> weakpt_idx;

        void clearEdges();
        
        void checkForTearingEvents();
        void processTearingEvents();

        void createNewTear(FabricEdge* e);
        void propagateTear(FabricEdge* e);

        bool isWeakPoint(POINT* p);

        //Testing Functions
        void checkForTearingEventsTest();
};


#endif
