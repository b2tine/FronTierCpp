#ifndef BVH_TESTER_H
#define BVH_TESTER_H

#include <BVH.h>


class BVH_Tester : public BVH
{
    private:
        

    public:

        BVH_Tester() = default;
        ~BVH_Tester() = default;

        BVH_Tester(const BVH_Tester&) = delete;
        BVH_Tester& operator=(const BVH_Tester&) = delete;
        BVH_Tester(BVH_Tester&&) = delete;
        BVH_Tester& operator=(BVH_Tester&&) = delete;


        void buildFromHseVector(std::vector<Hse*> hseList)
        {
            constructLeafNodes(hseList);
            buildHeirarchy();
        }

};



#endif
