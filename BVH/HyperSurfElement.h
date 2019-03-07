#ifndef HYPER_SURF_ELEMENT_H
#define HYPER_SURF_ELEMENT_H

#include <algorithm>
#include <iostream>
#include <limits>
#include <exception>
#include <cassert>

#include <FronTier.h>

//TODO: Add method to check if a Hse (HsPoint, HsBond, or HsTri)
//      is incident to an instance of HsTri.


enum class HseTag {FABRIC,STRING,RIGIDBODY,NONE};


class Hse
{
    private:

        HseTag tag{HseTag::NONE};

    public:

        Hse(HseTag);
        Hse() = default;
        virtual ~Hse() = default;

        Hse(const Hse&) = delete;
        Hse& operator=(const Hse&) = delete;
        Hse(Hse&&) = delete;
        Hse& operator=(Hse&&) = delete;

        virtual double max_coord(int) const = 0;
        virtual double min_coord(int) const = 0;
        virtual POINT* Point_of_hse(int) const = 0;
        virtual int num_pts() const = 0;

        //void setTag(HseTag Tag);
        const HseTag getTag() const noexcept;
};


//Wrapper for FronTier POINT
class HsPoint : public Hse
{
    private:

        POINT* point{nullptr};

    public:

        explicit HsPoint(POINT*);
        HsPoint(POINT*, HseTag);
        ~HsPoint() = default;

        HsPoint() = delete;
        HsPoint(const HsPoint&) = delete;
        HsPoint& operator=(const HsPoint&) = delete;
        HsPoint(HsPoint&&) = delete;
        HsPoint& operator=(HsPoint&&) = delete;

        POINT* Point_of_hse(int i = 0) const override;
        double min_coord(int) const override;
        double max_coord(int) const override;
        int num_pts() const noexcept override { return 1; }
};


//Wrapper for FronTier BOND
class HsBond : public Hse
{
    private:

        BOND* bond{nullptr};

    public:

        explicit HsBond(BOND*);
        HsBond(BOND*, HseTag);
        ~HsBond() = default;

        HsBond() = delete;
        HsBond(const HsBond&) = delete;
        HsBond& operator=(const HsBond&) = delete;
        HsBond(HsBond&&) = delete;
        HsBond& operator=(HsBond&&) = delete;
        
        POINT* Point_of_hse(int) const override;
        double min_coord(int) const override;
        double max_coord(int) const override;
        int num_pts() const noexcept override { return 2; }
};


//Wrapper for FronTier TRI
class HsTri : public Hse
{
    private:

        TRI* tri{nullptr};

    public:

        explicit HsTri(TRI*);
        HsTri(TRI*, HseTag);
        ~HsTri() = default;

        HsTri() = delete;
        HsTri(const HsTri&) = delete;
        HsTri& operator=(const HsTri&) = delete;
        HsTri(HsTri&&) = delete;
        HsTri& operator=(HsTri&&) = delete;
        
        POINT* Point_of_hse(int) const override;
        double min_coord(int) const override;
        double max_coord(int) const override;
        int num_pts() const noexcept override { return 3; }
};




#endif
