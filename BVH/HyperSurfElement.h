#ifndef HYPER_SURF_ELEMENT_H
#define HYPER_SURF_ELEMENT_H

#include <algorithm>
#include <memory>
#include <iostream>
#include <limits>
#include <exception>
#include <cassert>

#include <FronTier.h>


enum class HseTag
{
    FABRIC,
    STRING,
    RIGIDBODY,
    NONE
};


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

        virtual int num_pts() const = 0;
        virtual double max_coord(int i) const = 0;
        virtual double min_coord(int i) const = 0;
        virtual const POINT* const Point_of_hse(int i) const = 0;
        
        //void setTag(HseTag Tag);
        const HseTag getTag() const noexcept;
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
        
        int num_pts() const noexcept override { return 2; }
        double min_coord(int dim) const override;
        double max_coord(int dim) const override;
        const POINT* const Point_of_hse(int i) const override;
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
        
        int num_pts() const noexcept override { return 3; }
        double min_coord(int dim) const override;
        double max_coord(int dim) const override;
        const POINT* const Point_of_hse(int i) const override;
};

        

const bool areAdjacentHse(const Hse* const A, const Hse* const B);





#endif
