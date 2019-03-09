#include "FrontWrapper.h"

//Initialization


//TODO: The function pointers that FT_StartUp() hooks up
//      should be virtual methods with default implementations.
//      The user can then inherit from the FrontWrapper class
//      and override the default methods.
void FrontWrapper::StartUp(F_BASIC_DATA* ft_basic)
{
    ResetTime();
    FT_StartUp(&cfront,ft_basic);
}

void FrontWrapper::ResetTime()
{
    FT_ResetTime(&cfront);
}

void FrontWrapper::InitIntfc(LEVEL_FUNC_PACK* level_func_pack)
{
    FT_InitIntfc(&cfront,level_func_pack);
    ClipIntfcToSubdomain();
}

//TODO: Why is this never used for 3d?
//      Also, documention in fapi.h suggests
//      this is only for the initial interface;
//      that's not obvious from the name.
void FrontWrapper::ClipIntfcToSubdomain()
{
    if( Dim() < 3 )
        FT_ClipIntfcToSubdomain(&cfront);
}

/*
void FrontWrapper::ReadTimeControl(std::string in_name)
{
    FT_ReadTimeControl(in_name.c_str());
}
*/

void FrontWrapper::InitFrontVeloFunc(VELO_FUNC_PACK* velo_func_pack)
{
    FT_InitFrontVeloFunc(&cfront,velo_func_pack);
}


//Temporary manual setters

//TODO: These should go in a separate testing subclass.
//      Production code should parse from input file only.
void FrontWrapper::setMaxTime(double time)
{
    cfront.max_time = time;
}

void FrontWrapper::setMaxStep(double step)
{
    cfront.max_step = step;
}

void FrontWrapper::setPrintTimeInterval(double interval)
{
    cfront.print_time_interval = interval;
}

void FrontWrapper::setMovieFrameInterval(double interval)
{
    cfront.movie_frame_interval = interval;
}

void FrontWrapper::setFrequencyOfRedistribution(double val, int i)
{
    cfront.Redist.freq_redist[i] = val;
}
// END TODO


//Getters

const int FrontWrapper::Dim() const
{
    return current_interface()->dim;
}

const double FrontWrapper::TimeStepFactor() const
{
    return cfront.Tstep.time_step_factor;
}


//Output

const bool FrontWrapper::IsSaveTime()
{
    if( FT_IsSaveTime(&cfront) )
        return true;
    else
        return false;
}

void FrontWrapper::Save()
{
    FT_Save(&cfront);
}

const bool FrontWrapper::IsDrawTime()
{
    if( FT_IsDrawTime(&cfront) )
        return true;
    else
        return false;
}

void FrontWrapper::Draw()
{
    FT_Draw(&cfront);
}

void FrontWrapper::PrintTimeStamp()
{
    FT_PrintTimeStamp(&cfront);
}


//Propagation

void FrontWrapper::Propagate()
{
    FT_Propagate(&cfront);
}

void FrontWrapper::SetNextTimeStep()
{
    FT_SetTimeStep(&cfront);
}

void FrontWrapper::SetOutputCounter()
{
    FT_SetOutputCounter(&cfront);
}

void FrontWrapper::TimeControlFilter()
{
    FT_TimeControlFilter(&cfront);
}

void FrontWrapper::AddTimeStepToCounter()
{
    FT_AddTimeStepToCounter(&cfront);
}

const bool FrontWrapper::TimeLimitReached()
{
    if( FT_TimeLimitReached(&cfront) )
        return true;
    else
        return false;
}


//Memory

void FrontWrapper::FreeMainIntfc()
{
    FT_FreeMainIntfc(&cfront);
}


//Mesh/Interface
        
void FrontWrapper::RedistMesh()
{
    FT_RedistMesh(&cfront);
}

