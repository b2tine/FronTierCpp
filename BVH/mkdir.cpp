#include <sys/types.h>
#include <sys/stat.h>
#include <string>


static void mkdirTree(std::string sub, std::string dir)
{
    if(sub.length() == 0)
        return;

    int i = 0;
    for( i; i < sub.length(); i++)
    {
        dir += sub[i];
        if (sub[i] == '/')
            break;
    }

    mkdir(dir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

    if(i+1 < sub.length())
        mkdirTree(sub.substr(i+1), dir);
}


void create_directory(std::string new_dir)
{
    struct stat st;
    int status = stat(new_dir.c_str(), &st);
    if( status != 0 && !S_ISDIR(st.st_mode) )
        mkdirTree(new_dir, "");
    /*else
        std::cout << "Directory already exists.\n";*/
}
