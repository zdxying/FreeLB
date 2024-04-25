#include "freelb.h"
#include "freelb.hh"

int main()
{
    DirCreator::Create_Dir("./test/error"); // Error Creating Directory: ./test/error
    DirCreator::Create_Dir("./test0/", "./test1/", "./test2/");
    DirCreator::Create_Dir("./test/", "./test/sub/", "./test/sub/sub0/","./test/sub1/");
    // DirCreator_::Create_Dir("test");
    // DirCreator_::Create_Dir("test/sub");
    // DirCreator_::Create_Dir("test0", "test1", "test2");
    return 0;
}