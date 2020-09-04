#include "os_service.h"
#include <stdio.h>
#include <unistd.h>
#include <sys/types.h>
#include <signal.h>
#include <sstream>
#include <string>
#include <iostream>

namespace os
{
void os_service::print_callstack()
{
    std::string cmd = "pstack";
    pid_t pid = getpid();
    std::stringstream ss;
    ss << cmd;
    ss << " ";
    ss << pid;
    ss << "| tee sig.rpt";
    bool res = system( ss.str().c_str() );
    if ( res == true )
    {
        std::cout << "Fail to dump callstack";
    }
    /*std::string cmd2 = "sed \'s/^/!/\' sig.rpt -i";
    system( cmd2.c_str() );*/
}

void os_service::initialize_handle()
{
    struct sigaction act;
    act.sa_sigaction = _handler;
    sigemptyset( &act.sa_mask );
    act.sa_flags = SA_SIGINFO;
    sigaction( SIGINT, &act, 0 );
    sigaction( SIGABRT, &act, 0 );
    sigaction( SIGFPE, &act, 0 );
    sigaction( SIGSEGV, &act, 0 );
}

void os_service::_handler( int, siginfo_t*, void* )
{
    print_callstack();
    exit( 0 );
}

}
