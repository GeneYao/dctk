#ifndef OS_SERVICE_H
#define OS_SERVICE_H
#include <signal.h>
namespace os
{

class os_service
{
public:
    static void print_callstack();
    static void initialize_handle();
private:
    static void _handler( int, siginfo_t*, void* );
};
}
#endif // SYSTEM_SERVICE_H
