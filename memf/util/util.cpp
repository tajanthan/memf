
#ifdef _WIN32
#include "windows.h"
#include "psapi.h"	// order important
#elif __unix__
#include <unistd.h>
#include <sys/resource.h>
//#include <cstdio>
//#include <cstdlib>
//#include <cstring>
//#include <fstream>
#endif

#include "util.h"

using namespace std;

#ifdef _WIN32

size_t getMemoryUsage()
{
	PROCESS_MEMORY_COUNTERS pmc;
	GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(pmc));
	return (size_t)pmc.PeakPagefileUsage;
}

#elif __unix__

//// memory usage
//long long int parseLine(char* line)
//{
//    int i = strlen(line);
//    while (*line < '0' || *line > '9') line++;
//    line[i-3] = '\0';
//    //i = atoi(line);
//    return atoll(line);
//}

//int getVmValue()
//{ //Note: this value is in KB!
//    FILE* file = fopen("/proc/self/status", "r");
//    int result = -1;
//    char line[128];
//
//
//    while (fgets(line, 128, file) != NULL){
//        if (strncmp(line, "VmSize:", 7) == 0){
//            result = parseLine(line);
//            break;
//        }
//    }
//    fclose(file);
//    return result;
//}

//int getRamValue()
size_t getMemoryUsage()
{ //Note: this value is in KB!
//    FILE* file = fopen("/proc/self/status", "r");
//    size_t result = 0;
//    char line[128];
//
//
//    while (fgets(line, 128, file) != NULL){
//        if (strncmp(line, "VmRSS:", 6) == 0){
//            result = parseLine(line);
//            break;
//        }
//    }
//    fclose(file);
//    return result * 1024;

	struct rusage rusage;
	getrusage( RUSAGE_SELF, &rusage );
	return (size_t)(rusage.ru_maxrss * 1024L);
}
//
#endif

