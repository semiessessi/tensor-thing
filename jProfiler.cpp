/* jProfiler class by Semi Essessi
 *
 * (class description goes here)
 *
 */

#include "jProfiler.h"

#include <windows.h>

FILE* jProfiler::f = 0;
int jProfiler::instCount = 0;

jProfiler::jProfiler(const wchar_t* msg)
{
	// constructor code for menuVar
	int i = (int)wcslen(msg)+1;
	str = new wchar_t[i];
	memcpy(str, msg, sizeof(wchar_t)*i);
	str[i-1] =  0;

	QueryPerformanceFrequency(&tps);
	QueryPerformanceCounter(&start);

	instCount++;
	Initialise();
}

jProfiler::~jProfiler()
{
	// destructor code for menuVar
	QueryPerformanceCounter(&buf);
   	// work out change in time
   	double dt=((float)buf.QuadPart - (float)start.QuadPart)/(float)tps.QuadPart;
	fwprintf(f, L"%s : %.20f\r\n", str, dt);
	
	if(str) delete[] str;
	instCount--;
	Shutdown();
}

void jProfiler::Initialise()
{
	if(!f)
	{
		f = _wfopen(L"profilerlog.txt", L"wb");
		unsigned short a = 0xFEFF;
		fwrite(&a, sizeof(unsigned short), 1, f);
	}
}

void jProfiler::Shutdown()
{
	if(instCount==0) if(f) fclose(f);
}