/* jProfiler class by Semi Essessi
 *
 * (class description goes here)
 *
 */

#ifndef __JPROFILER_H
#define __JPROFILER_H

#include <stdio.h>
#include <windows.h>

class jProfiler
{
private:
	wchar_t* 		str;
	LARGE_INTEGER	start;
	LARGE_INTEGER	tps;
	LARGE_INTEGER	buf;

	static FILE*	f;
	static int		instCount;

	static void		Initialise();
	static void		Shutdown();
public:
	jProfiler(const wchar_t* msg);
	~jProfiler();
};

#endif