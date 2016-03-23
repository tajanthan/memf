#pragma once

#include <cassert>
#include <ctime>
#include <fstream>
#include <iostream>
#include <string>

class Timer
{
public:

	Timer()
	{
	}

	~Timer()
	{
	}

	void startTimer()
	{
		startTime = clock();
		incrTime = startTime;
		totalTime = 0;
	}

	void printData(std::ofstream& fout, const std::string& str)
	{
		assert(startTime);
		clock_t endTime = clock();
		double iTime = ((double)endTime - (double)incrTime) / CLOCKS_PER_SEC;
		totalTime = ((double)endTime - (double)startTime) / CLOCKS_PER_SEC;
		std::cout << "Time: [" << iTime << "s, " << totalTime << "s]\t" << str << std::endl;
		fout << "Time: [" << iTime << "s, " << totalTime << "s]\t" << str << std::endl;
		incrTime = endTime;
	}

	void pause()
	{
		assert(incrTime);
		clock_t endTime = clock();
		totalTime += ((double)endTime - (double)incrTime) / CLOCKS_PER_SEC;
	}

	void resume()
	{
		incrTime = clock();
	}

	double getTotalTime()
	{
		return totalTime;
	}

	void stopTimer()
	{
		startTime = 0;
		incrTime = 0;
		totalTime = 0;
	}

private:
	clock_t startTime, incrTime;
	double totalTime;
};

