#pragma once
//timers are in milliseconds
struct Timer
{
private:
	std::chrono::time_point<std::chrono::high_resolution_clock> start_frame_time;
public:
	void Begin()
	{
		start_frame_time = std::chrono::high_resolution_clock::now();
	}

	bool TimeElapsed(int duration_milliseconds)
	{
		auto endTimepoint = std::chrono::high_resolution_clock::now();

		auto start = std::chrono::time_point_cast<std::chrono::milliseconds>(start_frame_time).time_since_epoch().count();
		auto end = std::chrono::time_point_cast<std::chrono::milliseconds>(endTimepoint).time_since_epoch().count();
		auto duration = end - start;

		if (duration >= duration_milliseconds)
			return true;
		else
			return false;
	}

	float getTime()
	{
		auto endTimepoint = std::chrono::high_resolution_clock::now();

		auto start = std::chrono::time_point_cast<std::chrono::milliseconds>(start_frame_time).time_since_epoch().count();
		auto end = std::chrono::time_point_cast<std::chrono::milliseconds>(endTimepoint).time_since_epoch().count();
		auto duration = end - start;
		return duration;
	}

	float getTimeMicro()
	{
		auto endTimepoint = std::chrono::high_resolution_clock::now();

		auto start = std::chrono::time_point_cast<std::chrono::microseconds>(start_frame_time).time_since_epoch().count();
		auto end = std::chrono::time_point_cast<std::chrono::microseconds>(endTimepoint).time_since_epoch().count();
		auto duration = end - start;
		return duration;
	}

	float getTimeNano()
	{
		auto endTimepoint = std::chrono::high_resolution_clock::now();

		auto start = std::chrono::time_point_cast<std::chrono::nanoseconds>(start_frame_time).time_since_epoch().count();
		auto end = std::chrono::time_point_cast<std::chrono::nanoseconds>(endTimepoint).time_since_epoch().count();
		auto duration = end - start;
		return duration;
	}
};