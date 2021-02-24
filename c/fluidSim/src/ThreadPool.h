#pragma once
#include <mutex>
#include <condition_variable>
#include <functional>
#include <thread>
#include <iostream>
#include <vector>

class Worker {
	bool shouldExit = false;
	bool hasWork = false;
	std::condition_variable sleepCV;
	std::mutex sleepMux;
	std::function<void()> work;

public:
	std::thread thread;

	void threadLoop();
	void start(std::function<void()> work);
	void terminate();
};

class ThreadPool {
	Worker *workers = nullptr; // Has to be array since mutexes are not allowed to be copied or moved
	int numWorkers = 0;

public:
	~ThreadPool();
	void init(int numWorkers);
	void computeOnThreads(int size, std::function<void(int, int)> work);
};