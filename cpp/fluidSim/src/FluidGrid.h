#pragma once
#include <immintrin.h>

#include <chrono>
#include <map>
#include <string>

#include "ThreadPool.h"

enum class Direction { NONE, HORIZONTAL, VERTICAL };

class FluidGrid {
public:
	explicit FluidGrid(int size);

	void	step(float dt, int iterations, double diffusionRate, double viscosity, double fadeRate);
	void	addDensity(int x, int y, float amount, float dt);
	void    addVelocity(int x, int y, float velocityXAmount, float velocityYAmount, float dt);
	std::vector<float> &getDensity();

private:
	std::vector<float> density, prevDensity;
	std::vector<float> velocityX, prevVelocityX;
	std::vector<float> velocityY, prevVelocityY;
	std::vector<float> tmp;

	std::map<std::string, double> timers;
	std::chrono::steady_clock::time_point timerT0;
	uint32_t runs = 0;

	std::chrono::steady_clock::time_point printT0;
	double  microsSinceLastPrint = 0;

	int		size;
	ThreadPool threadPool;

	// AVX2 (256-bit )for float (32-bit)
	static constexpr int SIMD_STRIDE = 8;

	// Optimal number of threads is 6 for quad-core with hyperthreading (6700k)
	// except for resolutions under 200
	static constexpr int WORKER_THREADS = 6;

	void	advect(Direction direction, std::vector<float> &arr, std::vector<float> &prevArr, float dt);
	void	advectLoop(int startIndex, int endIndex, std::vector<float> &arr, std::vector<float> &prevArr, float dt);
	void	project(int iterations, std::vector<float> &divergence, std::vector<float> &pressure);
	void	projectDivergenceLoop(int startIndex, int endIndex, std::vector<float> &divergence, std::vector<float> &pressure);
	void	projectGradientSubtractLoop(int startIndex, int endIndex, std::vector<float> &pressure);
	void	diffuse(Direction direction, int iterations, std::vector<float> &arr, std::vector<float> &prevArr, float dt, double diffusion);
	void	linearSolve(Direction direction, int iterations, std::vector<float> &arr, std::vector<float> &prevArr, float neighborDiffusion, float scaling);
	void	linearSolveLoop(int startIndex, int endIndex, std::vector<float> &arr, std::vector<float> &prevArr, float neighborDiffusion, __m256 _neighborDiffusion, float reciprocalScaling, __m256 _reciprocalScaling);
	void	setBounds(Direction direction, std::vector<float> &arr);
	void	fadeDensity(float dt, double fadeRate);

	void startTimer();
	void endTimer(const std::string &timerName);
	void checkPrint();
};


