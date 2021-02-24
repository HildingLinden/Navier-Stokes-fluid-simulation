#pragma once
#include <chrono>
#include <map>
#include <string>
#include <immintrin.h>
#include "ThreadPool.h"

enum class Direction { NONE, HORIZONTAL, VERTICAL };

class FluidGrid {
public:
	~FluidGrid();
	FluidGrid(int size);

	void	step(float dt, int iterations, double diffusionRate, double viscosity, double fadeRate);
	void	addDensity(int x, int y, float amount, float dt);
	void    addVelocity(int x, int y, float velocityXAmount, float velocityYAmount, float dt);
	float *getDensity();

private:
	float *density, *prevDensity;
	float *velocityX, *prevVelocityX;
	float *velocityY, *prevVelocityY;
	float *tmp;

	std::map<std::string, double> timers;
	std::chrono::steady_clock::time_point timerT0;
	long	runs = 0;

	std::chrono::steady_clock::time_point printT0;
	double  microsSinceLastPrint = 0;

	int		size;
	ThreadPool threadPool;

	void	advect(Direction direction, float *arr, float *prevArr, float *velocityX, float *velocityY, float dt);
	void	project(int iterations, float *velocityX, float *velocityY, float *p, float *div);
	void	projectHeightMapLoop(int startIndex, int endIndex, float *velX, float *velY, float *p, float *div);
	void	projectMassConservLoop(int startIndex, int endIndex, float *p, float *div);
	void	diffuse(Direction direction, int iterations, float *arr, float *prevArr, float dt, double diffusion);
	void	linearSolve(Direction direction, int iterations, float *arr, float *prevArr, float neighborDiffusion, float scaling);
	void	linearSolveLoop(int startIndex, int endIndex, float *arr, float *prevArr, float neighborDiffusion, __m256 _neighborDiffusion, float reciprocalScaling, __m256 _reciprocalScaling);
	void	setBounds(Direction direction, float *arr);
	void	fadeDensity(float dt, double fadeRate);

	void startTimer();
	void endTimer(std::string timerName);
	void checkPrint();
};


