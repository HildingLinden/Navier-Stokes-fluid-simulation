#pragma once
#include <chrono>
#include <map>
#include <string>
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

	void	advect(Direction direction, float *arr, float *prevArr, float *velocityX, float *velocityY, float dt);
	void	project(int iterations, float *velocityX, float *velocityY, float *p, float *div);
	void	diffuse(Direction direction, int iterations, float *arr, float *prevArr, float dt, double diffusion);
	void	linearSolve(Direction direction, int iterations, float *arr, float *prevArr, float neighborDiffusion, float scaling);
	void	setBounds(Direction direction, float *arr);
	void	fadeDensity(float dt, double fadeRate);

	void startTimer();
	void endTimer(std::string timerName);
	void checkPrint();
};


