#include "FluidGrid.h"
#include <iostream>
#include <algorithm>

FluidGrid::FluidGrid(int size) {
	density = (float *)calloc(size * size, sizeof(float));
	prevDensity = (float *)calloc(size * size, sizeof(float));

	velocityX = (float *)calloc(size * size, sizeof(float));
	prevVelocityX = (float *)calloc(size * size, sizeof(float));

	velocityY = (float *)calloc(size * size, sizeof(float));
	prevVelocityY = (float *)calloc(size * size, sizeof(float));

	this->size = size;

	if (!(density && prevDensity && velocityX && prevVelocityX && velocityY && prevVelocityY)) {
		throw std::runtime_error("Failed to allocate memory for FluidGrid\n");
	}
}

FluidGrid::~FluidGrid() {
	free(density);
	free(prevDensity);
	free(velocityX);
	free(prevVelocityX);
	free(velocityY);
	free(prevVelocityY);
}

void FluidGrid::step(float dt, int iterations, float diffusionRate, float viscosity, float fadeRate) {
	// Diffuse X velocity
	std::swap(this->velocityX, this->prevVelocityX);
	diffuse(Direction::HORIZONTAL, iterations, this->velocityX, this->prevVelocityX, dt, viscosity);

	// Diffuse Y velocity
	std::swap(this->velocityY, this->prevVelocityY);
	diffuse(Direction::VERTICAL, iterations, this->velocityY, this->prevVelocityY, dt, viscosity);

	// Conserve mass of the velocity field
	project(iterations, this->velocityX, this->velocityY, this->prevVelocityX, this->prevVelocityY);

	// Self advection
	std::swap(this->velocityX, this->prevVelocityX);
	std::swap(this->velocityY, this->prevVelocityY);
	advect(Direction::HORIZONTAL, this->velocityX, this->prevVelocityX, this->prevVelocityX, this->prevVelocityY, dt);
	advect(Direction::VERTICAL, this->velocityY, this->prevVelocityY, this->prevVelocityX, this->prevVelocityY, dt);

	// Conserve mass of the velocity field
	project(10, this->velocityX, this->velocityY, this->prevVelocityX, this->prevVelocityY);

	// Diffuse density
	std::swap(this->density, this->prevDensity);
	diffuse(Direction::NONE, iterations, this->density, this->prevDensity, dt, diffusionRate);

	// Advect density
	std::swap(this->density, this->prevDensity);
	advect(Direction::NONE, this->density, this->prevDensity, this->velocityX, this->velocityY, dt);

	// Fade the density to avoid filling the volume
	fadeDensity(dt, fadeRate);
}

void FluidGrid::addVelocity(int x, int y, float amountX, float amountY, float dt) {
		this->velocityX[x + y * this->size] += amountX * dt * this->size;
		this->velocityY[x + y * this->size] += amountY * dt * this->size;
}

void FluidGrid::addDensity(int x, int y, float amount, float dt) {
	this->density[x + y * this->size] += amount * dt * this->size;
}

// Linear backtrace
void FluidGrid::advect(Direction direction, float *arr, float *prevArr, float *velX, float *velyY, float dt) {
	for (int y = 1; y < size-1; y++) {
		for (int x = 1; x < size-1; x++) {
			float prevX = x - velX[x+y*size] * dt * this->size;
			float prevY = y - velyY[x+y*size] * dt * this->size;

			prevX = std::clamp(prevX, 0.5f, size - 1.5f);
			prevY = std::clamp(prevY, 0.5f, size - 1.5f);

			int prevXInt = (int)prevX;
			int prevYInt = (int)prevY;

			float prevXDecimals = prevX - prevXInt;
			float prevXDecimalsRemainder = 1.0f - prevXDecimals;

			float prevYDecimals = prevY - prevYInt;
			float prevYDecimalsRemainder = 1.0f - prevYDecimals;

			arr[x + y * this->size] =
				prevXDecimalsRemainder * (prevYDecimalsRemainder * prevArr[prevXInt	    + prevYInt * this->size] + prevYDecimals * prevArr[prevXInt	   + (prevYInt+1) * this->size]) +
				prevXDecimals		   * (prevYDecimalsRemainder * prevArr[(prevXInt+1) + prevYInt * this->size] + prevYDecimals * prevArr[(prevXInt+1) + (prevYInt+1) * this->size]);
		}
	}

	setBounds(direction, arr);
}

// Hodge decomposition
void FluidGrid::project(int iterations, float *velocityX, float *velocityY, float *p, float *div) {
	float sizeReciprocal = 1.0f / this->size;

	// Compute height field (Poisson-pressure equation)
	for (int y = 1; y < this->size - 1; y++) {
		for (int x = 1; x < this->size - 1; x++) {
			div[x + y * size] =
				-0.5f * (
					velocityX[(x-1) + y		* this->size] - velocityX[(x+1) + y		* this->size] +
					velocityY[x		+ (y-1) * this->size] - velocityY[x		+ (y+1) * this->size]
				) * sizeReciprocal;
			p[x + y * size] = 0;
		}
	}
	setBounds(Direction::NONE, div);
	setBounds(Direction::NONE, p);
	linearSolve(Direction::NONE, iterations, p, div, 1, 4);

	// Compute mass conserving field (Velocity field - Height field)
	for (int y = 1; y < this->size - 1; y++) {
		for (int x = 1; x < this->size - 1; x++) {
			velocityX[x + y * this->size] -= 0.5f * (p[(x-1) + y		* this->size] - p[(x+1) + y		* this->size]) * this->size;
			velocityY[x + y * this->size] -= 0.5f * (p[x		+ (y-1) * this->size] - p[x		+ (y+1) * this->size]) * this->size;
		}
	}

	setBounds(Direction::HORIZONTAL, velocityX);
	setBounds(Direction::VERTICAL, velocityY);
}

// Mass conserving
void FluidGrid::diffuse(Direction direction, int iterations, float *arr, float *prevArr, float dt, float diffusion) {
	float neighborDiffusion = diffusion * dt * this->size * this->size;
	float scaling = 1 + 4 * neighborDiffusion;

	linearSolve(direction, iterations, arr, prevArr, neighborDiffusion, scaling);
}

// Using Gauss-Seidel relaxation
void FluidGrid::linearSolve(Direction direction, int iterations, float *arr, float *prevArr, float neighborDiffusion, float scaling) {
	float reciprocalScaling = 1.0f / scaling;

	for (int iteration = 0; iteration < iterations; iteration++) {
		for (int y = 1; y < this->size - 1; y++) {
			for (int x = 1; x < this->size - 1; x++) {
				float previous = prevArr[x + y * this->size];
				float neighbors = 
					neighborDiffusion * (
						arr[x	  + ((y-1) * this->size)] +
						arr[x	  + ((y+1) * this->size)] +
						arr[(x-1) + (y	   * this->size)] +
						arr[(x+1) + (y	   * this->size)]
					);
				arr[x + y * this->size] = (previous + neighbors) * reciprocalScaling;
			}
		}
		// Controlling boundary after every iterations
		setBounds(direction, arr);
	}
}

void FluidGrid::setBounds(Direction direction, float *arr) {
	for (int i = 1; i < this->size-1; i++) {
		// Left and right edge get the reverse velocity of the neighbor if in X direction
		arr[0		 + i * size] = (direction == Direction::HORIZONTAL) ? -arr[1		  + i * size] : arr[1		 + i * size];
		arr[(size-1) + i * size] = (direction == Direction::HORIZONTAL) ? -arr[(size-2) + i * size] : arr[(size-2) + i * size];

		// Top and bottom edge get the reverse velocity of the neighbor if in Y direction
		arr[i + 0		 * size] = (direction == Direction::VERTICAL) ? -arr[i + 1		* size] : arr[i + 1		   * size];
		arr[i + (size-1) * size] = (direction == Direction::VERTICAL) ? -arr[i + (size-2) * size] : arr[i + (size-2) * size];
	}

	// The corners get the average values of their neighbors
	arr[0		 + 0		* size] = 0.5f * (arr[1		  + 0		 * size] + arr[0		+ 1		   * size]);
	arr[(size-1) + 0		* size] = 0.5f * (arr[(size-2) + 0		 * size] + arr[(size-1) + 1		   * size]);
	arr[0		 + (size-1)	* size] = 0.5f * (arr[0		  + (size-2) * size] + arr[1		+ (size-1) * size]);
	arr[(size-1) + (size-1) * size] = 0.5f * (arr[(size-2) + (size-1) * size] + arr[(size-1) + (size-2) * size]);
}

void FluidGrid::fadeDensity(float dt, float fadeRate) {
	for (int y = 1; y < this->size - 1; y++) {
		for (int x = 1; x < this->size - 1; x++) {
			this->density[x + y * this->size] *= 1 - dt * fadeRate * this->size;
		}
	}
}

float *FluidGrid::getDensity() {
	return this->density;
}