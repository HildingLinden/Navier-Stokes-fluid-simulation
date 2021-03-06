#include "FluidGrid.h"

#include <algorithm>
#include <chrono>
#include <iostream>


FluidGrid::FluidGrid(const int size) {
	printT0 = std::chrono::steady_clock::now();

	density.resize(size * size);
	prevDensity.resize(size * size);

	velocityX.resize(size * size);
	prevVelocityX.resize(size * size);

	velocityY.resize(size * size);
	prevVelocityY.resize(size * size);

	tmp.resize(size * size);

	this->size = size;

	threadPool.init(WORKER_THREADS);
}

void FluidGrid::step(float dt, int iterations, double diffusionRate, double viscosity, double fadeRate) {
	startTimer();
	// Diffuse X velocity
	diffuse(Direction::HORIZONTAL, iterations, velocityX, prevVelocityX, dt, viscosity);

	// Diffuse Y velocity
	diffuse(Direction::VERTICAL, iterations, velocityY, prevVelocityY, dt, viscosity);
	endTimer("Diffuse");

	startTimer();
	// Conserve mass of the velocity field
	project(iterations, prevVelocityX, prevVelocityY);
	endTimer("Project");

	startTimer();
	// Self advection
	std::swap(velocityX, prevVelocityX);
	std::swap(velocityY, prevVelocityY);
	advect(Direction::HORIZONTAL, velocityX, prevVelocityX, dt);
	advect(Direction::VERTICAL, velocityY, prevVelocityY, dt);
	endTimer("Advect");

	startTimer();
	// Conserve mass of the velocity field
	project(iterations, prevVelocityX, prevVelocityY);
	endTimer("Project");

	startTimer();
	// Diffuse density
	std::swap(density, prevDensity);
	diffuse(Direction::NONE, iterations, density, prevDensity, dt, diffusionRate);
	endTimer("Diffuse");

	startTimer();
	// Advect density
	std::swap(velocityX, prevVelocityX);
	std::swap(velocityY, prevVelocityY);
	std::swap(density, prevDensity);
	advect(Direction::NONE, density, prevDensity, dt);
	endTimer("Advect");

	// Fade the density to avoid filling the volume
	fadeDensity(dt, fadeRate);

	checkPrint();
}

void FluidGrid::addVelocity(int x, int y, float velocityXAmount, float velocityYAmount, float dt) {
	prevVelocityX[x + y * size] += velocityXAmount * dt * static_cast<float>(size);
	prevVelocityY[x + y * size] += velocityYAmount * dt * static_cast<float>(size);
}

void FluidGrid::addDensity(int x, int y, float amount, float dt) {
	density[x + y * size] += amount * dt * static_cast<float>(size);
}

// Linear backtrace
void FluidGrid::advect(Direction direction, std::vector<float> &arr, std::vector<float> &prevArr, float dt) {
	threadPool.computeOnThreads(
		size,
		[&](int startIndex, int endIndex) {
			advectLoop(startIndex, endIndex, arr, prevArr, dt);
		}
	);

	setBounds(direction, arr);
}

void FluidGrid::advectLoop(int startIndex, int endIndex, std::vector<float> &arr, std::vector<float> &prevArr, float dt) {
	const float scaling = dt * static_cast<float>(size);
	const float maxClamp = static_cast<float>(size) - 1.5F;
	constexpr float minClamp = 0.5F;

	const __m256 _xIncrements = _mm256_set_ps(7.0F, 6.0F, 5.0F, 4.0F, 3.0F, 2.0F, 1.0F, 0.0F);
	const __m256 _scaling = _mm256_set1_ps(scaling);
	const __m256 _maxClamp = _mm256_set1_ps(maxClamp);
	const __m256 _minClamp = _mm256_set1_ps(0.5F);
	const __m256 _one = _mm256_set1_ps(1.0F);
	const __m256i _size = _mm256_set1_epi32(size);

	for (int y = startIndex; y < endIndex; y++) {
		int x = 1;

		const __m256 _y = _mm256_set1_ps(static_cast<float>(y));

		for (; x < size - (SIMD_STRIDE + 1); x += SIMD_STRIDE) {
			// Increment x by 0 to 7 since we are computing 8 iterations of x at once
			__m256 _x = _mm256_set1_ps(static_cast<float>(x));
			_x = _mm256_add_ps(_x, _xIncrements);

			// Subtract the current scaled velocity from y and incremented x
			__m256 _prevX = _mm256_load_ps(&prevVelocityX[x + y * size]);
			__m256 _prevY = _mm256_load_ps(&prevVelocityY[x + y * size]);

			_prevX = _mm256_mul_ps(_prevX, _scaling);
			_prevY = _mm256_mul_ps(_prevY, _scaling);
			_prevX = _mm256_sub_ps(_x, _prevX);
			_prevY = _mm256_sub_ps(_y, _prevY);

			// Clamp indices to minClamp and maxClamp
			_prevX = _mm256_max_ps(_prevX, _minClamp);
			_prevX = _mm256_min_ps(_prevX, _maxClamp);
			_prevY = _mm256_max_ps(_prevY, _minClamp);
			_prevY = _mm256_min_ps(_prevY, _maxClamp);

			// Get the integer part and decimal part of the indices
			// as well as 1 - decimal part
			const __m256 _prevXRound = _mm256_floor_ps(_prevX);
			const __m256 _prevYRound = _mm256_floor_ps(_prevY);

			const __m256i _prevXInt = _mm256_cvtps_epi32(_prevXRound);
			const __m256i _prevYInt = _mm256_cvtps_epi32(_prevYRound);

			const __m256 _prevXDecimals = _mm256_sub_ps(_prevX, _prevXRound);
			const __m256 _prevYDecimals = _mm256_sub_ps(_prevY, _prevYRound);

			const __m256 _prevXRemainingDecimals = _mm256_sub_ps(_one, _prevXDecimals);
			const __m256 _prevYRemainingDecimals = _mm256_sub_ps(_one, _prevYDecimals);

			// Multiply y by size and add x to get the actual indices
			__m256i _prevIndices = _mm256_mullo_epi32(_prevYInt, _size);
			_prevIndices = _mm256_add_epi32(_prevIndices, _prevXInt);

			// Get the middle, right, down and downRight cells by changing the base address of the gather
			const __m256 _middle = _mm256_i32gather_ps(&prevArr[0], _prevIndices, sizeof(float));
			const __m256 _right = _mm256_i32gather_ps(&prevArr[1], _prevIndices, sizeof(float));
			const __m256 _down = _mm256_i32gather_ps(&prevArr[size], _prevIndices, sizeof(float));
			const __m256 _downRight = _mm256_i32gather_ps(&prevArr[size + 1], _prevIndices, sizeof(float));

			// Do the linear interpolation (See below for explanation)
			const __m256 _middleScaling = _mm256_mul_ps(_prevXRemainingDecimals, _prevYRemainingDecimals);
			const __m256 _downScaling = _mm256_mul_ps(_prevXRemainingDecimals, _prevYDecimals);
			const __m256 _rightScaling = _mm256_mul_ps(_prevXDecimals, _prevYRemainingDecimals);
			const __m256 _downRightScaling = _mm256_mul_ps(_prevXDecimals, _prevYDecimals);

			__m256 _result = _mm256_mul_ps(_middle, _middleScaling);
			_result = _mm256_fmadd_ps(_down, _downScaling, _result);
			_result = _mm256_fmadd_ps(_right, _rightScaling, _result);
			_result = _mm256_fmadd_ps(_downRight, _downRightScaling, _result);

			_mm256_storeu_ps(&arr[x + y * size], _result);
		}

		for (; x < size - 1; x++) {
			// Subtract velocity from current cell to get previous cell indices
			float prevX = static_cast<float>(x) - prevVelocityX[x + y * size] * scaling;
			float prevY = static_cast<float>(y) - prevVelocityY[x + y * size] * scaling;

			// Clamp previous cells in case they are outside of the domain
			prevX = std::clamp(prevX, minClamp, maxClamp);
			prevY = std::clamp(prevY, minClamp, maxClamp);

			/*
			Casting the index of the previous cell to integers truncates the index which means that
			the index is moved up and to the left
			*/
			const int prevXInt = static_cast<int>(prevX);
			const int prevYInt = static_cast<int>(prevY);

			/*
			To interpolate we then:
				multiply the decimals of x and y with the cell below and to the right
				multiply the decimals of x and remaining decimals of y with the cell to the right
				multiply the decimals of y and remaining decimals of x with the cell above
				multiply the remaining decimals of x and y with the current cell
			and add them all together
			*/
			float prevXDecimals = prevX - static_cast<float>(prevXInt);
			float prevXRemainingDecimals = 1.0F - prevXDecimals;

			float prevYDecimals = prevY - static_cast<float>(prevYInt);
			float prevYRemainingDecimals = 1.0F - prevYDecimals;

			arr[x + y * size] =
				prevXRemainingDecimals * (prevYRemainingDecimals * prevArr[prevXInt + prevYInt * size] + prevYDecimals * prevArr[prevXInt + (prevYInt + 1) * size]) +
				prevXDecimals * (prevYRemainingDecimals * prevArr[(prevXInt + 1) + prevYInt * size] + prevYDecimals * prevArr[(prevXInt + 1) + (prevYInt + 1) * size]);
		}
	}
}

// Hodge decomposition
void FluidGrid::project(int iterations, std::vector<float> &divergence, std::vector<float> &pressure) {
	/* 
	Compute height field (Poisson-pressure equation on threads
	*/
	threadPool.computeOnThreads(size, [&](int startIndex, int endIndex) {projectDivergenceLoop(startIndex, endIndex, divergence, pressure); });

	setBounds(Direction::NONE, divergence);
	setBounds(Direction::NONE, pressure);
	linearSolve(Direction::NONE, iterations, pressure, divergence, 1, 4);

	/*
	Compute mass conserving field (Velocity field - Height field) on threads
	*/
	threadPool.computeOnThreads(size, [&](int startIndex, int endIndex) {projectGradientSubtractLoop(startIndex, endIndex, pressure); });

	setBounds(Direction::HORIZONTAL, velocityX);
	setBounds(Direction::VERTICAL, velocityY);
}

void FluidGrid::projectDivergenceLoop(int startIndex, int endIndex, std::vector<float> &divergence, std::vector<float> &pressure) {
	float sizeReciprocal = 1.0F / static_cast<float>(size);

	const __m256 _sizeReciprocal = _mm256_set1_ps(sizeReciprocal);
	const __m256 _minusHalf = _mm256_set1_ps(-0.5F);
	const __m256 _zero = _mm256_setzero_ps();

	for (int y = startIndex; y < endIndex; y++) {
		int x = 1;

		for (; x < size - (SIMD_STRIDE + 1); x += SIMD_STRIDE) {
			/*
			x = x[left] - x[right]
			y = y[up] - y[down]
			cell = x + y
			cell *= -0.5
			cell *= sizeReciprocal
			*/
			const __m256 _up = _mm256_loadu_ps(&velocityY[x + (y - 1) * size]);
			const __m256 _down = _mm256_loadu_ps(&velocityY[x + (y + 1) * size]);
			const __m256 _left = _mm256_loadu_ps(&velocityX[(x - 1) + y * size]);
			const __m256 _right = _mm256_loadu_ps(&velocityX[(x + 1) + y * size]);

			const __m256 _x = _mm256_sub_ps(_left, _right);
			const __m256 _y = _mm256_sub_ps(_up, _down);
			__m256 _cell = _mm256_add_ps(_x, _y);
			_cell = _mm256_mul_ps(_cell, _minusHalf);
			_cell = _mm256_mul_ps(_cell, _sizeReciprocal);
			_mm256_storeu_ps(&divergence[x + y * size], _cell);
			_mm256_storeu_ps(&pressure[x + y * size], _zero);
		}

		constexpr float neighborScaling = 0.5F;

		// Take care of rest of the cells if x-2 is not evenly divisible by 8
		for (; x < size - 1; x++) {
			divergence[x + y * size] =
				-neighborScaling * (
					velocityX[(x - 1) + y * size] - velocityX[(x + 1) + y * size] +
					velocityY[x + (y - 1) * size] - velocityY[x + (y + 1) * size])
				* sizeReciprocal;
			pressure[x + y * size] = 0;
		}
	}
}

void FluidGrid::projectGradientSubtractLoop(int startIndex, int endIndex, std::vector<float> &pressure) {
	const __m256 _size = _mm256_set1_ps(static_cast<float>(size));
	const __m256 _minusHalf = _mm256_set1_ps(-0.5F);

	for (int y = startIndex; y < endIndex; y++) {
		int x = 1;

		for (; x < size - (SIMD_STRIDE + 1); x += SIMD_STRIDE) {
			/*
			val = p[left] - p[right];
			val *= size;
			oldVal += -0.5F * val;
			*/
			const __m256 _left = _mm256_loadu_ps(&pressure[(x - 1) + y * size]);
			const __m256 _right = _mm256_loadu_ps(&pressure[(x + 1) + y * size]);
			__m256 _val = _mm256_sub_ps(_left, _right);
			_val = _mm256_mul_ps(_val, _size);
			__m256 _oldVal = _mm256_loadu_ps(&velocityX[x + y * size]);
			_oldVal = _mm256_fmadd_ps(_val, _minusHalf, _oldVal);
			_mm256_storeu_ps(&velocityX[x + y * size], _oldVal);

			/*
			val = p[up] - p[down];
			val *= size;
			oldVal += -0.5F * val;
			*/
			const __m256 _up = _mm256_loadu_ps(&pressure[x + (y - 1) * size]);
			const __m256 _down = _mm256_loadu_ps(&pressure[x + (y + 1) * size]);
			_val = _mm256_sub_ps(_up, _down);
			_val = _mm256_mul_ps(_val, _size);
			_oldVal = _mm256_loadu_ps(&velocityY[x + y * size]);
			_oldVal = _mm256_fmadd_ps(_val, _minusHalf, _oldVal);
			_mm256_storeu_ps(&velocityY[x + y * size], _oldVal);
		}

		constexpr float neighborScaling = 0.5F;

		// Take care of rest of the cells if x-2 is not evenly divisible by SIMD_STRIDE
		for (; x < size - 1; x++) {
			velocityX[x + y * size] -= neighborScaling * (pressure[(x - 1) + y * size] - pressure[(x + 1) + y * size]) * static_cast<float>(size);
			velocityY[x + y * size] -= neighborScaling * (pressure[x + (y - 1) * size] - pressure[x + (y + 1) * size]) * static_cast<float>(size);
		}
	}
}

// Mass conserving
void FluidGrid::diffuse(Direction direction, int iterations, std::vector<float> &arr, std::vector<float> &prevArr, float dt, double diffusion) {
	const float neighborDiffusion = static_cast<float>(diffusion * dt * size * size);
	const float scaling = 1 + 4 * neighborDiffusion;

	linearSolve(direction, iterations, arr, prevArr, neighborDiffusion, scaling);
}

// Using Jacobi relaxation
void FluidGrid::linearSolve(Direction direction, int iterations, std::vector<float> &arr, std::vector<float> &prevArr, float neighborDiffusion, float scaling) {
	const float reciprocalScaling = 1.0F / scaling;
	const __m256 _reciprocalScaling = _mm256_set1_ps(reciprocalScaling);
	const __m256 _neighborDiffusion = _mm256_set1_ps(neighborDiffusion);

	for (int iteration = 0; iteration < iterations; iteration++) {
		/*
		Run on threads
		*/
		threadPool.computeOnThreads(size, [&](int startIndex, int endIndex) {
				linearSolveLoop(startIndex, endIndex, arr, prevArr, neighborDiffusion,  _neighborDiffusion, reciprocalScaling, _reciprocalScaling);
		});

		std::swap(tmp, arr);

		// Controlling boundary after every iterations
		setBounds(direction, arr);
	}
}

void FluidGrid::linearSolveLoop(int startIndex, int endIndex, std::vector<float> &arr, std::vector<float> &prevArr, float neighborDiffusion, __m256 _neighborDiffusion, float reciprocalScaling, __m256 _reciprocalScaling) {
	for (int y = startIndex; y < endIndex; y++) {
		int x = 1;

		for (; x < size -  (SIMD_STRIDE + 1); x += SIMD_STRIDE) {
			const __m256 _up = _mm256_loadu_ps(&arr[x + (y - 1) * size]);
			const __m256 _left = _mm256_loadu_ps(&arr[(x - 1) + y * size]);
			const __m256 _right = _mm256_loadu_ps(&arr[(x + 1) + y * size]);
			const __m256 _down = _mm256_loadu_ps(&arr[x + (y + 1) * size]);

			/*
			cell = previous;
			cell +=	(up * neighborDiffusion);
			cell += (left * neighborDiffusion);
			cell += (right * neighborDiffusion);
			cell += (down * neighborDiffusion);
			cell *= reciprocalScaling
			*/
			__m256 _cell = _mm256_loadu_ps(&prevArr[x + y * size]);
			_cell = _mm256_fmadd_ps(_up, _neighborDiffusion, _cell);
			_cell = _mm256_fmadd_ps(_left, _neighborDiffusion, _cell);
			_cell = _mm256_fmadd_ps(_right, _neighborDiffusion, _cell);
			_cell = _mm256_fmadd_ps(_down, _neighborDiffusion, _cell);
			_cell = _mm256_mul_ps(_cell, _reciprocalScaling);

			// Store result
			_mm256_storeu_ps(&tmp[x + y * size], _cell);
		}

		// Take care of rest of the cells if x-2 is not evenly divisible by SIMD_STRIDE
		for (; x < size - 1; x++) {
			float neighbors =
				neighborDiffusion * (
					arr[x + ((y - 1) * size)] +
					arr[x + ((y + 1) * size)] +
					arr[(x - 1) + (y * size)] +
					arr[(x + 1) + (y * size)]);
			float previous = prevArr[x + y * size];
			tmp[x + y * size] = (previous + neighbors) * reciprocalScaling;
		}
	}
}

void FluidGrid::setBounds(Direction direction, std::vector<float> &arr) {
	for (int i = 1; i < size-1; i++) {
		// Left and right edge get the reverse velocity of the neighbor if in X direction
		arr[0		 + i * size] = (direction == Direction::HORIZONTAL) ? -arr[1		  + i * size] : arr[1		 + i * size];
		arr[(size-1) + i * size] = (direction == Direction::HORIZONTAL) ? -arr[(size-2) + i * size] : arr[(size-2) + i * size];

		// Top and bottom edge get the reverse velocity of the neighbor if in Y direction
		arr[i + 0		 * size] = (direction == Direction::VERTICAL) ? -arr[i + 1		* size] : arr[i + 1		   * size];
		arr[i + (size-1) * size] = (direction == Direction::VERTICAL) ? -arr[i + (size-2) * size] : arr[i + (size-2) * size];
	}

	constexpr float neighborScaling = 0.5F;

	// The corners get the average values of their neighbors
	arr[0		 + 0		* size] = neighborScaling * (arr[1		  + 0		 * size] + arr[0		+ 1		   * size]);
	arr[(size-1) + 0		* size] = neighborScaling * (arr[(size-2) + 0		 * size] + arr[(size-1) + 1		   * size]);
	arr[0		 + (size-1)	* size] = neighborScaling * (arr[0		  + (size-2) * size] + arr[1		+ (size-1) * size]);
	arr[(size-1) + (size-1) * size] = neighborScaling * (arr[(size-2) + (size-1) * size] + arr[(size-1) + (size-2) * size]);
}

void FluidGrid::fadeDensity(float dt, double fadeRate) {
	const double scaledFadeRate = 1 - dt * fadeRate * size;
	for (int y = 1; y < size - 1; y++) {
		for (int x = 1; x < size - 1; x++) {
			density[x + y * size] = static_cast<float>(scaledFadeRate * density[x + y * size]);
		}
	}
}

void FluidGrid::startTimer() {
	timerT0 = std::chrono::steady_clock::now();
}

void FluidGrid::endTimer(const std::string &timerName) {
	const auto t1 = std::chrono::steady_clock::now();
	const std::chrono::duration<double, std::micro> microSeconds = t1 - timerT0;

	const auto iterator = timers.find(timerName);
	if (iterator != timers.end()) {
		iterator->second += microSeconds.count();
	} else {
		timers.insert(std::pair<std::string, double>(timerName, microSeconds.count()));
	}
}

void FluidGrid::checkPrint() {
	runs++;
	const auto printT1 = std::chrono::steady_clock::now();
	const std::chrono::duration<double, std::micro> printMicroSeconds = printT1 - printT0;
	microsSinceLastPrint += printMicroSeconds.count();
	constexpr double ONE_SECOND_IN_MICROSECONDS = 1000000.0;
	if (microsSinceLastPrint > ONE_SECOND_IN_MICROSECONDS) {
		for (std::pair<const std::string, double> &timer : timers) {
			std::cout << timer.first << ": " << timer.second / runs << " microseconds\n";
			timer.second = 0;
		}
		std::cout << std::endl;
		microsSinceLastPrint = 0;
		runs = 0;
	}

	printT0 = printT1;
}

std::vector<float> &FluidGrid::getDensity() {
	return density;
}
