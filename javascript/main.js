"use strict";

function main() {
	let canvas = document.getElementById("canvas");
	let context = canvas.getContext("2d");
	context.imageSmoothingEnabled = false;

	let physicsTimeText = document.getElementById("physicsTime");
	let drawTimeText = document.getElementById("drawTime");
	let timeStepText = document.getElementById("timeStep");

	let resolution = 100;
	let SCALE = canvas.width/resolution;

	let diffusion = 0.001;
	let viscosity = 0.00001;
	let fadeRate = 0.001;
	let iterations = 4;
	let density = 200;
	let velocity = 500;
	let timeScale = 1.0;

	let grid = new FluidGrid(resolution, diffusion, viscosity, fadeRate, iterations);

	let velX = 0;
	let velY = 0;

	// Create inital frame and set alpha to 255
	let frame = context.createImageData(resolution, resolution);
	for (let y = 0; y < resolution; y++) {
		for (let x = 0; x < resolution; x++) {
			frame.data[x*4 + y*4*resolution + 3] = 255;
		}
	}

	let timeStep0 = performance.now()

	function loop() {
		let timeStep1 = performance.now();
		let dt = (timeStep1 - timeStep0) / (1000 / timeScale);
		timeStep0 = timeStep1;

		timeStepText.innerHTML = "Time between frames " + (dt*(1000 / timeScale)).toFixed(3) + "ms";

		// Add velocity to the middle of the canvas in the direction of the mouse
		grid.addVelocity(resolution/2, resolution/2, velX, velY, dt);

		// Add density to the middle of the canvas
		grid.addDensity(resolution/2, resolution/2, density, dt);

		// Step the physics forward a time step
		let stepTime = grid.step(dt);

		physicsTimeText.innerHTML = "Physics time " + stepTime.toFixed(3) + "ms";

		let t0 = performance.now();

		// Reset canvas
		context.clearRect(0, 0, canvas.width, canvas.height);

		// Convert density into grayscale pixels in a frame
		for (let y = 0; y < resolution; y++) {
			for (let x = 0; x < resolution; x++) {
				let color = 255 * grid.density[x+y*resolution];
				frame.data[x*4 + y*4*resolution] = color;
				frame.data[x*4 + y*4*resolution + 1] = color;
				frame.data[x*4 + y*4*resolution + 2] = color;
			}
		}

		// Render frame on canvas then scale it up to the correct size
		context.putImageData(frame, 0, 0);
		context.drawImage(canvas, 0, 0, resolution, resolution, 0, 0, canvas.width, canvas.height);

		let t1 = performance.now();
		let drawTime = t1-t0;
		drawTimeText.innerHTML = "Draw time " + drawTime.toFixed(3) + "ms";

		window.requestAnimationFrame(loop);
	}

	window.requestAnimationFrame(loop);

	canvas.addEventListener("mousemove", (event) => {
		let x = event.clientX - canvas.width/2 - 10;
		let y = event.clientY - canvas.height/2 - 10;
		let angle = Math.atan2(y, x);
		velX = velocity * Math.cos(angle);
		velY = velocity * Math.sin(angle);
	});

	let intervals = new Array(0, 0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1);
	let intervalLength = 100/intervals.length;

	function getLogSliderValue(value) {
		let closestIndex = 0;
		let closestDiff = Infinity;

		for (const [index, interval] of intervals.entries()) {
			let diff = Math.abs(value-interval);
			if (diff < closestDiff) {
				closestDiff = diff;
				closestIndex = index;
			}
		}

		return (closestIndex+1) * intervalLength;
	}

	// Logarithmic
	let diffusionOutput = document.getElementById("diffusionOutput");
	let diffusionSlider = document.getElementById("diffusionSlider");
	diffusionOutput.value = diffusion;
	diffusionSlider.value = getLogSliderValue(diffusion);
	diffusionSlider.addEventListener("change", (event) => {
		grid.diffusion = intervals[Math.max(Math.ceil(event.target.value/intervalLength)-1,0)];
		diffusion = intervals[Math.max(Math.ceil(event.target.value/intervalLength)-1,0)];
	});
	diffusionSlider.addEventListener("input", (event) => {
		diffusionOutput.value = intervals[Math.max(Math.ceil(event.target.value/intervalLength)-1,0)];
	});

	// Logarithmic
	let viscosityOutput = document.getElementById("viscosityOutput");
	let viscositySlider = document.getElementById("viscositySlider");
	viscosityOutput.value = viscosity;
	viscositySlider.value = getLogSliderValue(viscosity);
	viscositySlider.addEventListener("change", (event) => {
		grid.viscosity = intervals[Math.max(Math.ceil(event.target.value/intervalLength)-1,0)];
		viscosity = intervals[Math.max(Math.ceil(event.target.value/intervalLength)-1,0)];
	});
	viscositySlider.addEventListener("input", (event) => {
		viscosityOutput.value = intervals[Math.max(Math.ceil(event.target.value/intervalLength)-1,0)];
	});

	// Logarithmic
	let fadeRateOutput = document.getElementById("fadeRateOutput");
	let fadeRateSlider = document.getElementById("fadeRateSlider");
	fadeRateOutput.value = fadeRate;
	fadeRateSlider.value = getLogSliderValue(fadeRate);
	fadeRateSlider.addEventListener("change", (event) => {
		grid.fadeRate = intervals[Math.max(Math.ceil(event.target.value/intervalLength)-1,0)];
		fadeRate = intervals[Math.max(Math.ceil(event.target.value/intervalLength)-1,0)];
	});
	fadeRateSlider.addEventListener("input", (event) => {
		fadeRateOutput.value = intervals[Math.max(Math.ceil(event.target.value/intervalLength)-1,0)];
	});

	let iterationsOutput = document.getElementById("iterationsOutput");
	let iterationsSlider = document.getElementById("iterationsSlider");
	iterationsOutput.value = iterations;
	iterationsSlider.value = iterations;
	iterationsSlider.addEventListener("change", (event) => {
		grid.iterations = event.target.value;
		iterations = event.target.value;
	});
	iterationsSlider.addEventListener("input", (event) => {
		iterationsOutput.value = event.target.value;
	});

	let densityOutput = document.getElementById("densityOutput");
	let densitySlider = document.getElementById("densitySlider");
	densityOutput.value = density;
	densitySlider.value = density;
	densitySlider.addEventListener("change", (event) => {
		density = event.target.value;
	});
	densitySlider.addEventListener("input", (event) => {
		densityOutput.value = event.target.value;
	});

	let velocityOutput = document.getElementById("velocityOutput");
	let velocitySlider = document.getElementById("velocitySlider");
	velocityOutput.value = velocity;
	velocitySlider.value = velocity;
	velocitySlider.addEventListener("change", (event) => {
		velocity = event.target.value;
	});
	velocitySlider.addEventListener("input", (event) => {
		velocityOutput.value = event.target.value;
	});

	let resolutionOutput = document.getElementById("resolutionOutput");
	let resolutionSlider = document.getElementById("resolutionSlider");
	resolutionOutput.value = resolution;
	resolutionSlider.value = resolution;
	resolutionSlider.addEventListener("change", (event) => {
		resolution = event.target.value;

		// Recreate frame and grid
		frame = context.createImageData(resolution, resolution);
			for (let y = 0; y < resolution; y++) {
				for (let x = 0; x < resolution; x++) {
					frame.data[x*4 + y*4*resolution + 3] = 255;
				}
			}
		grid = new FluidGrid(resolution, diffusion, viscosity, fadeRate, iterations);
	});
	resolutionSlider.addEventListener("input", (event) => {
		resolutionOutput.value = event.target.value;
	});

	let timeScaleOutput = document.getElementById("timeScaleOutput");
	let timeScaleSlider = document.getElementById("timeScaleSlider");
	timeScaleOutput.value = timeScale;
	timeScaleSlider.value = timeScale;
	timeScaleSlider.addEventListener("change", (event) => {
		timeScale = event.target.value;
	});
	timeScaleSlider.addEventListener("input", (event) => {
		timeScaleOutput.value = event.target.value;
	});



	function scale(x, in_min, in_max, out_min, out_max) {
		return(x - in_min) * (out_max - out_min) / (in_max - in_min) + out_min;
	}
}

