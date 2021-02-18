"use strict";

function main() {
	let canvas = document.getElementById("canvas");
	let context = canvas.getContext("2d");
	context.imageSmoothingEnabled = false;

	let stepTimeText = document.getElementById("stepTime");
	let drawTimeText = document.getElementById("drawTime");

	let N = 200;
	let SCALE = canvas.width/N;

	let diffusion = 0.001;
	let viscosity = 0.00001;
	let fadeRate = 0.00;
	let iterations = 4;

	let grid = new FluidGrid(N, diffusion, viscosity, fadeRate, iterations);

	let velX = 0;
	let velY = 0;

	// Create inital frame and set alpha to 255
	let frame = context.createImageData(N, N);
	for (let y = 0; y < N; y++) {
		for (let x = 0; x < N; x++) {
			frame.data[x*4 + y*4*N + 3] = 255;
		}
	}

	let dt = 0;

	function loop() {
		// Add velocity to the middle of the canvas in the direction of the mouse
		grid.addVelocity(N/2, N/2, velX, velY, dt);

		// Add density to the middle of the canvas
		grid.addDensity(N/2, N/2, 100, dt);

		// Step the physics forward a time step
		let stepTime = grid.step(dt);

		stepTimeText.innerHTML = "Physics step time " + stepTime.toFixed(3) + "ms";

		let t0 = performance.now();

		// Reset canvas
		context.clearRect(0, 0, canvas.width, canvas.height);

		// Convert density into grayscale pixels in a frame
		for (let y = 0; y < N; y++) {
			for (let x = 0; x < N; x++) {
				let color = 255 * grid.density[x+y*N];
				frame.data[x*4 + y*4*N] = color;
				frame.data[x*4 + y*4*N + 1] = color;
				frame.data[x*4 + y*4*N + 2] = color;
			}
		}

		// Render frame on canvas then scale it up to the correct size
		context.putImageData(frame, 0, 0);
		context.drawImage(canvas, 0, 0, N, N, 0, 0, canvas.width, canvas.height);

		let t1 = performance.now();
		let drawTime = t1-t0;
		drawTimeText.innerHTML = "Draw time " + drawTime.toFixed(3) + "ms";

		dt = (stepTime + drawTime) / 1000;

		window.requestAnimationFrame(loop);
	}

	window.requestAnimationFrame(loop);

	canvas.addEventListener("mousemove", (event) => {
		let x = event.clientX - canvas.width/2 - 10;
		let y = event.clientY - canvas.height/2 - 10;
		let angle = Math.atan2(y, x);
		velX = 80 * Math.cos(angle);
		velY = 80 * Math.sin(angle);

	});

	function scale(x, in_min, in_max, out_min, out_max) {
		return(x - in_min) * (out_max - out_min) / (in_max - in_min) + out_min;
	}
}

