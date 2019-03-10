"use strict";

function main() {
	let canvas = document.getElementById("canvas");
	let context = canvas.getContext("2d");

	let N = 50;
	let SCALE = 8;

	// n, diffusion, viscosity, dt
	let grid = new FluidGrid(N, 0.0001, 0.01, 0.1);

	function draw() {
		context.clearRect(0, 0, canvas.width, canvas.height);
		for (let i = 0; i < N; i++) {
			for (let j = 0; j < N; j++) {
				let color = 255 * grid.density[i+j*N];
				context.fillStyle = "rgba("+color+","+color+","+color+",255)";
				context.fillRect(i*SCALE, j*SCALE, SCALE, SCALE);
			}
		}
		grid.stepForward();
		window.requestAnimationFrame(draw);
	}

	window.requestAnimationFrame(draw);

	canvas.addEventListener("mousemove", (event) => {
		grid.addDensity(Math.floor(event.clientX/SCALE)-1, Math.floor(event.clientY/SCALE)-1, 10);
		//grid.addVelocity(Math.floor(event.clientX/SCALE)-1, Math.floor(event.clientY/SCALE)-1, 10);
	});
	/*
	let imageData = context.createImageData(400, 400);
	let data = imageData.data;
	let steps = 0;
	function step() {
		context.clearRect(0, 0, canvas.width, canvas.height);
		context.putImageData(imageData, steps%400, 100);
		steps += 2;
		window.requestAnimationFrame(step);
	}

	window.requestAnimationFrame(step);
	*/
}

