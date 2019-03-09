"use strict";

function main() {
	let canvas = document.getElementById("canvas");
	let context = canvas.getContext("2d");

	let N = 40;
	let SCALE = 10;

	// n, diffusion, viscosity, dt
	let grid = new FluidGrid(N, 0.01, 0.01, 0.001);
	grid.addDensity(20,20, 1000);
	//grid.addVelocity(20,20, 1000, 0);
	//console.log(grid.density[21+20*N]);

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

