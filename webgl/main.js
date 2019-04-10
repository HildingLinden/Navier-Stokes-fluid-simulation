"use strict";

let vertexShaderCode = `#version 300 es
in vec2 a_position;
in vec2 a_textureCoord;

out vec2 v_textureCoord;

uniform vec2 u_resolution;

void main() {
	vec2 clipSpace = ((a_position/u_resolution) * 2.0) - 1.0;
	gl_Position = vec4(clipSpace, 0, 1);

	v_textureCoord = a_textureCoord;
}`;

let fragmentShaderCode = `#version 300 es
precision mediump float;

in vec2 v_textureCoord;

out vec4 color;

uniform sampler2D u_texture;

void main() {
	// color = texelFetch(u_texture, ivec2(v_textureCoord), 0);
	color =
		(texelFetch(u_texture, ivec2(v_textureCoord) + ivec2( 0,-1), 0) +
		 texelFetch(u_texture, ivec2(v_textureCoord) + ivec2( 1, 0), 0) +
		 texelFetch(u_texture, ivec2(v_textureCoord) + ivec2( 0, 1), 0) +
		 texelFetch(u_texture, ivec2(v_textureCoord) + ivec2(-1, 0), 0));
}
`;

let N = 10;

function loadImage() {
	let img = new Image();
	img.addEventListener("load", () => {
		console.log("Image loaded");
		main(img);
	}, false);
	img.src = "test.png";
}

function main(img) {
	// Canvas and context setup
	let canvas = document.getElementById("canvas");
	let context = canvas.getContext("webgl2",
		{antialias: false, powerPreference: "high-performance"});

	// Enable 32-bit floating-point values for the color buffer
	context.getExtension("EXT_color_buffer_float");
	context.pixelStorei(context.UNPACK_ALIGNMENT, 1);

	context.clearColor(0, 0, 0, 0);

	// Shader and program setup
	let vertexShader = createShader(context, context.VERTEX_SHADER, vertexShaderCode);
	let fragmentShader = createShader(context, context.FRAGMENT_SHADER, fragmentShaderCode);

	let program = linkShaders(context, vertexShader, fragmentShader);
	context.useProgram(program);

	// Attribute/Uniform setup
	let vao = context.createVertexArray();
	context.bindVertexArray(vao);

	let positionAttributeLocation = context.getAttribLocation(program, "a_position");
	let positionBuffer = context.createBuffer();
	context.bindBuffer(context.ARRAY_BUFFER, positionBuffer);
	let positions = new Float32Array([
		 1, 1,
		 1,101,
		101,101,
		101,101,
		101, 1,
		 1, 1
	]);
	context.bufferData(context.ARRAY_BUFFER, positions, context.STATIC_DRAW);
	context.enableVertexAttribArray(positionAttributeLocation);
	context.vertexAttribPointer(positionAttributeLocation, 2, context.FLOAT, false, 0, 0);

	let textureCoordAttributeLocation = context.getAttribLocation(program, "a_textureCoord");
	let textureCoordBuffer = context.createBuffer();
	context.bindBuffer(context.ARRAY_BUFFER, textureCoordBuffer);
	let textureCoord = new Float32Array([
		0, 0,
		0, N,
		N, N,
		N, N,
		N, 0,
		0, 0
	]);
	context.bufferData(context.ARRAY_BUFFER, textureCoord, context.STATIC_DRAW);
	context.enableVertexAttribArray(textureCoordAttributeLocation);
	context.vertexAttribPointer(textureCoordAttributeLocation, 2, context.FLOAT, false, 0, 0);

	let resolutionUniformLocation = context.getUniformLocation(program, "u_resolution");
	context.uniform2f(resolutionUniformLocation, context.canvas.width, context.canvas.height);

	/*
	 * Textures/Data
	 */
	// Initial data
	let initialData = textureSetup(context);
	let textureArr = new Uint8Array(N*N*4);
	textureArr[(N*N/2+N/2)*4] = 255;
	context.texImage2D(context.TEXTURE_2D, 0, context.RGBA, N, N, 0, context.RGBA, context.UNSIGNED_BYTE, textureArr);
	// R32F, RED, FLOAT

	let textures = [];
	let framebuffers = [];
	// For each intermediate step
	for (let i = 0; i < 2; i++) {
		let texture = textureSetup(context);
		textures.push(texture);

		// Empty texture
		context.texImage2D(context.TEXTURE_2D, 0, context.RGBA, N, N, 0, context.RGBA, context.UNSIGNED_BYTE, null);

		let fb = context.createFramebuffer();
		framebuffers.push(fb);
		context.bindFramebuffer(context.FRAMEBUFFER, fb);

		let attachment = context.COLOR_ATTACHMENT0;
		context.framebufferTexture2D(context.FRAMEBUFFER, attachment, context.TEXTURE_2D, texture, 0);
	}

	let textureUniformLocation = context.getUniformLocation(program, "u_texture");

	draw();

	function draw() {
		performance.mark("start");

		// Use and bind all the correct things
		context.useProgram(program);
		context.bindVertexArray(vao);
		context.activeTexture(context.TEXTURE0);
		context.bindTexture(context.TEXTURE_2D, initialData);
		context.uniform1i(textureUniformLocation, 0);

		let count = 0;
		// Draw to the first framebuffer with the initial data
		//context.clear(context.COLOR_BUFFER_BIT | context.DEPTH_BUFFER_BIT);
		context.bindFramebuffer(context.FRAMEBUFFER, framebuffers[0]);
		context.viewport(0, 0, N, N);
		context.drawArrays(context.TRIANGLES, 0, 6);

		// Then alternate textures and framebuffers
		// For all steps
		context.bindTexture(context.TEXTURE_2D, textures[count % 2]);
		context.bindFramebuffer(context.FRAMEBUFFER, framebuffers[++count % 2]);
		context.drawArrays(context.TRIANGLES, 0, 6);
		// End for

		// Finally draw to canvas
		context.bindTexture(context.TEXTURE_2D, textures[count % 2]);
		context.bindFramebuffer(context.FRAMEBUFFER, null);
		context.viewport(0, 0, context.canvas.width, context.canvas.height);
		context.drawArrays(context.TRIANGLES, 0, 6);

		// Only accurate to 2ms (FF60) or 20µs (FF59 or if privacy.reduceTimerPrecision is false)
		performance.measure("draw", "start");
		console.log("Time for draw: " + performance.getEntriesByType("measure")[0]["duration"] + "ms");

		// Read the drawn pixels for debugging
		// let results = new Uint8Array(context.drawingBufferWidth * context.drawingBufferHeight * 4);
		// context.readPixels(0, 0, context.drawingBufferWidth, context.drawingBufferHeight, context.RGBA, context.UNSIGNED_BYTE, results);
		// results.forEach((element, index) => {
		// 	if (element > 0) {
		// 		console.log(index);
		// 	}
		// });
		// console.log(results);
	}
}



function createShader(context, type, source) {
	let shader = context.createShader(type);
	context.shaderSource(shader, source);
	context.compileShader(shader);
	let success = context.getShaderParameter(shader, context.COMPILE_STATUS);
	if (success) {
		return shader;
	}

	let errorMsg = context.getShaderInfoLog(shader);
	context.deleteShader(shader);
	throw new Error("in shader, " + errorMsg + " \n" + source);
}

function linkShaders(context, vertexShader, fragmentShader) {
	let program = context.createProgram();
	context.attachShader(program, vertexShader);
	context.attachShader(program, fragmentShader);
	context.linkProgram(program);
	let success = context.getProgramParameter(program, context.LINK_STATUS);
	if (success) {
		return program;
	}

	let errorMsg = context.getProgramInfoLog(program);
	context.deleteProgram(program);
	throw new Error("in linkShaders: " + errorMsg);
}

function textureSetup(context) {
	let texture = context.createTexture();

	context.bindTexture(context.TEXTURE_2D, texture);

	context.texParameteri(context.TEXTURE_2D, context.TEXTURE_WRAP_S, context.CLAMP_TO_EDGE);
	context.texParameteri(context.TEXTURE_2D, context.TEXTURE_WRAP_T, context.CLAMP_TO_EDGE);
	context.texParameteri(context.TEXTURE_2D, context.TEXTURE_MIN_FILTER, context.NEAREST);
	context.texParameteri(context.TEXTURE_2D, context.TEXTURE_MAG_FILTER, context.NEAREST);

	return texture;
}