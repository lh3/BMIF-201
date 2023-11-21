#!/usr/bin/env k8

const data_str = // data
`RRAAAR.
AA.RRAR
ARARARA
.ARRAAR
.ARRRAA`;

function parse_data(data_str) { // parse input data string
	let data = [];
	for (const s of data_str.split(/[\n\s]+/)) { // traverse each line
		let t = [];
		for (let i = 0; i < s.length; ++i)
			t[i] = s[i] === "R"? 1 : s[i] === "A"? -1 : 0; // R=1, A=-1, .=0
		data.push(t);
		if (t.length != data[0].length)
			throw Error("Unequal lengths");
	}
	return data;
}

function random_delta(m) { // generate random delta (SNP phase)
	let delta = [];
	for (let i = 0; i < m; ++i)
		delta[i] = Math.random() < .5? -1 : 1;
	return delta;
}

function cal_sigma(data, delta) { // compute sigma (read phase) from delta
	let sigma = [], m0 = 0, m1 = 0;
	for (let k = 0; k < data.length; ++k) {
		const t = data[k];
		for (let i = 0; i < t.length; ++i) {
			if (t[i] == delta[i]) ++m0;
			else if (t[i] == -delta[i]) ++m1;
		}
		sigma[k] = m0 > m1? 1 : m0 < m1? -1 : Math.random() < 0.5? 1 : -1;
	}
	return sigma;
}

function flip_delta(data, sigma, delta) { // optimize delta
	let n_flip = 0, n_err = 0;
	for (let i = 0; i < delta.length; ++i) {
		let m0 = 0, m1 = 0;
		let x = [];
		for (let k = 0; k < sigma.length; ++k) {
			x.push(data[k][i]);
			if (data[k][i] == 0) continue;
			if (data[k][i] == delta[i] * sigma[k]) ++m0; // m0: number matches
			else ++m1; // m1: number of mismatches
		}
		if (m0 < m1 || (m0 == m1 && Math.random() < .5))
			delta[i] = -delta[i], n_flip += m1 - m0; // flip delta[i]
		n_err += m0 < m1? m0 : m1;
	}
	return [n_flip, n_err]; // return the number of flips and the number of remaining errors
}

function phase(data) {
	let delta = random_delta(data[0].length);
	let sigma = cal_sigma(data, delta);
	while (1) {
		const [n_flip, n_err] = flip_delta(data, sigma, delta);
		if (n_flip == 0) // no flip; local optimum
			return [n_err, delta]
		sigma = cal_sigma(data, delta);
	} 
}

const n_try = 100;
let data = parse_data(data_str);
for (let i = 0; i < n_try; ++i) {
	let [n_err, delta] = phase(data);
	print(n_err, delta.join(","));
}
