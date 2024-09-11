#!/usr/bin/env node

function aln_xdrop(a, b, pen, xdrop) {
	let max_sc = xdrop, max_i = -1, max_j = -1;
	let st = 0, H = [];
	for (let j = 0; j < b.length; ++j) { // initialize row -1
		const sc = xdrop - pen * (j + 1);
		if (sc <= 0) break;
		H[j] = sc;
	}
	for (let i = 0; i < a.length; ++i) {
		let H11 = st == 0? xdrop - pen * i : H[st - 1]; // H[i-1,j-1]
		let H01 = H11 - pen; // H[i,j-1] = H[i-1,j-1] - pen
		if (H01 <= 0) ++st; // ignored dropped cells
		let max_row_sc = 0; // max score on row i
		for (let j = st; j < b.length; ++j) {
			let H00 = H11 + (a[i] == b[j]? 1 : -pen); // H[i,j] = H[i-1,j-1] + s(a[i], b[j])
			if (j < H.length)
				H00 = H00 > H[j] - pen? H00 : H[j] - pen; // H[i,j] = max(H[i,j], H[i-1,j] - pen)
			H00 = H00 > H01 - pen? H00 : H01 - pen; // H[i,j] = max(H[i,j], H[i,j-1] - pen)
			if (j >= H.length && H00 <= max_sc - xdrop) break;
			max_row_sc = max_row_sc > H00? max_row_sc : H00; // keep the max row score
			if (H00 > max_sc) // keep the overall max score and the cell position
				max_sc = H00, max_i = i, max_j = j;
			H01 = H00;
			H11 = j < H.length? H[j] : 0;
			H[j] = H00;
		}
		if (max_row_sc <= max_sc - xdrop) break;
	}
	return [max_sc - xdrop, max_i, max_j];
}

let a =  "GTTGATGGTCTACAACGTTATCGTCACAGCCCATGCATTTGTAATAATCTTCTTCATAGTAATA";
let b = "GaTAGATGGTCTGAGCTATGATATCAATTGGCTTCCTAGGGTTTATCGTGTGAGCACACCATATT";
let [max_sc, max_i, max_j] = aln_xdrop(a, b, 2, 10);
console.log(max_sc, max_i, max_j, a.substr(0, max_i+1), b.substr(0, max_j+1));
