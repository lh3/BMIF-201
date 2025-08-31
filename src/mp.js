#!/usr/bin/env k8

"use strict";

class PKNode {
	constructor() {
		this.id = -1;
		this.parent = null, this.child = [];
		this.name = this.meta = this.seq = "";
		this.d = -1.0; // branch length
	}
}

class PKTree {
	constructor() {
		this.root = null, this.node = [];
		this.error = 0; // 1: unmatched ); 2: unmatched (; 4: unmatched ]; 8: duplicated leaves
	}
}

function pk_parse_nh(str) {
	const re = /(\(|((\)?[^,;:\[\]\(\)]+|\))(:[\d.eE\-]+)?(\[[^\[\]]*\])?))/g;
	const str_compact = str.replace(/[;\s\n]+/g, ""); // compacted string without SPACE, new line or ';'
	let tree = new PKTree(), m, stack = [], leaf_cnt = {};
	while ((m = re.exec(str_compact)) != null) { // each match returns '(' or a node
		let p = new PKNode();
		if (m[1] == '(') {
			p.id = -2; // a placehold for the left '('
		} else if (m[3]) {
			if (m[3][0] == ')') { // an internal node
				let q;
				while (stack.length > 0) { // look for the matching '('
					q = stack.pop();
					if (q.id == -2) break;
					p.child.push(q);
				}
				p.child.reverse();
				if (q.id != -2) {
					tree.error |= 1; // ERROR: unmatched ')'
					break;
				}
				for (const q of p.child) q.parent = p; // set parent node
				if (m[3].length > 1) p.name = m[3].substr(1); // set internal name if present
			} else { // a leaf
				p.name = m[3];
				if (leaf_cnt[p.name] == null) leaf_cnt[p.name] = 0;
				if (++leaf_cnt[p.name] > 1) tree.error |= 8; // ERROR: duplicated leaf name
			}
			if (m[4]) p.d = parseFloat(m[4].substr(1)); // set branch length if present
			p.id = tree.node.length;
			tree.node.push(p);
		}
		stack.push(p);
	}
	if (stack.length != 1) tree.error |= 2; // ERROR: unmatched '('
	tree.root = tree.node[tree.node.length - 1];
	return tree;
}

function pk_parse_msa_simple(tree, msa_str) {
	let h = {};
	for (let i = 0; i < tree.node.length; ++i) {
		const p = tree.node[i];
		if (p.child.length == 0)
			h[p.name] = i;
	}
	for (const line of msa_str.split("\n")) {
		let t = line.split(/\s+/);
		if (t.length < 2) return;
		if (h[t[0]] == null) return;
		tree.node[h[t[0]]].seq += t[1];
	}
}

PKTree.prototype.seq_len = function() {
	let len = -1;
	for (let i = 0; i < this.node.length; ++i) {
		const p = this.node[i];
		if (p.child.length == 0) {
			if (len < 0) len = p.seq.length;
			else if (len != p.seq.length) return -1;
		}
	}
	return len;
}

function pk_mp1(node, col) {
	let R = [], c = 0;
	for (let i = 0; i < node.length; ++i) R[i] = {};
	for (let i = 0; i < node.length; ++i) {
		const p = node[i];
		if (p.child.length == 0) { // leaf
			R[i][p.seq[col]] = 1;
		} else { // internal node
			let h = {};
			for (let j = 0; j < p.child.length; ++j) {
				const Rj = R[p.child[j].id];
				for (const x in Rj) {
					if (h[x] == null) h[x] = 0;
					++h[x];
				}
			}
			let inter = {}, n_inter = 0;
			for (const x in h)
				if (h[x] == p.child.length)
					inter[x] = 1, ++n_inter;
			if (n_inter > 0) R[i] = inter;
			else R[i] = h, ++c;
		}
	}
	let a = [];
	for (const x in R[node.length-1])
		a.push(x);
	return [c, a];
}

function main(args) {
	let msa =
`human TCCTGCCTCATCCTATTATTTATCGCACCTAC-GTTCAATATTACAGGCGAA-CATA-CTTACTAAAGTGTGTTAATTAATTAATGCTTGTAG
bonobo TCCTGCCCCATTACGTTATTTATCGCACCTAC-GTTCAATATTATTACCTAG-CATGATTTACTAAAGCGTGTTAATTAATTAATGCTTGTAG
chimp TCCTGCCCCATTGTATTATTTATCGCACCTAC-GTTCAATATTACGACCTAG-CATA-CCTACTAAAGTGTGTTGATTAATTAATGCTTGCAG
gorilla TCCTGCCCCATGCTACCATTTATCGCACCTAC-GTTCAATATTACAGCCGAG-CGCA-CAGTGTTCATGGTGTTAATTAATTCATGCTTGTTG
oran-pa TCCTACCTCATGCCATTATTAATCGCGCCTAATATCCAATATCCTAGCCCCACCCTC-AGTGTTTGAAGCTGCTATTTAATTTATGCTAG-AG
oran-pp TCCTGCCCCATGGCGTTATTGATCGCGCCTAACGTCCAATGTTCTAGCGCCC-CCTC-CCTATTGAAAGTTGTTATTTAATTTATGCTAG-AG
gibbon TTCTGACCCATCCTATTGTTGATCGCGCCTAC-GTTCAATATCCCAGCCGAG-CATA-CTTACACTAAGGTGTTAATTAATTCATGCTTGTTG`;
	let trees = [
		"((((human,(chimp,bonobo)),gorilla),(oran-pa,oran-pp)),gibbon)",
		"((((human,gorilla),(chimp,bonobo)),(oran-pa,oran-pp)),gibbon)",
		"((((human,bonobo),gorilla),(oran-pa,oran-pp)),(chimp,gibbon))",
		"(human,(((gibbon,(oran-pa,oran-pp)),gorilla),(chimp,bonobo)))"];
	for (let i = 0; i < trees.length; ++i) {
		let tree = pk_parse_nh(trees[i]);
		pk_parse_msa_simple(tree, msa);
		let tot = 0;
		for (let i = 0; i < tree.seq_len(); ++i) {
			let [c, a] = pk_mp1(tree.node, i);
			tot += c;
		}
		if (typeof print == "function") print(tot, trees[i]);
		else console.log(tot, trees[i]);
	}
}

main(arguments);
