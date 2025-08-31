#!/usr/bin/env k8

"use strict";

let PhyloKit = {};

class PKNode {
	constructor() {
		this.id = -1;
		this.parent = null;
		this.child = [];
		this.name = this.meta = this.seq = "";
		this.d = -1.0; // branch length
		this.hidden = false; // hide the subtree descended from the node
	}
}

class PKTree {
	constructor() {
		this.root = null;
		this.node = [];
		this.error = 0; // 1: unmatched ); 2: unmatched (; 4: unmatched ]; 8: duplicated leaves
	}
}

// parse a tree in the NH format
PhyloKit.parse_nh = function(str) {
	// helper function: parse a single node
	function parse_node(str, l, tree, x) {
		let beg, end = 0, i;
		let z = new PKNode();
		for (i = l, beg = l; i < str.length && str[i] != ',' && str[i] != ')'; ++i) {
			const c = str[i];
			if (c == '[') { // NH comment; skip
				const meta_beg = i;
				if (end == 0) end = i;
				do { ++i; } while (i < str.length && str.charAt(i) != ']');
				if (i == str.length) {
					tree.error |= 4;
					break;
				}
				z.meta = str.substr(meta_beg, i - meta_beg + 1);
			} else if (c == ':') { // end of an internal node
				if (end == 0) end = i;
				let j;
				for (j = ++i; i < str.length; ++i) {
					const cc = str[i];
					if ((cc < '0' || cc > '9') && cc != 'e' && cc != 'E' && cc != '+' && cc != '-' && cc != '.')
						break;
				}
				z.d = parseFloat(str.substr(j, i - j));
				--i;
			} else { // external node
				if (c < '!' && c > '~' && end == 0) end = i;
			}
		}
		if (end == 0) end = i;
		if (end > beg) {
			z.name = str.substr(beg, end - beg);
			if (z.name == ";") z.name = "";
		}
		tree.node.push(z);
		return i;
	}

	let tree = new PKTree();
	let stack = [];
	for (let l = 0; l < str.length;) {
		while (l < str.length && (str.charCodeAt(l) < 33 || str.charCodeAt(l) > 126)) ++l;
		if (l == str.length) break;
		let c = str[l];
		if (c == ',') ++l;
		else if (c == '(') {
			stack.push(-1); ++l;
		} else if (c == ')') {
			let i;
			for (i = stack.length - 1; i >= 0; --i)
				if (stack[i] < 0) break;
			if (i < 0) {
				tree.error |= 1; break;
			}
			const x = tree.node.length;
			let m = stack.length - 1 - i;
			l = parse_node(str, l + 1, tree, m);
			for (i = stack.length - 1, m = m - 1; m >= 0; --m, --i) {
				tree.node[x].child[m] = tree.node[stack[i]];
				tree.node[stack[i]].parent = tree.node[x];
			}
			stack.length = i;
			stack.push(x);
		} else {
			stack.push(tree.node.length);
			l = parse_node(str, l, tree, 0);
		}
	}
	if (stack.length > 1) tree.error |= 2;
	tree.root = tree.node[tree.node.length - 1];
	for (let i = 0; i < tree.node.length; ++i) tree.node[i].id = i;
	return tree;
}

PKTree.prototype.get_leaf_dict = function() {
	let h = {};
	for (let i = 0; i < this.node.length; ++i) {
		const p = this.node[i];
		if (p.child.length == 0)
			h[p.name] = i;
	}
	return h;
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

PKTree.prototype.get_msa_col = function(col) {
	let a = [];
	for (let i = 0; i < this.node.length; ++i)
		a.push(this.node[i].seq[col]);
	return a;
}

class PKMsaReader {
	constructor(tree) {
		this.tree = tree;
		this.name2leaf = tree.get_leaf_dict();
	}
	parse_line(line) {
		let t = line.split(/\s+/);
		if (t.length < 2) return;
		if (this.name2leaf[t[0]] == null) return;
		this.tree.node[this.name2leaf[t[0]]].seq += t[1];
	}
}

PhyloKit.mp1 = function(node, col) {
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
		"(human,(((gibbon,(oran-pa,oran-pp)),gorilla),(chimp,bonobo)))"];
	for (let i = 0; i < trees.length; ++i) {
		let tree = PhyloKit.parse_nh(trees[i]);
		let r = new PKMsaReader(tree);
		for (const line of msa.split("\n"))
			r.parse_line(line);
		let tot = 0;
		for (let i = 0; i < tree.seq_len(); ++i) {
			let [c, a] = PhyloKit.mp1(tree.node, i);
			tot += c;
		}
		print(tot, trees[i]);
	}
}

main(arguments);
