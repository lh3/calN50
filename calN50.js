#!/usr/bin/env k8

var version = "r4";

var getopt = function(args, ostr) {
	var oli; // option letter list index
	if (typeof(getopt.place) == 'undefined')
		getopt.ind = 0, getopt.arg = null, getopt.place = -1;
	if (getopt.place == -1) { // update scanning pointer
		if (getopt.ind >= args.length || args[getopt.ind].charAt(getopt.place = 0) != '-') {
			getopt.place = -1;
			return null;
		}
		if (getopt.place + 1 < args[getopt.ind].length && args[getopt.ind].charAt(++getopt.place) == '-') { // found "--"
			++getopt.ind;
			getopt.place = -1;
			return null;
		}
	}
	var optopt = args[getopt.ind].charAt(getopt.place++); // character checked for validity
	if (optopt == ':' || (oli = ostr.indexOf(optopt)) < 0) {
		if (optopt == '-') return null; //  if the user didn't specify '-' as an option, assume it means null.
		if (getopt.place < 0) ++getopt.ind;
		return '?';
	}
	if (oli+1 >= ostr.length || ostr.charAt(++oli) != ':') { // don't need argument
		getopt.arg = null;
		if (getopt.place < 0 || getopt.place >= args[getopt.ind].length) ++getopt.ind, getopt.place = -1;
	} else { // need an argument
		if (getopt.place >= 0 && getopt.place < args[getopt.ind].length)
			getopt.arg = args[getopt.ind].substr(getopt.place);
		else if (args.length <= ++getopt.ind) { // no arg
			getopt.place = -1;
			if (ostr.length > 0 && ostr.charAt(0) == ':') return ':';
			return '?';
		} else getopt.arg = args[getopt.ind]; // white space
		getopt.place = -1;
		++getopt.ind;
	}
	return optopt;
}

function parseNum(s) {
	var m, x = null;
	if ((m = /^(\d*\.?\d*)([mMgGkK]?)/.exec(s)) != null) {
		x = parseFloat(m[1]);
		if (m[2] == 'k' || m[2] == 'K') x *= 1000;
		else if (m[2] == 'm' || m[2] == 'M') x *= 1000000;
		else if (m[2] == 'g' || m[2] == 'G') x *= 1000000000;
	}
	return Math.floor(x + .499);
}

function main(args) {
	var c, step = 0.1, min_len = 0, tot_len = null, fn_fai = null;
	while ((c = getopt(args, "l:L:s:f:v")) != null) {
		if (c == 's') step = parseFloat(getopt.arg);
		else if (c == 'l') min_len = parseNum(getopt.arg);
		else if (c == 'L') tot_len = parseNum(getopt.arg);
		else if (c == 'f') fn_fai = getopt.arg;
		else if (c == 'v') {
			print(version);
			exit(0);
		}
	}

	if (args.length == 0) {
		print("Usage: calN50.js [options] <in.fa>|<in.gfa>|<in.fai>");
		print("Options:");
		print("  -l NUM     min length [0]");
		print("  -L NUM     total length for NGx []");
		print("  -f FILE    reference .fai file for NGx []");
		print("  -s FLOAT   N50 step size [" + step + "]");
		print("  -v         print version number");
		exit(0);
	}

	var file, buf = new Bytes();
	if (fn_fai) {
		file = new File(fn_fai);
		tot_len = 0;
		while (file.readline(buf) >= 0) {
			var t = buf.toString().split("\t");
			if (t.length >= 2)
				tot_len += parseInt(t[1]);
		}
		file.close();
	}

	file = args[getopt.ind] == '-'? new File() : new File(args[getopt.ind]);

	var is_fa = false, is_gfa = false, len = 0, name = null;
	var a = [];
	while (file.readline(buf) >= 0) {
		if (buf.length == 0) continue;
		var m, s = buf.toString();
		if (s[0] == '>') { // fasta header
			if ((m = /^>(\S+)/.exec(s)) != null) {
				if (name) a.push([name, len]);
				is_fa = true, name = m[1], len = 0;
			}
		} else if (is_fa) { // fasta sequence line
			len += s.length;
		} else { // gfa or length line
			if ((m = /^S\t(\S+)\t([a-zA-Z]+)|(\*.*\tLN:i:(\d+))/.exec(s)) != null) { // GFA S-line
				if (m[4] != null || m[2] != null) {
					is_gfa = true;
					if (m[4] != null) a.push([m[1], parseInt(m[4])]);
					else a.push([m[1], m[2].length]);
				}
			} else if (!is_gfa) {
				if ((m = /^(\S+)\t(\d+)/.exec(s)) != null) // length line
					a.push([m[1], parseInt(m[2])]);
			}
		}
	}
	if (is_fa && name && len)
		a.push([name, len]);

	file.close();
	buf.destroy();

	if (a.length == 0) {
		warn("ERROR: no sequences found");
		return 1;
	}

	a.sort(function(x,y) { return y[1]-x[1] });
	if (min_len > 0) {
		var j = a.length;
		for (var i = a.length - 1; i >= 0; --i)
			if (a[i][1] >= min_len) {
				j = i;
				break;
			}
		a.length = j + 1;
	}

	print("CC\tGS   genome_size_if_provided");
	print("CC\tSZ   total_sequence_length");
	print("CC\tNN   number_of_sequences");
	print("CC\tNL   x   Nx   Lx");
	print("CC\tAU   auN");
	print("CC");

	var sum = 0;
	for (var i = 0; i < a.length; ++i)
		sum += a[i][1];
	if (tot_len != null) print("GS", tot_len);
	print("SZ", sum);
	print("NN", a.length);
	if (tot_len != null) sum = tot_len;

	var n = 0, x = 0, next = 0, y = 0;
	for (var i = 0; i < a.length; ++i) {
		if (x >= sum) break;
		var l = x + a[i][1] <= sum? a[i][1] : sum - x;
		y += l * (l / sum);
		x += a[i][1], ++n;
		if (x > sum * next - 0.01) {
			do {
				print("NL", Math.floor(next*100.0+.499), a[i][1], n);
				next += step;
			} while (x > sum * next - 0.01);
		}
	}
	print("AU", y.toFixed(0));
	return 0;
}

var ret = main(arguments);
exit(ret);
