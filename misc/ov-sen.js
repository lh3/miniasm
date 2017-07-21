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

var c, min_len = 2000, min_mapq = 10;
while ((c = getopt(arguments, "l:q:")) != null)
	if (c == 'l') min_len = parseInt(getopt.arg);
	else if (c == 'q') min_mapq = parseInt(getopt.arg);

if (arguments.length - getopt.ind < 2) {
	warn("Usage: k8 ov-sen.js [options] <in.ref-sorted.paf> <in.ovlp.paf>");
	warn("Options:");
	warn("  -l INT     min overlap length [2000]");
	warn("  -q INT     min mapping quality [10]");
	exit(1);
}

var file, buf = new Bytes();

var h = {};

file = new File(arguments[getopt.ind]);
var a = [];
while (file.readline(buf) >= 0) {
	var t = buf.toString().split("\t");
	if (parseInt(t[11]) < min_mapq) continue;
	if (parseInt(t[10]) < min_len) continue;
	var st = parseInt(t[7]), en = parseInt(t[8]);

	var n_shift = 0;
	if (a.length > 0) {
		for (var i = 0; i < a.length; ++i) {
			if (t[5] != a[i][1]) {
				++n_shift;
			} else {
				var min_en = a[i][3] < en? a[i][3] : en;
				if (min_en - st >= min_len) break;
				++n_shift;
			}
		}
	}
	if (n_shift > 0) {
		for (var i = 0; i < n_shift; ++i)
			a.shift();
	}
	if (a.length > 0) {
		for (var i = 0; i < a.length; ++i) {
			if (t[5] != a[i][1]) continue;
			var min_en = a[i][3] < en? a[i][3] : en;
			if (min_en - st < min_len) continue;
			h[a[i][0] + "\t" + t[0]] = 0;
			//print(a[i][0], t[0], min_en - st);
		}
	}
	a.push([t[0], t[5], st, en]);
}
file.close();

file = new File(arguments[getopt.ind+1]);
while (file.readline(buf) >= 0) {
	var t = buf.toString().split("\t");
	var key = t[0] + "\t" + t[5];
	if (h[key] != null) ++h[key];
	else {
		key = t[5] + "\t" + t[0];
		if (h[key] != null) ++h[key];
	}
}
file.close();

buf.destroy();

var n_ovlp = 0, n_missed = 0;
for (var key in h) {
	++n_ovlp;
	if (h[key] == 0) ++n_missed;
}
print(n_ovlp + " overlaps");
print(n_missed + " missed");
print((1 - n_missed/n_ovlp).toFixed(4) + " sensitivity");
