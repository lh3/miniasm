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

function pafmask(a, mask_level)
{
	var k = 1;
	for (var i = 1; i < a.length; ++i) {
		var j, ai = a[i];
		for (j = 0; j < k; ++j) {
			var ol = 0, aj = a[j];
			if (ai[2] < aj[2]) {
				if (ai[3] > aj[2])
					ol = ai[3] - aj[2];
			} else {
				if (aj[3] > ai[2])
					ol = aj[3] - ai[2];
			}
			var min_l = ai[3] - ai[2] < aj[3] - aj[2]? ai[3] - ai[2] : aj[3] - aj[2];
			if (ol > min_l * mask_level)
				break;
		}
		if (j == k) a[k++] = ai;
	}
	a.length = k;
}

function pafmerge(a, max_gap)
{
	for (var i = 1; i < a.length; ++i) {
		var ai = a[i];
		for (var j = 0; j < i; ++j) {
			var aj = a[j];
			if (aj[4] != ai[4] || aj[5] != ai[5]) continue; // diff strand or chr
			var ts = [ai[7], aj[7]], te = [ai[8], aj[8]];
			var qs = [ai[2], aj[2]], qe = [ai[3], aj[3]];
			if (qs[0] > qs[1]) {
				qs = [aj[2], ai[2]], qe = [aj[3], ai[3]];
				ts = [aj[7], ai[7]], te = [aj[8], ai[8]];
				if (ai[4] == '-') {
					ts = [aj[6] - aj[8], ai[6] - ai[8]];
					te = [aj[6] - aj[7], ai[6] - ai[7]];
				}
			} else {
				if (ai[4] == '-') {
					ts = [ai[6] - ai[8], aj[6] - aj[8]];
					te = [ai[6] - ai[7], aj[6] - aj[7]];
				}
			}
			if (qe[0] > qe[1]) continue; // contained
			if (ts[0] > ts[1]) continue;
			var qg = qs[1] - qe[0], tg = ts[1] - te[0];
			if ((qg < 0 && tg < 0) || Math.abs(tg - qg) < max_gap) {
				//print("Merged: ["+ai[2]+","+ai[3]+") <=> ["+aj[2]+","+aj[3]+") "+ai[4]+" ["+ai[7]+","+ai[8]+") <=> ["+aj[7]+","+aj[8]+")");
				aj[2] = qs[0], aj[3] = qe[1];
				if (aj[4] == '+') {
					aj[7] = ts[0], aj[8] = te[1];
				} else {
					aj[7] = aj[6] - te[1], aj[8] = aj[6] - ts[0];
				}
				aj[9] += ai[9], aj[10] += ai[10];
				aj[11] = aj[11] > ai[11]? aj[11] : ai[11];
				a[i] = [];
				break;
			}
		}
	}
	var k = 0;
	for (var i = 0; i < a.length; ++i)
		if (a[i].length != 0) a[k++] = a[i];
	a.length = k;
}

function paftop(a, mask_level, max_gap)
{
	for (var i = 0; i < a.length; ++i) {
		for (var j = 1; j <= 3; ++j) a[i][j] = parseInt(a[i][j]);
		for (var j = 6; j <= 11; ++j) a[i][j] = parseInt(a[i][j]);
	}
	a.sort(function(x,y){return y[9]-x[9];});
	pafmask(a, mask_level);
	pafmerge(a, max_gap);
	pafmask(a, mask_level);
	for (var i = 0; i < a.length; ++i)
		if (a[i].length) print(a[i].join("\t"));
}

var c, mask_level = .5, max_gap = 1000;
while ((c = getopt(arguments, 'm:g:')) != null)
	if (c == 'm') mask_level = parseFloat(getopt.arg);
	else if (c == 'g') max_gap = parseInt(getopt.arg);

var file = arguments.length == getopt.ind? new File() : new File(arguments[getopt.ind]);
var buf = new Bytes();

var last = null, a = [];
while (file.readline(buf) >= 0) {
	var t = buf.toString().split("\t");
	if (t[0] != last) {
		if (a.length) paftop(a, mask_level, max_gap);
		a = [], last = t[0];
	}
	a.push(t);
}
if (a.length) paftop(a, mask_level, max_gap);

buf.destroy();
file.close();
