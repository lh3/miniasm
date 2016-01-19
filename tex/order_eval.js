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

var c, ws = 5, min_span = 2000;

while ((c = getopt(arguments, "w:s:")) != null)
	if (c == 'w') ws = parseInt(getopt.arg);
	else if (c == 's') min_span = parseInt(getopt.arg);

if (arguments.length - getopt.ind < 2) {
	print("Usage: k8 cmp_order.js <gfa2bed.bed> <paftop>");
	exit(1);
}

var b = new Bytes();

var bed = [], h = {}, end = {}, last_u = null, last_r = null, to_end = 0;
var f = new File(arguments[getopt.ind]);
while (f.readline(b) >= 0) {
	var t = b.toString().split("\t");
	var r = t[0] + ":" + (parseInt(t[1]) + 1) + "-" + t[2];
	h[r] = bed.length;
	if (to_end > 0) end[r] = 1, --to_end;
	if (last_u == null || t[3] != last_u) {
		end[r] = 1, to_end = ws - 1;
		if (last_r != null) {
			end[last_r] = 1;
			for (var j = bed.length - 1; j >= 0 && j >= bed.length - ws; --j)
				end[bed[j][2]] = 1;
		}
	}
	var center = Math.floor(parseInt(t[5]) + (parseInt(t[2]) - parseInt(t[1])) / 2);
	bed.push([t[3], t[4], r, center]);
	last_r = r; last_u = t[3];
}
end[last_r] = 1;
for (var j = bed.length - 1; j >= 0 && j >= bed.length - ws; --j)
	end[bed[j][2]] = 1;
f.close();

var paf = [];
f = new File(arguments[getopt.ind+1]);
while (f.readline(b) >= 0) {
	var t = b.toString().split("\t");
	if (parseInt(t[3]) - parseInt(t[2]) < min_span) continue; // a tiny hit
	if (paf.length && t[0] == paf[paf.length - 1][0]) continue; // dup
	var center;
	t[1] = parseInt(t[1]);
	t[2] = parseInt(t[2]); t[3] = parseInt(t[3]);
	t[8] = parseInt(t[8]); t[9] = parseInt(t[9]);
	if (t[4] == '+') {
		center = Math.floor(((t[7] - t[2]) + (t[8] + (t[1] - t[3]))) / 2);
	} else {
		center = Math.floor(((t[7] - (t[1] - t[3])) + (t[8] + t[2])) / 2);
	}
	paf.push([t[0], t[5], t[4], parseInt(t[7]), center]);
}
f.close();

paf.sort(function(x,y){return x[1]<y[1]?-1:x[1]>y[1]?1:x[3]-y[3]});

var chr_se = {}, start = 0;
for (var i = 1; i <= paf.length; ++i) {
	if (i == paf.length || paf[i][1] != paf[i-1][1]) {
		chr_se[paf[i-1][1]] = [start, i];
		start = i;
	}
}

var cnt = 0;
for (var k in chr_se) {
	var st = chr_se[k][0], en = chr_se[k][1];
	for (var i = st + ws + 1; i < en - ws - 1; ++i) {
		var j;
		for (j = i - 1; j >= 0; --j)
			if (paf[i][0] != paf[j][0])
				break;
		if (j < 0) continue; // the first read has multiple mappings
		if (paf[i][1] != paf[j][1]) continue; // different reference chr
		var hi = h[paf[i][0]], hj = h[paf[j][0]];
		var paf_diff = paf[i][4] - paf[j][4];
		var bed_diff = bed[hi][0] == bed[hj][0]? Math.abs(bed[hi][3] - bed[hj][3]) : '*';
		if (hi - hj > ws || hj - hi > ws || bed[hi][0] != bed[hj][0]) {
			if (end[paf[i][0]] != null && end[paf[j][0]] != null) continue;
			if (bed_diff != '*' && Math.abs(paf_diff - bed_diff) < min_span) continue;
			print("E", paf[j][1], bed[hi][0] != bed[hj][0]? '*' : hi-hj, paf_diff, bed_diff, bed[hj][0], bed[hi][0], paf[j][0], paf[i][0]);
			++cnt;
		}
	}
}

print("C", cnt);

b.destroy();
