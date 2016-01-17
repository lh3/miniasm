if (arguments.length != 2) {
	print("Usage: k8 paf_srtcmp.js <bwamem.srt.paf> <minimap.srt.paf>");
	exit(1);
}

function read1(f, buf, last)
{
	var a = [], l = null;
	if (last != null) a.push(last);
	while (f.readline(buf) >= 0) {
		var t = buf.toString().split("\t");
		for (var j = 1; j <= 3; ++j) t[j] = parseInt(t[j]);
		for (var j = 6; j <= 11; ++j) t[j] = parseInt(t[j]);
		if (last == null) {
			last = t;
			a.push(t);
		} else if (last[0] != t[0]) {
			l = t;
			break;
		} else a.push(t);
	}
//	if (a.length > 0) print(a[0][0], a.length);
	return [a, l];
}

var buf = new Bytes();
var fb = new File(arguments[0]);
var fm = new File(arguments[1]);

var tot = 0, matched = 0;

var sb = read1(fb, buf, null);
var sm = read1(fm, buf, null);

function sync()
{
//	print("here!");
	while (sb[0][0][0] != sm[0][0][0]) {
		if (sb[0][0][0] < sm[0][0][0]) {
			if (sb[0].length == 1) ++tot;
			sb = read1(fb, buf, sb[1]);
			if (sb[0].length == 0) break;
		} else if (sb[0][0][0] > sm[0][0][0]) {
			sm = read1(fm, buf, sm[1]);
			if (sm[0].length == 0) break;
		}
	}
}

while (1) {
	sync();
	if (sb[0].length == 0) break;
	if (sm[0].length == 0) {
		while (sb[0].length) {
			if (sb[0].length == 1) ++tot;
			sb = read1(fb, buf, sb[1]);
		}
		break;
	}
	if (sb[0].length == 1) {
		var end = sm[0].length, hit = 0;
		++tot;
		for (var j = 0; j < end; ++j) {
			if (sb[0][0][4] != sm[0][j][4] || sb[0][0][5] != sm[0][j][5]) continue;
			if (sb[0][0][8] > sm[0][j][7] && sm[0][j][8] > sb[0][0][7]) {
				var ol, ml;
				ol = sb[0][0][8] - sm[0][j][7];
				ml = sm[0][j][8] - sb[0][0][7];
				var r = ol < ml? ol / ml : ml / ol;
				if (r >= .3333) ++matched, hit = 1;
				break;
			}
		}
		if (hit == 0) print(sb[0][0].join("\t"));
	}
	sb = read1(fb, buf, sb[1]);
	sm = read1(fm, buf, sm[1]);
	if (sb[0].length == 0) break;
}

fb.close();
fm.close();
buf.destroy();

print(tot, matched, matched/tot);
