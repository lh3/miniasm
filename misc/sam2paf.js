var file = arguments.length? new File(arguments[0]) : new File();
var buf = new Bytes();
var re = /(\d+)([MIDSHN])/g;

var len = {}, lineno = 0;
while (file.readline(buf) >= 0) {
	var m, line = buf.toString();
	++lineno;
	if (line.charAt(0) == '@') {
		if (/^@SQ/.test(line)) {
			var name = (m = /\tSN:(\S+)/.exec(line)) != null? m[1] : null;
			var l = (m = /\tLN:(\d+)/.exec(line)) != null? parseInt(m[1]) : null;
			if (name != null && l != null) len[name] = l;
		}
		continue;
	}
	var t = line.split("\t");
	var flag = parseInt(t[1]);
	if (t[2] == '*' || (flag&4)) continue;
	var tlen = len[t[2]];
	if (tlen == null) throw Error("ERROR at line " + lineno + ": can't find the length of contig " + t[2]);
	var NM = (m = /\tNM:i:(\d+)/.exec(line)) != null? parseInt(m[1]) : null;
	if (NM == null) throw Error("ERROR at line " + lineno + ": no NM tag");
	var clip = [0, 0], I = [0, 0], D = [0, 0], M = 0, N = 0;
	while ((m = re.exec(t[5])) != null) {
		var l = parseInt(m[1]);
		if (m[2] == 'M') M += l;
		else if (m[2] == 'I') ++I[0], I[1] += l;
		else if (m[2] == 'D') ++D[0], D[1] += l;
		else if (m[2] == 'N') N += l;
		else if (m[2] == 'S' || m[2] == 'H')
			clip[M == 0? 0 : 1] = l;
	}
	if (NM < I[1] + D[1]) {
		warn("WARNING at line " + lineno + ": NM is less than the total number of gaps; skipped");
		continue;
	}
	var extra = ["mm:i:"+(NM-I[1]-D[1]), "io:i:"+I[0], "in:i:"+I[1], "do:i:"+D[0], "dn:i:"+D[1]];
	var match = M - (NM - I[1] - D[1]);
	var blen = M + I[1] + D[1];
	var qlen = M + I[1] + clip[0] + clip[1];
	var qs, qe;
	if (flag&16) qs = clip[1], qe = qlen - clip[0];
	else qs = clip[0], qe = qlen - clip[1];
	var ts = parseInt(t[3]) - 1, te = ts + M + D[1] + N;
	var a = [t[0], qlen, qs, qe, flag&16? '-' : '+', t[2], tlen, ts, te, match, blen, t[4]];
	print(a.join("\t"), extra.join("\t"));
}

buf.destroy();
file.close();
