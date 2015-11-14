var file = arguments.length? new File(arguments[0]) : new File();
var buf = new Bytes();

while (file.readline(buf) >= 0) {
	var x, t = buf.toString().split(/\s+/);
	for (var i = 5; i <= 7; ++i) t[i] = parseInt(t[i]);
	for (var i = 9; i <= 11; ++i) t[i] = parseInt(t[i]);
	var bl = t[6] - t[5] > t[10] - t[9]? t[6] - t[5] : t[10] - t[9];
	var r = parseFloat(t[2]);
	var ml = Math.floor((r <= 1.? bl * r : bl * r / 100) + .499);
	var cm = "cm:i:" + Math.floor(parseFloat(t[3]) + .499);
//	if (t[4] == '1') x = t[7] - t[5], t[5] = t[7] - t[6], t[6] = x;
//	if (t[8] == '1') x = t[11] - t[9], t[9] = t[11] - t[10], t[10] = x;
	var rev = t[4] == t[8]? '+' : '-';
	print([t[0], t[7], t[5], t[6], rev, t[1], t[11], t[9], t[10], ml, bl, 255, cm].join("\t"));
	print([t[1], t[11], t[9], t[10], rev, t[0], t[7], t[5], t[6], ml, bl, 255, cm].join("\t"));
}

buf.destroy();
file.close();
