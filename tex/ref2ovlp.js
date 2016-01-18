var min_l = 2000, min_q = 10;

var file = arguments.length? new File(arguments[0]) : new File();
var buf = new Bytes();

var a = [];
while (file.readline(buf) >= 0) {
	var t = buf.toString().split("\t");
	for (var j = 1; j <= 3; ++j) t[j] = parseInt(t[j]);
	for (var j = 6; j <= 11; ++j) t[j] = parseInt(t[j]);
	if (t[1] < min_l || t[11] < min_q) continue;
	var st = 0;
	for (var i = 0; i < a.length; ++i) {
		if (t[7] + min_l >= a[i][8]) {
			a[i] = null;
		} else if (t[8] <= a[i][8]) {
			print(t[0], a[i][0], -1);
		} else {
			print(t[0], a[i][0], a[i][8] - t[7]);
		}
	}
	var n = 0;
	for (var i = 0; i < a.length; ++i)
		if (a[i] != null) a[n++] = a[i];
	a.length = n;
	a.push(t);
}

buf.destroy();
file.close();
