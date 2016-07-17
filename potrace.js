/* Copyright (C) 2001-2013 Peter Selinger.
 *
 * A javascript port of Potrace (http://potrace.sourceforge.net).
 *
 * Licensed under the GPL
 *
 * Usage
 *   loadImageFromFile(file) : load image from File API
 *   loadImageFromUrl(url): load image from URL
 *     because of the same-origin policy, can not load image from another domain.
 *     input color/grayscale image is simply converted to binary image. no pre-
 *     process is performed.
 *
 *   setParameter({para1: value, ...}) : set parameters
 *     parameters:
 *        turnpolicy ("black" / "white" / "left" / "right" / "minority" / "majority")
 *          how to resolve ambiguities in path decomposition. (default: "minority")
 *        turdsize
 *          suppress speckles of up to this size (default: 2)
 *        optcurve (true / false)
 *          turn on/off curve optimization (default: true)
 *        alphamax
 *          corner threshold parameter (default: 1)
 *        opttolerance
 *          curve optimization tolerance (default: 0.2)
 *
 *   process(callback) : wait for the image be loaded, then run potrace algorithm,
 *                       then call callback function.
 *
 *   getSVG(size, opt_type) : return a string of generated SVG image.
 *                                    result_image_size = original_image_size * size
 *                                    optional parameter opt_type can be "curve"
 */
var Potrace;
(function (Potrace_1) {
    var Point = (function () {
        function Point(x, y) {
            this.x = x;
            this.y = y;
        }
        Point.prototype.copy = function () {
            return new Point(this.x, this.y);
        };
        return Point;
    }());
    var Bitmap = (function () {
        function Bitmap(width, height) {
            this.width = width;
            this.height = height;
            this.data = new Int8Array(width * height);
        }
        Bitmap.prototype.at = function (x, y) {
            return x >= 0 && x < this.width && y >= 0 && y < this.height &&
                this.data[this.width * y + x] === 1;
        };
        Bitmap.prototype.flip = function (x, y) {
            var i = this.width * y + x;
            this.data[i] = this.data[i] ? 0 : 1;
        };
        Bitmap.prototype.copy = function () {
            var bm = new Bitmap(this.width, this.height);
            for (var i = 0, len = this.data.length; i < len; ++i) {
                bm.data[i] = this.data[i];
            }
            return bm;
        };
        Bitmap.prototype.findNext = function (point) {
            for (var i = this.width * point.y + point.x, len = this.data.length; i < len; ++i) {
                if (this.data[i]) {
                    var y = Math.floor(i / this.width);
                    return new Point(i - y * this.width, y);
                }
            }
            return null;
        };
        Bitmap.createFromImage = function (src) {
            var canvas = document.createElement('canvas');
            canvas.width = src.width;
            canvas.height = src.height;
            canvas.getContext('2d').drawImage(src, 0, 0);
            var bm = new Bitmap(canvas.width, canvas.height);
            var data = canvas.getContext('2d').getImageData(0, 0, bm.width, bm.height).data;
            for (var i = 0, j = 0, l = data.length; i < l; i += 4, ++j) {
                bm.data[j] = 0.2126 * data[i] + 0.7153 * data[i + 1] + 0.0721 * data[i + 2] < 128 ? 1 : 0;
            }
            return bm;
        };
        Bitmap.createFromFunction = function (f, width, height) {
            var bm = new Bitmap(width, height);
            for (var i = 0, y = 0; y < height; ++y) {
                for (var x = 0; x < width; ++i, ++x) {
                    bm.data[i] = f(x, y) ? 1 : 0;
                }
            }
            return bm;
        };
        return Bitmap;
    }());
    var Path = (function () {
        function Path() {
            this.area = 0;
            this.len = 0;
            this.pt = [];
            this.minX = 100000;
            this.minY = 100000;
            this.maxX = -1;
            this.maxY = -1;
            this.sign = '+';
        }
        return Path;
    }());
    var Curve = (function () {
        function Curve(n) {
            this.n = n;
            this.alphaCurve = 0;
            this.tag = new Array(n);
            this.c = new Array(n * 3);
            this.vertex = new Array(n);
            this.alpha = new Array(n);
            this.alpha0 = new Array(n);
            this.beta = new Array(n);
        }
        return Curve;
    }());
    var Quad = (function () {
        function Quad() {
            this.data = [0, 0, 0, 0, 0, 0, 0, 0, 0];
        }
        Quad.prototype.at = function (x, y) {
            return this.data[x * 3 + y];
        };
        return Quad;
    }());
    var Sum = (function () {
        function Sum(x, y, xy, x2, y2) {
            this.x = x;
            this.y = y;
            this.xy = xy;
            this.x2 = x2;
            this.y2 = y2;
        }
        return Sum;
    }());
    var Opti = (function () {
        function Opti() {
            this.pen = 0;
            this.c = [new Point(0, 0), new Point(0, 0)];
            this.t = 0;
            this.s = 0;
            this.alpha = 0;
        }
        return Opti;
    }());
    function majority(bm1, x, y) {
        for (var i = 2; i < 5; i++) {
            var ct = 0;
            for (var a = -i + 1; a <= i - 1; a++) {
                ct += bm1.at(x + a, y + i - 1) ? 1 : -1;
                ct += bm1.at(x + i - 1, y + a - 1) ? 1 : -1;
                ct += bm1.at(x + a - 1, y - i) ? 1 : -1;
                ct += bm1.at(x - i, y + a) ? 1 : -1;
            }
            if (ct > 0) {
                return 1;
            }
            else if (ct < 0) {
                return 0;
            }
        }
        return 0;
    }
    function findPath(bm, infoTurnpolicy, bm1, point) {
        var path = new Path();
        var x = point.x, y = point.y, dirx = 0, diry = 1, tmp;
        path.sign = bm.at(x, y) ? '+' : '-';
        while (1) {
            path.pt.push(new Point(x, y));
            if (x > path.maxX) {
                path.maxX = x;
            }
            if (x < path.minX) {
                path.minX = x;
            }
            if (y > path.maxY) {
                path.maxY = y;
            }
            if (y < path.minY) {
                path.minY = y;
            }
            ++path.len;
            x += dirx;
            y += diry;
            path.area -= x * diry;
            if (x === point.x && y === point.y) {
                break;
            }
            var l = bm1.at(x + (dirx + diry - 1) / 2, y + (diry - dirx - 1) / 2);
            var r = bm1.at(x + (dirx - diry - 1) / 2, y + (diry + dirx - 1) / 2);
            if (r && !l) {
                if (infoTurnpolicy === 'right' ||
                    (infoTurnpolicy === 'black' && path.sign === '+') ||
                    (infoTurnpolicy === 'white' && path.sign === '-') ||
                    (infoTurnpolicy === 'majority' && majority(bm1, x, y)) ||
                    (infoTurnpolicy === 'minority' && !majority(bm1, x, y))) {
                    tmp = dirx;
                    dirx = -diry;
                    diry = tmp;
                }
                else {
                    tmp = dirx;
                    dirx = diry;
                    diry = -tmp;
                }
            }
            else if (r) {
                tmp = dirx;
                dirx = -diry;
                diry = tmp;
            }
            else if (!l) {
                tmp = dirx;
                dirx = diry;
                diry = -tmp;
            }
        }
        return path;
    }
    function xorPath(bm1, path) {
        var y1 = path.pt[0].y;
        var len = path.len;
        for (var i = 1; i < len; i++) {
            var x = path.pt[i].x;
            var y = path.pt[i].y;
            if (y !== y1) {
                var minY = y1 < y ? y1 : y;
                var maxX = path.maxX;
                for (var j = x; j < maxX; j++) {
                    bm1.flip(j, minY);
                }
                y1 = y;
            }
        }
    }
    // --------
    function mod(a, n) {
        return a >= n ? a % n : a >= 0 ? a : n - 1 - (-1 - a) % n;
    }
    function xprod(p1, p2) {
        return p1.x * p2.y - p1.y * p2.x;
    }
    function cyclic(a, b, c) {
        if (a <= c) {
            return a <= b && b < c;
        }
        else {
            return a <= b || b < c;
        }
    }
    function sign(i) {
        return i > 0 ? 1 : i < 0 ? -1 : 0;
    }
    function quadform(Q, w) {
        var v = [w.x, w.y, 1];
        var sum = 0.0;
        for (var i = 0; i < 3; i++) {
            for (var j = 0; j < 3; j++) {
                sum += v[i] * Q.at(i, j) * v[j];
            }
        }
        return sum;
    }
    function interval(lambda, a, b) {
        return new Point(a.x + lambda * (b.x - a.x), a.y + lambda * (b.y - a.y));
    }
    function dorth_infty(p0, p2) {
        return new Point(-sign(p2.y - p0.y), sign(p2.x - p0.x));
    }
    function ddenom(p0, p2) {
        var r = dorth_infty(p0, p2);
        return r.y * (p2.x - p0.x) - r.x * (p2.y - p0.y);
    }
    function dpara(p0, p1, p2) {
        var x1 = p1.x - p0.x;
        var y1 = p1.y - p0.y;
        var x2 = p2.x - p0.x;
        var y2 = p2.y - p0.y;
        return x1 * y2 - x2 * y1;
    }
    function cprod(p0, p1, p2, p3) {
        var x1 = p1.x - p0.x;
        var y1 = p1.y - p0.y;
        var x2 = p3.x - p2.x;
        var y2 = p3.y - p2.y;
        return x1 * y2 - x2 * y1;
    }
    function iprod(p0, p1, p2) {
        var x1 = p1.x - p0.x;
        var y1 = p1.y - p0.y;
        var x2 = p2.x - p0.x;
        var y2 = p2.y - p0.y;
        return x1 * x2 + y1 * y2;
    }
    function iprod1(p0, p1, p2, p3) {
        var x1 = p1.x - p0.x;
        var y1 = p1.y - p0.y;
        var x2 = p3.x - p2.x;
        var y2 = p3.y - p2.y;
        return x1 * x2 + y1 * y2;
    }
    function ddist(p, q) {
        return Math.sqrt((p.x - q.x) * (p.x - q.x) + (p.y - q.y) * (p.y - q.y));
    }
    function bezier(t, p0, p1, p2, p3) {
        var s = 1 - t;
        var s2 = s * s, t2 = t * t;
        var s3 = s2 * s, t3 = t2 * t;
        var s2t3 = 3 * s2 * t, t2s3 = 3 * t2 * s;
        return new Point(s3 * p0.x + s2t3 * p1.x + t2s3 * p2.x + t3 * p3.x, s3 * p0.y + s2t3 * p1.y + t2s3 * p2.y + t3 * p3.y);
    }
    function tangent(p0, p1, p2, p3, q0, q1) {
        var A = cprod(p0, p1, q0, q1);
        var B = cprod(p1, p2, q0, q1);
        var C = cprod(p2, p3, q0, q1);
        var a = A - 2 * B + C;
        var b = -2 * A + 2 * B;
        var c = A;
        var d = b * b - 4 * a * c;
        if (a === 0 || d < 0) {
            return -1.0;
        }
        var s = Math.sqrt(d);
        var r1 = (-b + s) / (2 * a);
        var r2 = (-b - s) / (2 * a);
        if (r1 >= 0 && r1 <= 1) {
            return r1;
        }
        else if (r2 >= 0 && r2 <= 1) {
            return r2;
        }
        else {
            return -1.0;
        }
    }
    function calcSums(path) {
        path.x0 = path.pt[0].x;
        path.y0 = path.pt[0].y;
        path.sums = [];
        var s = path.sums;
        s.push(new Sum(0, 0, 0, 0, 0));
        for (var i = 0; i < path.len; i++) {
            var x = path.pt[i].x - path.x0;
            var y = path.pt[i].y - path.y0;
            s.push(new Sum(s[i].x + x, s[i].y + y, s[i].xy + x * y, s[i].x2 + x * x, s[i].y2 + y * y));
        }
    }
    function calcLon(path) {
        var n = path.len, pt = path.pt;
        var dir, pivk = new Array(n), nc = new Array(n), ct = new Array(4);
        path.lon = new Array(n);
        var constraint = [new Point(0, 0), new Point(0, 0)], cur = new Point(0, 0), off = new Point(0, 0), dk = new Point(0, 0), foundk;
        for (var i = n - 1, k = 0; i >= 0; i--) {
            if (pt[i].x !== pt[k].x && pt[i].y !== pt[k].y) {
                k = i + 1;
            }
            nc[i] = k;
        }
        for (var i = n - 1; i >= 0; i--) {
            ct[0] = ct[1] = ct[2] = ct[3] = 0;
            dir = (3 + 3 * (pt[mod(i + 1, n)].x - pt[i].x) +
                (pt[mod(i + 1, n)].y - pt[i].y)) / 2;
            ct[dir]++;
            constraint[0].x = 0;
            constraint[0].y = 0;
            constraint[1].x = 0;
            constraint[1].y = 0;
            var k = nc[i];
            var k1 = i;
            while (1) {
                foundk = false;
                dir = (3 + 3 * sign(pt[k].x - pt[k1].x) +
                    sign(pt[k].y - pt[k1].y)) / 2;
                ct[dir]++;
                if (ct[0] && ct[1] && ct[2] && ct[3]) {
                    pivk[i] = k1;
                    foundk = true;
                    break;
                }
                cur.x = pt[k].x - pt[i].x;
                cur.y = pt[k].y - pt[i].y;
                if (xprod(constraint[0], cur) < 0 || xprod(constraint[1], cur) > 0) {
                    break;
                }
                if (Math.abs(cur.x) > 1 || Math.abs(cur.y) > 1) {
                    off.x = cur.x + ((cur.y >= 0 && (cur.y > 0 || cur.x < 0)) ? 1 : -1);
                    off.y = cur.y + ((cur.x <= 0 && (cur.x < 0 || cur.y < 0)) ? 1 : -1);
                    if (xprod(constraint[0], off) >= 0) {
                        constraint[0].x = off.x;
                        constraint[0].y = off.y;
                    }
                    off.x = cur.x + ((cur.y <= 0 && (cur.y < 0 || cur.x < 0)) ? 1 : -1);
                    off.y = cur.y + ((cur.x >= 0 && (cur.x > 0 || cur.y < 0)) ? 1 : -1);
                    if (xprod(constraint[1], off) <= 0) {
                        constraint[1].x = off.x;
                        constraint[1].y = off.y;
                    }
                }
                k1 = k;
                k = nc[k1];
                if (!cyclic(k, i, k1)) {
                    break;
                }
            }
            if (!foundk) {
                dk.x = sign(pt[k].x - pt[k1].x);
                dk.y = sign(pt[k].y - pt[k1].y);
                cur.x = pt[k1].x - pt[i].x;
                cur.y = pt[k1].y - pt[i].y;
                var a = xprod(constraint[0], cur);
                var b = xprod(constraint[0], dk);
                var c = xprod(constraint[1], cur);
                var d = xprod(constraint[1], dk);
                var j_1 = 10000000;
                if (b < 0) {
                    j_1 = Math.floor(a / -b);
                }
                if (d > 0) {
                    j_1 = Math.min(j_1, Math.floor(-c / d));
                }
                pivk[i] = mod(k1 + j_1, n);
            }
        }
        var j = pivk[n - 1];
        path.lon[n - 1] = j;
        for (var i = n - 2; i >= 0; i--) {
            if (cyclic(i + 1, pivk[i], j)) {
                j = pivk[i];
            }
            path.lon[i] = j;
        }
        for (var i = n - 1; cyclic(mod(i + 1, n), j, path.lon[i]); i--) {
            path.lon[i] = j;
        }
    }
    function penalty3(path, i, j) {
        var n = path.len, pt = path.pt, sums = path.sums;
        var r = 0;
        if (j >= n) {
            j -= n;
            r = 1;
        }
        var x, y, x2, xy, y2, k;
        if (r === 0) {
            x = sums[j + 1].x - sums[i].x;
            y = sums[j + 1].y - sums[i].y;
            x2 = sums[j + 1].x2 - sums[i].x2;
            xy = sums[j + 1].xy - sums[i].xy;
            y2 = sums[j + 1].y2 - sums[i].y2;
            k = j + 1 - i;
        }
        else {
            x = sums[j + 1].x - sums[i].x + sums[n].x;
            y = sums[j + 1].y - sums[i].y + sums[n].y;
            x2 = sums[j + 1].x2 - sums[i].x2 + sums[n].x2;
            xy = sums[j + 1].xy - sums[i].xy + sums[n].xy;
            y2 = sums[j + 1].y2 - sums[i].y2 + sums[n].y2;
            k = j + 1 - i + n;
        }
        var px = (pt[i].x + pt[j].x) / 2.0 - pt[0].x;
        var py = (pt[i].y + pt[j].y) / 2.0 - pt[0].y;
        var ey = (pt[j].x - pt[i].x);
        var ex = -(pt[j].y - pt[i].y);
        var a = ((x2 - 2 * x * px) / k + px * px);
        var b = ((xy - x * py - y * px) / k + px * py);
        var c = ((y2 - 2 * y * py) / k + py * py);
        var s = ex * ex * a + 2 * ex * ey * b + ey * ey * c;
        return Math.sqrt(s);
    }
    function bestPolygon(path) {
        var n = path.len;
        var pen = new Array(n + 1);
        var prev = new Array(n + 1);
        var clip0 = new Array(n);
        var clip1 = new Array(n + 1);
        var seg0 = new Array(n + 1);
        var seg1 = new Array(n + 1);
        for (var i_1 = 0; i_1 < n; i_1++) {
            var c = mod(path.lon[mod(i_1 - 1, n)] - 1, n);
            if (c === i_1) {
                c = mod(i_1 + 1, n);
            }
            if (c < i_1) {
                clip0[i_1] = n;
            }
            else {
                clip0[i_1] = c;
            }
        }
        for (var i_2 = 0, j_2 = 1; i_2 < n; i_2++) {
            while (j_2 <= clip0[i_2]) {
                clip1[j_2] = i_2;
                j_2++;
            }
        }
        var j = 0;
        for (var i_3 = 0; i_3 < n; j++) {
            seg0[j] = i_3;
            i_3 = clip0[i_3];
        }
        seg0[j] = n;
        var m = j;
        var i = n;
        for (var j_3 = m; j_3 > 0; j_3--) {
            seg1[j_3] = i;
            i = clip1[i];
        }
        seg1[0] = 0;
        pen[0] = 0;
        for (var j_4 = 1; j_4 <= m; j_4++) {
            for (var i_4 = seg1[j_4]; i_4 <= seg0[j_4]; i_4++) {
                var best = -1;
                for (var k = seg0[j_4 - 1]; k >= clip1[i_4]; k--) {
                    var thispen = penalty3(path, k, i_4) + pen[k];
                    if (best < 0 || thispen < best) {
                        prev[i_4] = k;
                        best = thispen;
                    }
                }
                pen[i_4] = best;
            }
        }
        path.m = m;
        path.po = new Array(m);
        for (var i_5 = n, j_5 = m - 1; i_5 > 0; j_5--) {
            i_5 = prev[i_5];
            path.po[j_5] = i_5;
        }
    }
    function pointslope(path, i, j, ctr, dir) {
        var n = path.len, sums = path.sums;
        var r = 0;
        while (j >= n) {
            j -= n;
            r += 1;
        }
        while (i >= n) {
            i -= n;
            r -= 1;
        }
        while (j < 0) {
            j += n;
            r -= 1;
        }
        while (i < 0) {
            i += n;
            r += 1;
        }
        var x = sums[j + 1].x - sums[i].x + r * sums[n].x;
        var y = sums[j + 1].y - sums[i].y + r * sums[n].y;
        var x2 = sums[j + 1].x2 - sums[i].x2 + r * sums[n].x2;
        var xy = sums[j + 1].xy - sums[i].xy + r * sums[n].xy;
        var y2 = sums[j + 1].y2 - sums[i].y2 + r * sums[n].y2;
        var k = j + 1 - i + r * n;
        ctr.x = x / k;
        ctr.y = y / k;
        var a = (x2 - x * x / k) / k;
        var b = (xy - x * y / k) / k;
        var c = (y2 - y * y / k) / k;
        var lambda2 = (a + c + Math.sqrt((a - c) * (a - c) + 4 * b * b)) / 2;
        a -= lambda2;
        c -= lambda2;
        var l;
        if (Math.abs(a) >= Math.abs(c)) {
            l = Math.sqrt(a * a + b * b);
            if (l !== 0) {
                dir.x = -b / l;
                dir.y = a / l;
            }
        }
        else {
            l = Math.sqrt(c * c + b * b);
            if (l !== 0) {
                dir.x = -c / l;
                dir.y = b / l;
            }
        }
        if (l === 0) {
            dir.x = dir.y = 0;
        }
    }
    function adjustVertices(path) {
        var m = path.m, po = path.po, n = path.len, pt = path.pt, x0 = path.x0, y0 = path.y0;
        var ctr = new Array(m), dir = new Array(m), q = new Array(m);
        var v = new Array(3);
        var s = new Point(0, 0);
        path.curve = new Curve(m);
        for (var i = 0; i < m; ++i) {
            var j = po[mod(i + 1, m)];
            j = mod(j - po[i], n) + po[i];
            ctr[i] = new Point(0, 0);
            dir[i] = new Point(0, 0);
            pointslope(path, po[i], j, ctr[i], dir[i]);
        }
        for (var i = 0; i < m; ++i) {
            q[i] = new Quad();
            var d = dir[i].x * dir[i].x + dir[i].y * dir[i].y;
            if (d === 0.0) {
                for (var j = 0; j < 3; ++j) {
                    for (var k = 0; k < 3; ++k) {
                        q[i].data[j * 3 + k] = 0;
                    }
                }
            }
            else {
                v[0] = dir[i].y;
                v[1] = -dir[i].x;
                v[2] = -v[1] * ctr[i].y - v[0] * ctr[i].x;
                for (var l = 0; l < 3; ++l) {
                    for (var k = 0; k < 3; ++k) {
                        q[i].data[l * 3 + k] = v[l] * v[k] / d;
                    }
                }
            }
        }
        for (var i = 0; i < m; ++i) {
            var Q = new Quad();
            var w = new Point(0, 0);
            s.x = pt[po[i]].x - x0;
            s.y = pt[po[i]].y - y0;
            var j = mod(i - 1, m);
            for (var l = 0; l < 3; ++l) {
                for (var k = 0; k < 3; ++k) {
                    Q.data[l * 3 + k] = q[j].at(l, k) + q[i].at(l, k);
                }
            }
            while (true) {
                var det = Q.at(0, 0) * Q.at(1, 1) - Q.at(0, 1) * Q.at(1, 0);
                if (det !== 0.0) {
                    w.x = (-Q.at(0, 2) * Q.at(1, 1) + Q.at(1, 2) * Q.at(0, 1)) / det;
                    w.y = (Q.at(0, 2) * Q.at(1, 0) - Q.at(1, 2) * Q.at(0, 0)) / det;
                    break;
                }
                if (Q.at(0, 0) > Q.at(1, 1)) {
                    v[0] = -Q.at(0, 1);
                    v[1] = Q.at(0, 0);
                }
                else if (Q.at(1, 1)) {
                    v[0] = -Q.at(1, 1);
                    v[1] = Q.at(1, 0);
                }
                else {
                    v[0] = 1;
                    v[1] = 0;
                }
                var d = v[0] * v[0] + v[1] * v[1];
                v[2] = -v[1] * s.y - v[0] * s.x;
                for (var l = 0; l < 3; ++l) {
                    for (var k = 0; k < 3; ++k) {
                        Q.data[l * 3 + k] += v[l] * v[k] / d;
                    }
                }
            }
            var dx = Math.abs(w.x - s.x);
            var dy = Math.abs(w.y - s.y);
            if (dx <= 0.5 && dy <= 0.5) {
                path.curve.vertex[i] = new Point(w.x + x0, w.y + y0);
                continue;
            }
            var min = quadform(Q, s);
            var xmin = s.x;
            var ymin = s.y;
            if (Q.at(0, 0) !== 0.0) {
                for (var z = 0; z < 2; ++z) {
                    w.y = s.y - 0.5 + z;
                    w.x = -(Q.at(0, 1) * w.y + Q.at(0, 2)) / Q.at(0, 0);
                    dx = Math.abs(w.x - s.x);
                    var cand = quadform(Q, w);
                    if (dx <= 0.5 && cand < min) {
                        min = cand;
                        xmin = w.x;
                        ymin = w.y;
                    }
                }
            }
            if (Q.at(1, 1) !== 0.0) {
                for (var z = 0; z < 2; ++z) {
                    w.x = s.x - 0.5 + z;
                    w.y = -(Q.at(1, 0) * w.x + Q.at(1, 2)) / Q.at(1, 1);
                    dy = Math.abs(w.y - s.y);
                    var cand = quadform(Q, w);
                    if (dy <= 0.5 && cand < min) {
                        min = cand;
                        xmin = w.x;
                        ymin = w.y;
                    }
                }
            }
            for (var l = 0; l < 2; ++l) {
                for (var k = 0; k < 2; ++k) {
                    w.x = s.x - 0.5 + l;
                    w.y = s.y - 0.5 + k;
                    var cand = quadform(Q, w);
                    if (cand < min) {
                        min = cand;
                        xmin = w.x;
                        ymin = w.y;
                    }
                }
            }
            path.curve.vertex[i] = new Point(xmin + x0, ymin + y0);
        }
    }
    function reverse(path) {
        var curve = path.curve, m = curve.n, v = curve.vertex;
        for (var i = 0, j = m - 1; i < j; ++i, --j) {
            var tmp = v[i];
            v[i] = v[j];
            v[j] = tmp;
        }
    }
    function smooth(path, infoAlphamax) {
        var m = path.curve.n, curve = path.curve;
        for (var i = 0; i < m; ++i) {
            var j = mod(i + 1, m);
            var k = mod(i + 2, m);
            var p4 = interval(1 / 2.0, curve.vertex[k], curve.vertex[j]);
            var denom = ddenom(curve.vertex[i], curve.vertex[k]);
            var alpha = void 0;
            if (denom !== 0.0) {
                var dd = Math.abs(dpara(curve.vertex[i], curve.vertex[j], curve.vertex[k]) / denom);
                alpha = (dd > 1 ? (1 - 1.0 / dd) : 0) / 0.75;
            }
            else {
                alpha = 4 / 3.0;
            }
            curve.alpha0[j] = alpha;
            if (alpha >= infoAlphamax) {
                curve.tag[j] = 'CORNER';
                curve.c[3 * j + 1] = curve.vertex[j];
                curve.c[3 * j + 2] = p4;
            }
            else {
                if (alpha < 0.55) {
                    alpha = 0.55;
                }
                else if (alpha > 1) {
                    alpha = 1;
                }
                curve.tag[j] = 'CURVE';
                curve.c[3 * j + 0] = interval(0.5 + 0.5 * alpha, curve.vertex[i], curve.vertex[j]);
                curve.c[3 * j + 1] = interval(0.5 + 0.5 * alpha, curve.vertex[k], curve.vertex[j]);
                curve.c[3 * j + 2] = p4;
            }
            curve.alpha[j] = alpha;
            curve.beta[j] = 0.5;
        }
        curve.alphaCurve = 1;
    }
    function opti_penalty(path, i, j, res, opttolerance, convc, areac) {
        var m = path.curve.n, curve = path.curve, vertex = curve.vertex;
        if (i === j) {
            return true;
        }
        var k = i;
        var i1 = mod(i + 1, m);
        var k1 = mod(k + 1, m);
        var conv = convc[k1];
        if (conv === 0) {
            return true;
        }
        var d = ddist(vertex[i], vertex[i1]);
        for (var k_1 = k1; k_1 !== j; k_1 = k1) {
            k1 = mod(k_1 + 1, m);
            var k2 = mod(k_1 + 2, m);
            if (convc[k1] !== conv) {
                return true;
            }
            if (sign(cprod(vertex[i], vertex[i1], vertex[k1], vertex[k2])) !== conv) {
                return true;
            }
            if (iprod1(vertex[i], vertex[i1], vertex[k1], vertex[k2]) <
                d * ddist(vertex[k1], vertex[k2]) * -0.999847695156) {
                return true;
            }
        }
        var p0 = curve.c[mod(i, m) * 3 + 2].copy();
        var p1 = vertex[mod(i + 1, m)].copy();
        var p2 = vertex[mod(j, m)].copy();
        var p3 = curve.c[mod(j, m) * 3 + 2].copy();
        var area = areac[j] - areac[i];
        area -= dpara(vertex[0], curve.c[i * 3 + 2], curve.c[j * 3 + 2]) / 2;
        if (i >= j) {
            area += areac[m];
        }
        var A1 = dpara(p0, p1, p2);
        var A2 = dpara(p0, p1, p3);
        var A3 = dpara(p0, p2, p3);
        if (A2 === A1) {
            return true;
        }
        var A4 = A1 + A3 - A2;
        var t = A3 / (A3 - A4);
        var s = A2 / (A2 - A1);
        var A = A2 * t / 2.0;
        if (A === 0.0) {
            return true;
        }
        var R = area / A;
        var alpha = 2 - Math.sqrt(4 - R / 0.3);
        res.c[0] = interval(t * alpha, p0, p1);
        res.c[1] = interval(s * alpha, p3, p2);
        res.alpha = alpha;
        res.t = t;
        res.s = s;
        p1 = res.c[0].copy();
        p2 = res.c[1].copy();
        res.pen = 0;
        for (var k_2 = mod(i + 1, m), k1_1; k_2 !== j; k_2 = k1_1) {
            k1_1 = mod(k_2 + 1, m);
            var t_1 = tangent(p0, p1, p2, p3, vertex[k_2], vertex[k1_1]);
            if (t_1 < -0.5) {
                return true;
            }
            var pt = bezier(t_1, p0, p1, p2, p3);
            var d_1 = ddist(vertex[k_2], vertex[k1_1]);
            if (d_1 === 0.0) {
                return true;
            }
            var d1 = dpara(vertex[k_2], vertex[k1_1], pt) / d_1;
            if (Math.abs(d1) > opttolerance) {
                return true;
            }
            if (iprod(vertex[k_2], vertex[k1_1], pt) < 0 ||
                iprod(vertex[k1_1], vertex[k_2], pt) < 0) {
                return true;
            }
            res.pen += d1 * d1;
        }
        for (var k_3 = i, k1_2; k_3 !== j; k_3 = k1_2) {
            k1_2 = mod(k_3 + 1, m);
            var t_2 = tangent(p0, p1, p2, p3, curve.c[k_3 * 3 + 2], curve.c[k1_2 * 3 + 2]);
            if (t_2 < -0.5) {
                return true;
            }
            var pt = bezier(t_2, p0, p1, p2, p3);
            var d_2 = ddist(curve.c[k_3 * 3 + 2], curve.c[k1_2 * 3 + 2]);
            if (d_2 === 0.0) {
                return true;
            }
            var d1 = dpara(curve.c[k_3 * 3 + 2], curve.c[k1_2 * 3 + 2], pt) / d_2;
            var d2 = dpara(curve.c[k_3 * 3 + 2], curve.c[k1_2 * 3 + 2], vertex[k1_2]) / d_2;
            d2 *= 0.75 * curve.alpha[k1_2];
            if (d2 < 0) {
                d1 = -d1;
                d2 = -d2;
            }
            if (d1 < d2 - opttolerance) {
                return true;
            }
            if (d1 < d2) {
                res.pen += (d1 - d2) * (d1 - d2);
            }
        }
        return false;
    }
    function optiCurve(path, infoOpttolerance) {
        var curve = path.curve, m = curve.n, vert = curve.vertex, pt = new Array(m + 1), pen = new Array(m + 1), len = new Array(m + 1), opt = new Array(m + 1);
        var convc = new Array(m), areac = new Array(m + 1);
        for (var i = 0; i < m; ++i) {
            if (curve.tag[i] === 'CURVE') {
                convc[i] = sign(dpara(vert[mod(i - 1, m)], vert[i], vert[mod(i + 1, m)]));
            }
            else {
                convc[i] = 0;
            }
        }
        var area = 0.0;
        areac[0] = 0.0;
        var p0 = curve.vertex[0];
        for (var i = 0; i < m; ++i) {
            var i1 = mod(i + 1, m);
            if (curve.tag[i1] === 'CURVE') {
                var alpha = curve.alpha[i1];
                area += 0.3 * alpha * (4 - alpha) *
                    dpara(curve.c[i * 3 + 2], vert[i1], curve.c[i1 * 3 + 2]) / 2;
                area += dpara(p0, curve.c[i * 3 + 2], curve.c[i1 * 3 + 2]) / 2;
            }
            areac[i + 1] = area;
        }
        pt[0] = -1;
        pen[0] = 0;
        len[0] = 0;
        var o = new Opti();
        for (var j = 1; j <= m; ++j) {
            pt[j] = j - 1;
            pen[j] = pen[j - 1];
            len[j] = len[j - 1] + 1;
            for (var i = j - 2; i >= 0; --i) {
                var r = opti_penalty(path, i, mod(j, m), o, infoOpttolerance, convc, areac);
                if (r) {
                    break;
                }
                if (len[j] > len[i] + 1 ||
                    (len[j] === len[i] + 1 && pen[j] > pen[i] + o.pen)) {
                    pt[j] = i;
                    pen[j] = pen[i] + o.pen;
                    len[j] = len[i] + 1;
                    opt[j] = o;
                    o = new Opti();
                }
            }
        }
        var om = len[m];
        var ocurve = new Curve(om);
        var s = new Array(om);
        var t = new Array(om);
        for (var i = om - 1, j = m; i >= 0; --i) {
            if (pt[j] === j - 1) {
                ocurve.tag[i] = curve.tag[mod(j, m)];
                ocurve.c[i * 3 + 0] = curve.c[mod(j, m) * 3 + 0];
                ocurve.c[i * 3 + 1] = curve.c[mod(j, m) * 3 + 1];
                ocurve.c[i * 3 + 2] = curve.c[mod(j, m) * 3 + 2];
                ocurve.vertex[i] = curve.vertex[mod(j, m)];
                ocurve.alpha[i] = curve.alpha[mod(j, m)];
                ocurve.alpha0[i] = curve.alpha0[mod(j, m)];
                ocurve.beta[i] = curve.beta[mod(j, m)];
                s[i] = t[i] = 1.0;
            }
            else {
                ocurve.tag[i] = 'CURVE';
                ocurve.c[i * 3 + 0] = opt[j].c[0];
                ocurve.c[i * 3 + 1] = opt[j].c[1];
                ocurve.c[i * 3 + 2] = curve.c[mod(j, m) * 3 + 2];
                ocurve.vertex[i] = interval(opt[j].s, curve.c[mod(j, m) * 3 + 2], vert[mod(j, m)]);
                ocurve.alpha[i] = opt[j].alpha;
                ocurve.alpha0[i] = opt[j].alpha;
                s[i] = opt[j].s;
                t[i] = opt[j].t;
            }
            j = pt[j];
        }
        for (var i = 0; i < om; ++i) {
            ocurve.beta[i] = s[i] / (s[i] + t[mod(i + 1, om)]);
        }
        ocurve.alphaCurve = 1;
        path.curve = ocurve;
    }
    // --------
    function convertSVG(width, height, pathlist, scale, opt_type) {
        var w = width * scale, h = height * scale;
        var svg = [("<svg id=\"svg\" version=\"1.1\" width=\"" + w + "\" height=\"" + h + "\" xmlns=\"http://www.w3.org/2000/svg\">")];
        svg.push('<path d="');
        for (var i = 0, len = pathlist.length; i < len; ++i) {
            var curve = pathlist[i].curve, n = curve.n;
            svg.push('M' + (curve.c[(n - 1) * 3 + 2].x * scale).toFixed(3) +
                ' ' + (curve.c[(n - 1) * 3 + 2].y * scale).toFixed(3) + ' ');
            for (var i_6 = 0; i_6 < n; ++i_6) {
                if (curve.tag[i_6] === 'CURVE') {
                    svg.push('C ' + (curve.c[i_6 * 3 + 0].x * scale).toFixed(3) + ' ' +
                        (curve.c[i_6 * 3 + 0].y * scale).toFixed(3) + ',');
                    svg.push((curve.c[i_6 * 3 + 1].x * scale).toFixed(3) + ' ' +
                        (curve.c[i_6 * 3 + 1].y * scale).toFixed(3) + ',');
                    svg.push((curve.c[i_6 * 3 + 2].x * scale).toFixed(3) + ' ' +
                        (curve.c[i_6 * 3 + 2].y * scale).toFixed(3) + ' ');
                }
                else if (curve.tag[i_6] === 'CORNER') {
                    svg.push('L ' + (curve.c[i_6 * 3 + 1].x * scale).toFixed(3) + ' ' +
                        (curve.c[i_6 * 3 + 1].y * scale).toFixed(3) + ' ');
                    svg.push((curve.c[i_6 * 3 + 2].x * scale).toFixed(3) + ' ' +
                        (curve.c[i_6 * 3 + 2].y * scale).toFixed(3) + ' ');
                }
            }
        }
        if (opt_type === 'curve') {
            svg.push('" stroke="black" fill="none"/>');
        }
        else {
            svg.push('" stroke="none" fill="black" fill-rule="evenodd"/>');
        }
        svg.push('</svg>');
        return svg.join('');
    }
    function fromImage(src) {
        return new Potrace(Bitmap.createFromImage(src));
    }
    Potrace_1.fromImage = fromImage;
    function fromFunction(f, width, height) {
        return new Potrace(Bitmap.createFromFunction(f, width, height));
    }
    Potrace_1.fromFunction = fromFunction;
    var Potrace = (function () {
        function Potrace(bm) {
            this.pathlist = [];
            this.turnPolicy = 'minority';
            this.turdSize = 2;
            this.optCurve = true;
            this.alphaMax = 1;
            this.optTolerance = 0.2;
            this.width = bm.width;
            this.height = bm.height;
            // bitmap to pathlist
            var pathlist = this.pathlist;
            var bm1 = bm.copy();
            var currentPoint = new Point(0, 0);
            while (currentPoint = bm1.findNext(currentPoint)) {
                var path = findPath(bm, this.turnPolicy, bm1, currentPoint);
                xorPath(bm1, path);
                if (path.area > this.turdSize) {
                    pathlist.push(path);
                }
            }
            // process path
            for (var i = 0; i < pathlist.length; ++i) {
                var path = pathlist[i];
                calcSums(path);
                calcLon(path);
                bestPolygon(path);
                adjustVertices(path);
                if (path.sign === '-') {
                    reverse(path);
                }
                smooth(path, this.alphaMax);
                if (this.optCurve) {
                    optiCurve(path, this.optTolerance);
                }
            }
        }
        Potrace.prototype.getSVG = function (scale, optType) {
            return convertSVG(this.width, this.height, this.pathlist, scale, optType);
        };
        return Potrace;
    }());
})(Potrace || (Potrace = {}));
