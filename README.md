# potrace

A TypeScript port of [Potrace](http://potrace.sourceforge.net).<br>
It is based on [kilobtye/potrace](https://github.com/kilobtye/potrace).

## USAGE

### Example 1: Open URL and write SVG

```javascript
var img = new Image();
img.crossOrigin = 'anonymous';
img.src = 'https://www.gravatar.com/avatar/ea4d591101f572e45312cf75901032b4?s=512';
img.onload = function(){
   // Hint: You can also use a canvas as an image.
   document.body.innerHTML = potrace.fromImage(img).toSVG(1); // 1 == scale
}
```

### Example 2: Open custom data and stroke to the Canvas

```javascript
function getPixel(x, y){
   return x % 50 > 25 != y % 50 < 25;
}

var canvas = document.createElement('canvas');
canvas.width = 200;
canvas.height = 200;

var ctx = canvas.getContext('2d');
ctx.beginPath();
potrace.fromFunction(getPixel, canvas.width, canvas.height).strokePath(ctx);
ctx.fill(); // or ctx.stroke();
document.body.appendChild(canvas);
```

### Example 3: Open URL and get Paths

```javascript
var img = new Image();
img.crossOrigin = 'anonymous';
img.src = 'https://www.gravatar.com/avatar/ea4d591101f572e45312cf75901032b4?s=512';
img.onload = function(){
   var o = potrace.fromImage(img).simplify();

   var canvas = document.createElement('canvas');
   canvas.width = o.width;
   canvas.height = o.height;

   var ctx = canvas.getContext('2d');
   ctx.beginPath();
   for (var i = 0; i < o.paths.length; ++i) {
      var path = o.paths[i];
      switch(path.length) {
      case 2:
         ctx.moveTo(path[0], path[1]);
         break;
      case 4:
         ctx.lineTo(path[0], path[1]);
         ctx.lineTo(path[2], path[3]);
         break;
      case 6:
         ctx.bezierCurveTo(
            path[0], path[1],
            path[2], path[3],
            path[4], path[5]);
         break;
      }
   }
   ctx.closePath();
   ctx.fill();
   document.body.appendChild(canvas);
}
```



## LICENSE

GPLv2

This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307, USA.
