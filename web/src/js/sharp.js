// const images = require('../../content/features/*.png?{"outputs":["thumbnail"]}'); // one image

// Create thumbs for all images!
var context = require.context("../../content", true, /\.(gif|jpe?g|png|svg|tiff)(\?.*)?$/);
var obj = {};
context.keys().forEach(function (key) {
    obj[key] = context(key);
});
