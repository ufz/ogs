var $ = require("jquery");
var slick = require("slick-carousel");
require('style-loader!slick-carousel/slick/slick.css');
// Do not require here, import and override in styles.css
// require('style-loader!slick-carousel/slick/slick-theme.css');
require("@fancyapps/fancybox");
require('style-loader!@fancyapps/fancybox/dist/jquery.fancybox.min.css');

$(document).ready(function(){
  $('.carousel').slick({
    dots: true,
    infinite: true,
    fade: true,
    cssEase: 'linear',
    autoplaySpeed: 3000,
  });
  $('.carousel-play').slick('slickPlay');
});
