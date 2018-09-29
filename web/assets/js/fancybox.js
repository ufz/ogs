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
