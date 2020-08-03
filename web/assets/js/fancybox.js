$(document).ready(function(){
  $('.carousel').slick({
    dots: true,
    infinite: true,
    fade: true,
    cssEase: 'linear',
    autoplaySpeed: 3000,
  });
  $('.carousel-play').slick('slickPlay');

  // wrap all img elements in fancybox wrapper
  var imgs = document.querySelectorAll("img");

  for(var index=0; index < imgs.length; index++) {
    var img = imgs[index];
    // Check for other fancybox
    if(img.parentNode.hasAttribute("data-fancybox")) {
      continue;
    }
    if(img.src.includes("#")) {
      continue;
    }
    var wrapper = document.createElement("a");
    wrapper.setAttribute("data-fancybox", "");
    if(img.hasAttribute("alt")) {
      wrapper.setAttribute("data-caption", img.getAttribute("alt"));
    }
    wrapper.setAttribute("href", img.getAttribute("src"));
    img.parentNode.insertBefore(wrapper, img);
    wrapper.appendChild(img);
  }
});
