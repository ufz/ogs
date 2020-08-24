$('.section-docs > h2').add('.section-docs > h3').add('.section-docs > h4').each(function () {
  $('<a class="headerlink"> <i class="far fa-link"></i></a>').
    attr('href', '#' + this.id).
    appendTo(this);
});
