$('.section-docs, .section-docs > h2,h3,h4').each(function () {
  $('<a class="headerlink"> <i class="far fa-link"></i></a>').
    attr('href', '#' + this.id).
    appendTo(this);
});
