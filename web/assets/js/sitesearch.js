$(document).ready(function () {
  MicroModal.init();

  $('#nav-link-search').click(function (ev) {
    ev.preventDefault();

    MicroModal.show('modal-2', {
      onClose: function () { $('.nav-link-contact').blur(); },
      disableFocus: true
    });

    document.querySelector('.pagefind-ui__search-input').focus();
  });

});
